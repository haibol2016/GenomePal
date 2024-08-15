#' Create a BSgenome package
#'
#' Starting with a multi-fasta file for a reference genome, a BSgenome package
#' is created and installed.
#'
#' @param genome_fasta A character(1) specifying a path for a multi-fasta file
#'   for a reference genome. It can be a gzip-compressed or uncompressed fasta
#'   file.
#' @param dest_dir A character(1) specifying a path for output split, compressed
#'   fasta files
#' @param prefix A character(1), specifying chromosome prefix. It can be an
#'   empty string or "chr". Default is "".
#' @param latin_name A character(1), specifying the Latin name of the species of
#'   the reference genome, such as "Homo sapiens" for the human.
#' @param common_name A character(1), specifying the common name of the species
#'   of the reference genome, such as "Human".
#' @param genome_version A character(1), specifying the genome build, such as
#'   "GRCh38".
#' @param seed_file_name A character(1), specifying the path to a seed file to
#'   be created for building a BSgenome package.
#' @param fasta_url A character(1), specifying the URL from where the reference
#'   genome fasta file is downloaded.
#' @param release_date  A character(1), specifying the genome assembly release
#'   date, such as "August 2020".
#' @param source A character(1), specifying the source of the genome fasta
#'   file, such as "Ensembl".
#' @param version A character(1), specifying the version of a BSgenome package
#'   to be build, such as "1.0.0".
#'
#' @return A BSgenome package name. Users have to install the package
#'   to use it.
#' @export
#'
#' @importFrom BSgenomeForge forgeBSgenomeDataPkg
#'
#' @examples
#' if (TRUE) {
#'     dest_dir <- tempdir()
#'     genome_fasta <- file.path(dest_dir, "toy.example.fa")
#'     in_fasta <- file(genome_fasta, open = "w")
#'     writeLines(c(">chr1\nATCGCTGCGGATGCGG",
#'               ">chr2\nCCCCCCCCCCCGGAAA",
#'               ">chrM\nATACTACTGGA"), in_fasta)
#'     close(in_fasta)
#'
#'     make_BSgenome_pkg(
#'         genome_fasta = genome_fasta,
#'         dest_dir = file.path(dest_dir, "multifasta"),
#'         prefix = "none",
#'         latin_name = "Homo sapiens",
#'         common_name = "Human",
#'         genome_version = "GRCh38",
#'         seed_file_name = file.path(
#'             dest_dir,
#'             "human.genome.seed.txt"
#'         ),
#'         fasta_url = "http://ftp.ensembl.org/xxx.fa",
#'         release_date = "August 2007",
#'         source = "Ensembl",
#'         version = "1.0.0"
#'     )
#' }
make_BSgenome_pkg <-
    function(
        genome_fasta = NULL,
        dest_dir = tempdir(),
        prefix = c("", "chr"),
        latin_name = NULL,
        common_name = NULL,
        genome_version = NULL,
        seed_file_name = NULL,
        fasta_url = NULL,
        release_date = "August 2007",
        source = "Ensembl",
        version = "1.0.0") {
        if (!dir.exists(dest_dir)) {
            dir.create(dest_dir, recursive = TRUE)
        }

        multifasta_path <-
            generate_multifasta(
                genome_fasta = genome_fasta,
                dest_dir = dest_dir,
                prefix = prefix
            )

        seed_file <- generate_seed_file(
            multifasta_path = multifasta_path,
            latin_name = latin_name,
            common_name = common_name,
            genome_version = genome_version,
            seed_file_name = seed_file_name,
            fasta_url = fasta_url,
            release_date = release_date,
            source = source,
            version = version
        )

        pkgname <- forge_BSgenome(
            seed_file = seed_file[2],
            dest_dir = dest_dir
        )
        pkgname
    }


#' Split a mutli-fasta file by chromosome/scaffold into individual fasta file
#'
#' Given a multi-fasta file for a reference genome, create a gzip-compressed
#' fasta file for each chromosome/scaffold, such as chr1.fa.gz, chr2.fa.gz,...
#' A prefix "chr" may be added if necessary.
#'
#' @param genome_fasta A character(1) specifying a path for a multi-fasta file
#'   for a reference genome. It can be a gzip-compressed or uncompressed fasta
#'   file.
#' @param dest_dir A character(1) specifying a path for output split, compressed
#'   fasta files
#' @param prefix A character(1), specifying chromosome prefix. It can be "none"
#'   or "chr". Default is "none".
#'
#' @return A character string, which is th path to the directory where split,
#'   compressed fasta files are located.
#' @noRd
#'
#' @examples
#' dest_dir <- tempdir()
#' genome_fasta <- file.path(dest_dir, "toy.example.fa")
#' in_fasta <- file(genome_fasta, open = "w")
#' writeLines(c(">chr1\nATCGCTGCGGATGCGG",
#'               ">chr2\nCCCCCCCCCCCGGAAA",
#'               ">chrM\nATACTACTGGA"), in_fasta)
#' close(in_fasta)
#'
#' fa_dir <- generate_multifasta(
#'     genome_fasta = genome_fasta,
#'     dest_dir = dest_dir,
#'     prefix = "none"
#' )
#'
generate_multifasta <- function(genome_fasta = NULL,
                                dest_dir = NULL,
                                prefix = c("none", "chr")) {
    if (!file.exists(genome_fasta)) {
        stop(
            "A path to a reference genome fasta file, genome_fasta,",
            " are required, but it doesn't exist"
        )
    }

    if (is.null(dest_dir)) {
        stop("dest_dir is required")
    }

    if (!dir.exists(dest_dir)) {
        dir.create(dest_dir, recursive = TRUE)
    }

    prefix <-
        match.arg(prefix,
                  choices = c("none", "chr"),
                  several.ok = FALSE
        )
    if (prefix == "none") {
        prefix <- ""
    }

    if (grepl(".(fa|fasta).gz$", genome_fasta)) {
        in_fasta <- gzfile(genome_fasta, open = "rt")
    } else if (grepl(".(fa|fasta)$", genome_fasta)) {
        in_fasta <- file(genome_fasta, open = "r")
    } else {
        stop(
            "It seems the genome sequence file is not a valid fasta file ",
            "which should with an extension .fa, .fasta, .fa.gz, or .fasta.gz"
        )
    }

    f <- ""
    line_num <- 1
    while (length({
        line <- readLines(in_fasta, n = 1, warn = FALSE)
    }) > 0) {
        line <- trimws(line) # remove leading and trailing white spaces
        if (length(line) == 0) {
            next
        }
        if (grepl("^>", line)) {
            if (line_num > 1) {
                close(f)
            }
            f <- gzfile(
                file.path(
                    dest_dir,
                    gsub(
                        "^>(chr)?([^\\s]+).*",
                        paste0(prefix, "\\2.fa.gz"),
                        line,
                        perl = TRUE
                    )
                ),
                "w"
            )
            writeLines(gsub(
                "^>(chr)?([^\\s]+).*",
                paste0(">", prefix, "\\2"),
                line,
                perl = TRUE
            ), f)
        } else {
            writeLines(line, f)
        }
        line_num <- line_num + 1
    }
    close(f)
    close(in_fasta)
    dest_dir
}

#' Prepare a seed file for building a BSgenome package
#'
#'
#' @param multifasta_path A character(1), specifying the path to a directory
#'   containing fasta.gz files for individual chromosome/scaffolds.
#' @param latin_name A character(1), specifying the Latin name of the species of
#'   the reference genome, such as "Homo sapiens" for the human.
#' @param common_name A character(1), specifying the common name of the species
#'   of the reference genome, such as "Human".
#' @param genome_version A character(1), specifying the genome build, such as
#'   "GRCh38".
#' @param seed_file_name A character(1), specifying the path to a seed file to
#'   be created for building a BSgenome package.
#' @param fasta_url A character(1), specifying the URL from where the reference
#'   genome fasta file is downloaded.
#' @param release_date  A character(1), specifying the genome assembly release
#'   date, such as "August 2020".
#' @param source A character(1), specifying the source of the genome fasta
#'   file, such as "Ensembl".
#' @param circ_seq A character(n), specifying identifiers of circular sequences
#'   in a genome, such as mitochondrial genome.
#' @param version A character(1), specifying the version of a BSgenome package
#'   to be build, such as "1.0.0".
#'
#' @return A character(2), the package name and the path to the seed file.
#' @noRd
#'
#' @examples
#' dest_dir <- tempdir()
#' genome_fasta <- file.path(dest_dir, "toy.example.fa")
#' in_fasta <- file(genome_fasta, open = "w")
#' writeLines(c(">chr1\nATCGCTGCGGATGCGG",
#'              ">chr2\nCCCCCCCCCCCGGAAA",
#'              ">chrM\nATACTACTGGA"), in_fasta)
#' close(in_fasta)
#'
#' fa_dir <- generate_multifasta(
#'     genome_fasta = genome_fasta,
#'     dest_dir = dest_dir,
#'     prefix = "none"
#' )
#'
#' seed_file <- generate_seed_file(
#'     multifasta_path = fa_dir,
#'     latin_name = "Homo sapiens",
#'     common_name = "Human",
#'     genome_version = "GRCh38",
#'     seed_file_name = file.path(
#'         dest_dir,
#'         "human.genome.seed.txt"
#'     ),
#'     fasta_url = "http://ftp.ensembl.org/xxx.fa",
#'     release_date = "August 2007",
#'     source = "Ensembl",
#'     version = "1.0.0"
#' )
generate_seed_file <-
    function(
        multifasta_path = NULL,
        latin_name = NULL,
        common_name = NULL,
        genome_version = NULL,
        seed_file_name = NULL,
        fasta_url = NULL,
        release_date = "August 2007",
        source = "Ensembl",
        circ_seqs = c("chrM", "MT", "Pltd", "chrPltd"),
        version = "1.0.0") {
        if (is.null(multifasta_path) || is.null(latin_name) ||
            is.null(common_name) || is.null(genome_version) ||
            is.null(seed_file_name) || is.null(fasta_url)) {
            stop("All arguments except source and version are required")
        }
        if (!dir.exists(multifasta_path)) {
            stop("Path to multifasta ", multifasta_path, " doesn't exist")
        }
        chr_fa_files <- dir(multifasta_path, ".fa.gz$")
        if (length(chr_fa_files) < 1) {
            stop(
                "There are not multiple fasta files in the directory ",
                multifasta_path
            )
        }
        seed_dir <- dirname(seed_file_name)
        if (!dir.exists(seed_dir)) {
            dir.create(seed_dir, recursive = TRUE)
        }
        seed_fh <- file(seed_file_name, open = "w")
        latin_name <- trimws(latin_name, which = "both")
        BSgenomeObjname <- gsub("^(.).*\\s+(.+)", "\\1\\2", latin_name)
        package_name <- paste("BSgenome", BSgenomeObjname, source,
                              genome_version,
                              sep = "."
        )
        writeLines(paste0("Package: ", package_name), con = seed_fh)
        writeLines(paste(
            "Title: Full genome sequences for",
            latin_name,
            paste0("(", source),
            " version ",
            paste0(genome_version, ")")
        ), con = seed_fh)
        writeLines(
            paste(
                "Description: Full genome sequences for",
                latin_name,
                paste0("(", common_name, ")"),
                "as provided by",
                source,
                paste0("(", genome_version, ")"),
                "and stored in Biostrings objects."
            ), con = seed_fh
        )
        writeLines(paste0("Version: ", version), con = seed_fh)
        writeLines(paste0("organism: ", latin_name), con = seed_fh)
        writeLines(paste0("common_name: ", common_name), con = seed_fh)
        writeLines(paste0("provider: ", source), con = seed_fh)
        writeLines(paste0("release_date: ", release_date), con = seed_fh)
        writeLines(paste0("genome: ", genome_version), con = seed_fh)
        writeLines(paste0("source_url: ", fasta_url), con = seed_fh)
        writeLines(paste0("BSgenomeObjname: ", BSgenomeObjname), con = seed_fh)
        writeLines(paste0(
            "organism_biocview: ",
            gsub("\\s+", "_", latin_name, perl = TRUE)), con = seed_fh)
        chromosome_names <-
            gsub(".fa.gz$", "", dir(multifasta_path, "fa.gz$"))
        writeLines(paste0(
            'seqnames: c("',
            paste(chromosome_names, collapse = '","'),
            '")'
        ), con = seed_fh)

        circ_seqs <- circ_seqs[circ_seqs %in% chromosome_names]
        if (length(circ_seqs) >=1)
        {
            writeLines(paste0('circ_seqs: c("',
                              paste(circ_seqs, collapse = '","'), '")'),
                       con = seed_fh)
        }

        writeLines(paste0("seqs_srcdir: ", multifasta_path),
                   con = seed_fh)
        writeLines(paste0("seqfiles_suffix: .fa.gz"), con = seed_fh)
        close(seed_fh)
        c(package_name, seed_file_name)
    }

#' Create and install a BSgenome package
#'
#' Based on a seed file, a BSgenome package is created and installed
#'
#' @param seed_file A character(1), specifying a path to a seed file
#' @param dest_dir A character(1), specifying a directory where a tar-ball
#'   of a BSgenome package is created.
#'
#' @return A BSgenome package name
#' @noRd
#'
#' @examples
#' if (TRUE) {
#'     dest_dir <- tempdir()
#'     genome_fasta <- file.path(dest_dir, "toy.example.fa")
#'     in_fasta <- file(genome_fasta, open = "w")
#' writeLines(c(">chr1\nATCGCTGCGGATGCGG",
#'               ">chr2\nCCCCCCCCCCCGGAAA",
#'               ">chrM\nATACTACTGGA"), in_fasta)
#'     close(in_fasta)
#'
#'     fa_dir <- generate_multifasta(
#'         genome_fasta = genome_fasta,
#'         dest_dir = dest_dir,
#'         prefix = "none"
#'     )
#'
#'     seed_file <- generate_seed_file(
#'         multifasta_path = fa_dir,
#'         latin_name = "Homo sapiens",
#'         common_name = "Human",
#'         genome_version = "GRCh38",
#'         seed_file_name = file.path(
#'             dest_dir,
#'             "human.genome.seed.txt"
#'         ),
#'         fasta_url = "http://ftp.ensembl.org/xxx.fa",
#'         release_date = "August 2007",
#'         source = "Ensembl",
#'         version = "1.0.0"
#'     )
#'
#'     forge_BSgenome(seed_file[2], dest_dir = dest_dir)
#' }
forge_BSgenome <-
    function(
        seed_file = NULL,
        dest_dir = tempdir()) {
        if (!file.exists(seed_file)) {
            stop("seed_file are required")
        }
        if (!dir.exists(dest_dir)) {
            dir.create(dest_dir, recursive = TRUE)
        }

        pkgname <- gsub(
            "Package:\\s*([^\\s]+)",
            "\\1",
            trimws(readLines(seed_file, n = 1))
        )
        if (dir.exists(file.path(dest_dir, pkgname))) {
            unlink(file.path(dest_dir, pkgname), recursive = TRUE)
        }

        ## crate a BSgenome package from a seed file
        BSgenomeForge::forgeBSgenomeDataPkg(seed_file, destdir = dest_dir)

        ## check, build and install package
        ## OR install using command line
        # R CMD build BSgenome.Hsapiens.Ensembl.GRCh38
        # R CMD check BSgenome.Hsapiens.Ensembl.GRCh38_1.0.0.tar.gz
        # devtools::check(file.path(dest_dir, pkgname))
        # devtools::build(file.path(dest_dir, pkgname), path = dest_dir)
        # devtools::install(file.path(dest_dir, pkgname))
        pkgname
    }
