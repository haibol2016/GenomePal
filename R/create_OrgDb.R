#' Create an OrgDb package
#'
#' Create an OrgDb package based on a data frame containing `gene ID` and
#' `gene symbol`.
#'
#' @inheritParams AnnotationForge::makeOrgPackage
#' @param geneID_symbol A data frame with two columns named "gene_id" and
#'   "symbol"
#' @param id_type A character(1), specifying the type of gene identifiers in
#'   the `geneID_symbol` data frame. It should be "GENEID" for NCBI Entrz gene
#'   identifiers, such as '8359'  or "ENSEMBL" for Ensembl gene identifiers
#'   such as 'ENSG00000196176'.
#' @param outDir A character(1), specifying a path where the package source
#'   should be assembled.
#'
#' @return A character(1), the name of and OrgDb package.
#' @importFrom AnnotationForge makeOrgPackage
#' @export
#'
#' @examples
#' gtf <- system.file("extdata", "example.gtf.gz",
#'                    package = "GenomePal",
#'                    mustWork = TRUE)
#' geneID_symbol <- get_geneID_symbol(gtf = gtf, unique_symbol = TRUE)
#' outDir <- tempdir()
#' files <- dir(outDir, pattern = "org.*.eg.*",
#'              full.name = TRUE)
#' if (length(files) >= 1) {
#'     unlink(files, force = TRUE, recursive = TRUE)
#' }
#' try(OrgDb <- create_OrgDb(geneID_symbol,
#'                       id_type = "ENSEMBL",
#'                       outDir = outDir))
create_OrgDb <-
    function(geneID_symbol,
             id_type = c("GENEID", "ENSEMBL"),
             goTable = NULL,
             outDir = ".",
             tax_id = "9606",
             genus = "Homo",
             species = "sapeins",
             version = "0.1",
             maintainer = "Some One <so@someplace.org>",
             author = "Some One <so@someplace.org>")
    {
        if (!is.data.frame(geneID_symbol) ||
            !setequal(colnames(geneID_symbol), c("gene_id", "symbol"))) {
            stop(
                "'geneID_symbol' is not a dataframe or the colnames are not ",
                "'gene_id' and 'symbol'"
            )
        }
        if (any(!is.character(c(
            tax_id, genus, species, version, maintainer, author
        )))) {
            stop(
                "All 'tax_id', 'genus', 'species', 'version',
             'maintainer', and 'author' must be character"
            )
        }

        id_type <- match.arg(id_type, choices = c("GENEID", "ENSEMBL"))

        fSym <- unique(geneID_symbol[, c("gene_id", "symbol")])
        colnames(fSym) <- c("GID", "SYMBOL")
        if (all(grepl("_", fSym$SYMBOL, perl = TRUE)))
        {
            fSym$SYMBOL <- gsub("_.+$", "", fSym$SYMBOL, perl = TRUE)
        }

        ensembl <- unique(geneID_symbol[, c(1, 1)])
        colnames(ensembl) <- {
            if (id_type == "ENSEMBL") {
                c("GID", "ENSEMBL")
            } else {
                c("GID", "GENEID")
            }
        }
        if (!dir.exists(outDir)) {
            dir.create(outDir, recursive = TRUE)
        }

        AnnotationForge::makeOrgPackage(
            gene_info = fSym,
            ensembl = ensembl,
            version = version,
            maintainer = maintainer,
            author = author,
            outputDir = outDir,
            tax_id = tax_id,
            genus = genus,
            species = species,
            goTable = goTable
        )
        package_name <-
            paste("org",
                  paste0(gsub("^(.).*", "\\1", genus, perl = TRUE),
                         tolower(species)), "eg.db",
                  sep = ".")
        invisible(package_name)
    }
