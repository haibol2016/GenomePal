#' Extract gene identifiers and symbols from GTF
#'
#' @param gtf A character(1), specifying a path to an Ensembl GTF file. Many
#'   UCSC GTF files are malformatted, thus should not be used.
#' @param unique_symbol A logical(1), whether only return genes with unique
#'   symbols. Default is FALSE.
#' @return A data frame with two columns: `gene_id` and `symbol`.
#'
#' @importFrom collections ordered_dict
#' @importFrom utils read.delim globalVariables
#' @export
#'
#' @examples
#' gtf <- system.file("extdata", "example.gtf.gz",
#'                    package = "GenomePal",
#'                    mustWork = TRUE)
#' geneID_symbol <- get_geneID_symbol(gtf = gtf)
#'
get_geneID_symbol <-
    function(gtf = NULL,
             unique_symbol = FALSE) {
        if (is.null(gtf) || !file.exists(gtf))
        {
            stop("An existing GTF file is required")
        }
        if (grepl(".gtf.gz$", gtf))
        {
            in_gtf <- gzfile(gtf, open = "rt")
        } else if (grepl(".gtf$", gtf)) {
            in_gtf <- file(gtf, open = "r")
        } else {
            stop(
                "It seems the GTF file is not a GTF file ",
                "which should with an extension '.gtf', or '.gtf.gz'"
            )
        }

        ## this can also be done using rtracklayer::import
        id2symbol_dict <- ordered_dict()
        gtf_full <- read.delim(
            in_gtf,
            header = FALSE,
            as.is = TRUE,
            comment.char = "#",
            quote = ""
        )
        close(in_gtf)

        gtf_attr <- gtf_full[gtf_full[, 3] == "gene", 9]
        if (length(gtf_attr) < 1)
        {
            message("There is no entries for genes in the GTF\n",
                    "Use transcript entries")
            gtf <- gtf_full[gtf_full[, 3] == "transcript", 9]
            if (length(gtf) < 1) {
                stop("gtf is malformed, please doulbe check it")
            }
        }
        null <- lapply(gtf_attr, function(.x) {
            if (grepl("gene_id", .x)) {
                gene_id <- gsub('.*gene_id\\s+"([^".]+).+',
                                "\\1",
                                .x,
                                perl = TRUE)
                if (grepl("gene_name", .x))
                {
                    gene_symbol <- gsub('.*gene_name\\s+"([^"]+).*',
                                        "\\1",
                                        .x,
                                        perl = TRUE)
                    if (grepl("^ENS.*?[GT]\\d+", gene_symbol, perl = TRUE) ||
                        grepl("^(WBGene|FBgn|ENSDARG)\\d+",
                              gene_symbol, perl = TRUE) ||
                        grepl("^\\d+$", gene_symbol, perl = TRUE))
                    {
                        gene_symbol <- NA_character_
                    }
                } else {
                    gene_symbol <- NA_character_
                }
                id2symbol_dict$set(gene_id, gene_symbol)
            }
        })

        ## check gene_id type
        gene_ids <- unlist(id2symbol_dict$keys())
        id_type <- {
            if (grepl("^(ENS.+?|WBGene|FBgn|ENSDARG)\\d+",
                      gene_ids[1], perl = TRUE)) {
                "ensembl_gene_id"
            }
            else if (grepl("^\\d+$", gene_ids[1], perl = TRUE)) {
                "entrezgene_id"
            }
            else if (!grepl("^ENS.*?T\\d+", gene_ids[1], perl = TRUE) &&
                     grepl("[a-zA-Z0-9]+", gene_ids[1], perl = TRUE)) {
                "gene_name"
            }
            else {
                "unknown"
            }
        }

        if (id_type == "unknown")
        {
            stop("Unknown gene ID type!")
        } else if (id_type == "gene_name") {
            id2symbol <- data.frame(gene_id = gene_ids,
                                    symbol = gene_ids)
        } else {
            id2symbol <- data.frame(
                gene_id = unlist(id2symbol_dict$keys()),
                symbol = unlist(id2symbol_dict$values())
            )
            id2symbol$symbol <- ifelse(is.na(id2symbol$symbol),
                                       id2symbol$gene_id,
                                       id2symbol$symbol)
        }
        if (unique_symbol) {
            id2symbol <- id2symbol[!duplicated(id2symbol$symbol), ]
        }
        invisible(id2symbol)
    }
