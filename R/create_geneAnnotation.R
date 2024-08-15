#' Title
#'
#' @param TxDb An object of [GenomicFeatures::TxDb-class]. Chromosome lengths
#'   must be available, which is used to get rid of out-of-bound GRanges
#'   expanding/slopping.
#' @param geneID_symbol A data frame with two columns named `gene_id`
#'   and `symbol`.
#' @param flank  An integer(1), specifying how many bases to expanded in both
#'   directions on both sides of TSSs.
#' @param promoterRange A integer(2), specifying how many bases to expanded
#'   upstream and downstream of TSSs to define the promoter regions of genes.
#' @param filterChr A character(n), specifying the chromosomes/scaffolds to be
#'   excluded. `TSSs` in the `geneAnnotation` will be excluded if on the
#'   `filterChr`.
#'
#' @return An object of SimpleList with three elements as described below.
#'   \describe{
#'       \item{genes}{A GRanges defineing gene regions}
#'       \item{exons}{A GRanges defineing exonic regions}
#'       \item{TSS}{A GRanges defining TSSs}
#'   }
#' @export
#' @importFrom GenomicFeatures genes exonsBy transcripts transcriptsBy
#'
#' @importFrom plyranges remove_names mutate select filter stretch anchor_5p
#'   anchor_3p
#' @importFrom GenomicRanges resize mcols
#' @importFrom GenomeInfoDb seqlevels seqlevels<- sortSeqlevels
#' @importFrom S4Vectors SimpleList
#' @importFrom magrittr  %>%
#'
#' @examples
#' library("txdbmaker")
#' gtf <- system.file("extdata", "example.gtf.gz",
#'                    package = "GenomePal",
#'                    mustWork = TRUE)
#' geneID_symbol <- get_geneID_symbol(gtf = gtf)
#'
#' chrominfo <- data.frame(chrom = c('21','MT'),
#'                         length=c(46709983, 16569))
#' TxDb <- makeTxDbFromGFF(file = gtf,
#'                         format = "gtf",
#'                         chrominfo = chrominfo,
#'                         organism = "Homo sapiens",
#'                         taxonomyId = 9606,
#'                         circ_seqs = "MT")
#'
#' ## Or using seqinfo of aBSgenome
#' # library("BSgenome.Hsapiens.UCSC.hg38")
#' # seqlevelsStyle(BSgenome.Hsapiens.UCSC.hg38) <- "Ensembl"
#' # seqinfo_hg38 <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
#' # seqinfo_hg38 <- seqinfo_hg38[c("21", "MT")]
#' # TxDb <- makeTxDbFromGFF(file = gtf,
#' #                        format = "gtf",
#' #                        chrominfo = seqinfo_hg38,
#' #                        organism = "Homo sapiens",
#' #                        taxonomyId = 9606)
#' geneAnno <- create_geneAnnotation(TxDb = TxDb,
#'                                   geneID_symbol = geneID_symbol,
#'                                   filterChr = "MT")
create_geneAnnotation <-
    function(TxDb = NULL,
             geneID_symbol,
             filterChr = c("chrM", "chrY", "MT", "Y", "Pltd",
                           "chrMT", "chrPltd"),
             flank = 2000,
             promoterRange =
                 c(upstream = 2000, downstream = 2000))
{
    if (is.null(TxDb) || !is(TxDb, "TxDb")) {
        stop("'TxDb' must ba a TxDb object")
    }
    if (any(is.na(seqlengths(TxDb)))) {
        stop("Chromosome lengths are not available in the TxDb")
    }

    if (missing(geneID_symbol)) {
        stop("'geneID_symbol' is required")
    } else if (!is.data.frame(geneID_symbol) ||
        !setequal(colnames(geneID_symbol), c("gene_id", "symbol"))) {
        stop("geneID_symbol must be a dataframe with colnames: ",
        "'gene_id', 'symbol'")
    }

    ## geneID_symbol: col1: gene_id; col2: gene_symbol
    symbol <- geneID_symbol$symbol
    names(symbol) <- geneID_symbol$gene_id

    ## filter seqnames of no interest, such as mitochondrial genome
    seqlevels_all <- seqlevels(TxDb)
    seqlevels(TxDb) <- seqlevels_all[!seqlevels_all %in% filterChr]

    ## GRanges for genes
    genes <- GenomicFeatures::genes(TxDb) %>%
        plyranges::remove_names()
    genes$symbol <- symbol[genes$gene_id]  # more stable
    genes <- sort(sortSeqlevels(genes), ignore.strand = TRUE)

    ## get all transcripts
    tx <- unlist(transcriptsBy(TxDb, by = "gene")) %>%
        plyranges::mutate(gene = names(.)) %>%
        plyranges::remove_names() %>%
        plyranges::select(-c("tx_id")) %>%
        data.frame()
    tx_gene <- tx$gene
    names(tx_gene) <-tx$tx_name
    rm("tx")

    ## Create GRanges for exons
    exons <- unlist(exonsBy(TxDb,
                            by = "tx",
                            use.names = TRUE)) %>%
        plyranges::mutate(tx_name = names(.)) %>%
        plyranges::remove_names()
    exons$gene_id <- tx_gene[exons$tx_name]
    exons <- exons %>%
        plyranges::filter(!is.na(gene_id))
    exons$symbol <-  symbol[exons$gene_id]
    exons <- exons %>%
        plyranges::select(-c("exon_id",
                             "exon_name",
                             "exon_rank",
                             "gene_id",
                             "tx_name"))
    exons <- sort(sortSeqlevels(exons), ignore.strand = TRUE)

    ## Create GRanges for TSS
    TSS <- unique(resize(GenomicFeatures::transcripts(TxDb),
                         width = 1,
                         fix = "start")) %>%
        plyranges::select(-c("tx_id"))

    ## remove genes whose promoters are close to the chromosome end
    ## (promoter regions in upstream and downstream [2000, 2000])
    gene_start <-  resize(genes, width = 1, fix ="start")
    gene_start_downstream <- plyranges::stretch(anchor_5p(gene_start),
                                     extend = promoterRange[2])
    gene_start_upstream <- plyranges::stretch(anchor_3p(gene_start),
                                   extend = promoterRange[1])
    downstream_out_of_bound_index <-
        .get_out_of_bound_index(gene_start_downstream)
    upstream_out_of_bound_index <-
        .get_out_of_bound_index(gene_start_upstream)
    gene_out_of_bound_index <- c(downstream_out_of_bound_index,
                                 upstream_out_of_bound_index)

    if (length(gene_out_of_bound_index) > 0)
    {
        genes <- genes[-c(gene_out_of_bound_index)]
    }

    ## remove exons with genes removed due to out of bound
    exons <- exons[mcols(exons)$symbol %in% mcols(genes)$symbol]

    ## remove TSSs which are close to the chromosome end (<=2000 bp)
    TSS_2kb_flank <- resize(TSS,
                            width = 2 * flank + 1,
                            fix = "center")
    TSS_out_of_bound_index <-
        .get_out_of_bound_index(TSS_2kb_flank)
    if (length(TSS_out_of_bound_index) > 0)  # otherwise get empty TSS
    {
        TSS <- TSS[-c(TSS_out_of_bound_index)]
    }

    ## drop unused seqlevels
    seqlevels(genes) <- seqlevelsInUse(genes)
    seqlevels(exons) <- seqlevelsInUse(exons)
    seqlevels(TSS)   <- seqlevelsInUse(TSS)
    geneAnnotation <- SimpleList(genes = genes,
                                 exons = exons,
                                 TSS = TSS)
    geneAnnotation
}

### copied from GenomicRanges
### Returns index of out-of-bound ranges located on non-circular sequences
### whose length is not NA. Works on a GenomicRanges or GAlignments object.

#' Get indices of out of bound GRanges
#'
#' @param x An object of [GenomicRanges::GRanges-class] with seqinfo.
#'
#' @return A vector of integer indices of out-of-bound GRanges
#' @noRd
#' @importFrom GenomeInfoDb isCircular
#' @importFrom GenomicRanges start end
#'
.get_out_of_bound_index <- function(x)
{
    if (length(x) == 0L)
        return(integer(0))
    x_seqnames_id <- as.integer(seqnames(x))
    x_seqlengths <- unname(seqlengths(x))
    seqlevel_is_circ <- unname(isCircular(x)) %in% TRUE
    seqlength_is_na <- is.na(x_seqlengths)
    seqlevel_has_bounds <- !(seqlevel_is_circ | seqlength_is_na)
    which(seqlevel_has_bounds[x_seqnames_id] &
              (start(x) < 1L | end(x) > x_seqlengths[x_seqnames_id]))
}

