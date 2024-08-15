#' Create genomeAnnotation for ArchR analysis
#'
#' Create a genomeAnnotation object for customize ArchR analysis
#'
#' @param BSgenome An object of [BSgenome::BSgenome-class].
#' @param geneAnnotation A object of SimpleList with three sublist neamed as
#'   "genes", "exons", and "TSS", such as an object outputted by the
#'   [create_geneAnnotation()] function. TSS
#' @param blacklistBed A character(1), tab-delimited BED file without a header
#'   line.
#' @param header A logical (1), indicating whether the file contains the
#'   names of the variables as its first line.
#'
#' @return An object of SimpleList with three elements as described below.
#'   \describe{
#'       \item{genome}{the name of a BSgenome package}
#'       \item{chromSizes}{A GRanges storing the chromosomal sizes}
#'       \item{blacklist}{A GRanges storing a blacklist for genomic regions to
#'              be excluded}
#'   }
#' @export
#'
#' @importFrom methods is
#' @importFrom GenomicRanges makeGRangesFromDataFrame GRanges
#' @importFrom GenomeInfoDb seqnames seqlevels seqlevelsInUse seqlevels<-
#'   seqlengths
#' @importFrom S4Vectors SimpleList
#' @examples
#' library("txdbmaker")
#' library("BSgenome.Hsapiens.UCSC.hg38")
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
#'
#' seqlevelsStyle(BSgenome.Hsapiens.UCSC.hg38) <- "ENSEMBL"
#' genomeAnno <- create_genomeAnnotation(BSgenome = BSgenome.Hsapiens.UCSC.hg38,
#'                                       geneAnnotation = geneAnno,
#'                                       blacklistBed = NULL,
#'                                       header = FALSE)
#'
create_genomeAnnotation <-
    function(BSgenome = NULL,
             geneAnnotation = NULL,
             blacklistBed = NULL,
             header = FALSE)
{
    if (is.null(BSgenome) || !is(BSgenome, "BSgenome")) {
        stop("'BSgenome' must be a BSgenome object")
    }
    if (is.null(geneAnnotation) ||
        !is(geneAnnotation, "SimpleList") ||
        !setequal(names(geneAnnotation), c("genes", "exons", "TSS"))) {
        stop("'geneAnnotation' must be a SimpleList ",
             "containing GRanges for genes, exons and TSSs")
    }

    chrom_len <- seqlengths(BSgenome)
    chr_df <- data.frame(seqnames = names(chrom_len),
                         start = 1,
                         end = unname(chrom_len))
    chromSizes <- makeGRangesFromDataFrame(chr_df)

    ## filter extra chromosomes/scaffolds so no ArrowFile generated for them
    tss_chr <- unique(as.character(seqnames(geneAnnotation$TSS)))

    ## filter chromSizes
    seqlevels(chromSizes, pruning.mode="coarse") <- tss_chr
    seqlevels(chromSizes) <- seqlevelsInUse(chromSizes)

    if (!is.null(blacklistBed))
    {
        if (grepl(".bed.gz$", blacklistBed))
        {
            blacklist <- gzfile(blacklistBed, open ="rt")
        } else if (grepl(".bed$", blacklistBed)) {
            blacklist <- file(blacklistBed, open = "r")
        } else {
            stop("It seems 'blacklist' is not a bed file ",
                 "which should end with an extension '.bed', or '.bed.gz'")
        }
        blacklist_df <- read.delim(blacklist,
                                   header = header,
                                   comment.char = "#",
                                   as.is = TRUE)
        close(blacklist)

        if (ncol(blacklist_df) < 3 ||
            any(!blacklist_df[, 1] %in% names(chrom_len)) ||
            any(!is.numeric(blacklist_df[, 2])) ||
            any(!is.numeric(blacklist_df[, 3])) ||
            any(blacklist_df[, 2] >= blacklist_df[, 3]))
        {
            stop("blacklistBed is not a valid BED file\n",
                 "Please make sure the chromosome names of ",
                 "the blacklist are a subset of those of the BSgenome.")
        }
        colnames(blacklist_df)[1:3] <- c("seqnames", "start", "end")
        blacklist <- makeGRangesFromDataFrame(blacklist_df,
                                              starts.in.df.are.0based = TRUE)
        ## filter blacklist
        seqlevels(blacklist, pruning.mode="coarse") <- tss_chr
        seqlevels(blacklist) <- seqlevelsInUse(blacklist)
    } else {
        blacklist <- GRanges()
    }

    genomeAnnotation <- SimpleList(genome = BSgenome@pkgname,
                                   chromSizes = chromSizes,
                                   blacklist = blacklist)
    genomeAnnotation
}
