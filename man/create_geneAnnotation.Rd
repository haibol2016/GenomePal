% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_geneAnnotation.R
\name{create_geneAnnotation}
\alias{create_geneAnnotation}
\title{Title}
\usage{
create_geneAnnotation(
  TxDb = NULL,
  geneID_symbol,
  filterChr = c("chrM", "chrY", "MT", "Y", "Pltd", "chrMT", "chrPltd"),
  flank = 2000,
  promoterRange = c(upstream = 2000, downstream = 2000)
)
}
\arguments{
\item{TxDb}{An object of \link[GenomicFeatures:TxDb-class]{GenomicFeatures::TxDb}. Chromosome lengths
must be available, which is used to get rid of out-of-bound GRanges
expanding/slopping.}

\item{geneID_symbol}{A data frame with two columns named \code{gene_id}
and \code{symbol}.}

\item{filterChr}{A character(n), specifying the chromosomes/scaffolds to be
excluded. \code{TSSs} in the \code{geneAnnotation} will be excluded if on the
\code{filterChr}.}

\item{flank}{An integer(1), specifying how many bases to expanded in both
directions on both sides of TSSs.}

\item{promoterRange}{A integer(2), specifying how many bases to expanded
upstream and downstream of TSSs to define the promoter regions of genes.}
}
\value{
An object of SimpleList with three elements as described below.
\describe{
\item{genes}{A GRanges defineing gene regions}
\item{exons}{A GRanges defineing exonic regions}
\item{TSS}{A GRanges defining TSSs}
}
}
\description{
Title
}
\examples{
library("txdbmaker")
gtf <- system.file("extdata", "example.gtf.gz",
                   package = "GenomePal",
                   mustWork = TRUE)
geneID_symbol <- get_geneID_symbol(gtf = gtf)

chrominfo <- data.frame(chrom = c('21','MT'),
                        length=c(46709983, 16569))
TxDb <- makeTxDbFromGFF(file = gtf,
                        format = "gtf",
                        chrominfo = chrominfo,
                        organism = "Homo sapiens",
                        taxonomyId = 9606,
                        circ_seqs = "MT")

## Or using seqinfo of aBSgenome
# library("BSgenome.Hsapiens.UCSC.hg38")
# seqlevelsStyle(BSgenome.Hsapiens.UCSC.hg38) <- "Ensembl"
# seqinfo_hg38 <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
# seqinfo_hg38 <- seqinfo_hg38[c("21", "MT")]
# TxDb <- makeTxDbFromGFF(file = gtf,
#                        format = "gtf",
#                        chrominfo = seqinfo_hg38,
#                        organism = "Homo sapiens",
#                        taxonomyId = 9606)
geneAnno <- create_geneAnnotation(TxDb = TxDb,
                                  geneID_symbol = geneID_symbol,
                                  filterChr = "MT")
}
