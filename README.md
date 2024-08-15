
# GenomePal

<!-- badges: start -->
<!-- badges: end -->

The goal of GenomePal is to extract gene IDs and symbols from an Ensembl GTF
files, build BSgeome from a multifasta file, generate OrgDb, geneAnnotation,
and genomeAnnotation for ArchR-based scATAC-seq data analysis using a custom
genome.

## Installation

You can install the development version of GenomePal like so:

``` r
library(devtools)
devtools::install_github("haibol2016/GenomePal")
```

## Example 1: how to generate a BSgenome package
```r
dest_dir <- tempdir()

# For a real example, download C. elegans genome fasta sequences from Ensembl
# url <- paste0("https://ftp.ensembl.org/pub/release-112/fasta/",
                "caenorhabditis_elegans/dna/Caenorhabditis_elegans.",
                "WBcel235.dna.toplevel.fa.gz")
# options(timeout = max(3000, getOption("timeout")))                
# genome_fasta <- file.path(dest_dir, "WBcel235.dna.toplevel.fa.gz")
# download.file(url, destile = genome_fasta),
                mode = "wb")

# For a toy example                
genome_fasta <- file.path(dest_dir, "toy.example.fa")
in_fasta <- file(genome_fasta, open = "w")
writeLines(c(">chr1\nATCGCTGCGGATGCGG",
          ">chr2\nCCCCCCCCCCCGGAAA",
          ">chrM\nATACTACTGGA"), in_fasta)
close(in_fasta)

make_BSgenome_pkg(
    genome_fasta = genome_fasta,
    dest_dir = file.path(dest_dir, "multifasta"),
    prefix = "none",
    latin_name = "Homo sapiens",
    common_name = "Human",
    genome_version = "GRCh38",
    seed_file_name = file.path(
        dest_dir,
        "human.genome.seed.txt"
    ),
    fasta_url = "http://ftp.ensembl.org/xxx.fa",
    release_date = "August 2007",
    source = "Ensembl",
    version = "1.0.0"
)
```

## Example 2: extract gene IDs and symbol from an Ensembl GTF file and 
create a OrgDb

```r
gtf <- system.file("extdata", "example.gtf.gz",
                   package = "GenomePal",
                   mustWork = TRUE)
geneID_symbol <- get_geneID_symbol(gtf = gtf, unique_symbol = TRUE)
outDir <- tempdir()
files <- dir(outDir, pattern = "org.*.eg.*",
             full.name = TRUE)
if (length(files) >= 1) {
    unlink(files, force = TRUE, recursive = TRUE)
}
try(OrgDb <- create_OrgDb(geneID_symbol,
                      id_type = "ENSEMBL",
                      outDir = outDir))
```


## Example 3
This is a basic example which shows you how to create geneAnnotation and 
genomeAnnotation for ArchR-based scATAC-seq data analysis using a 
custom genome.

``` r
library("txdbmaker")
library("BSgenome.Hsapiens.UCSC.hg38")
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
#'
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
#'
seqlevelsStyle(BSgenome.Hsapiens.UCSC.hg38) <- "ENSEMBL"
genomeAnno <- create_genomeAnnotation(BSgenome = BSgenome.Hsapiens.UCSC.hg38,
                                      geneAnnotation = geneAnno,
                                      blacklistBed = NULL,
                                      header = FALSE)
```

