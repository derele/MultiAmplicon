## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----buildup, message=FALSE---------------------------------------------------
library(MultiAmplicon)
## we'll also use some dada2 functions directly
library(dada2)

fastq.dir <- system.file("extdata", "fastq", package = "MultiAmplicon")
fastq.files <- list.files(fastq.dir, full.names=TRUE)

Ffastq.file <- fastq.files[grepl("F_filt", fastq.files)]
Rfastq.file <- fastq.files[grepl("R_filt", fastq.files)]

names(Rfastq.file) <- names(Ffastq.file) <-
    gsub("_F_filt\\.fastq\\.gz", "", basename(Ffastq.file))

PRF <- PairedReadFileSet(Ffastq.file, Rfastq.file)

primerF <- c(Amp1F = "AGAGTTTGATCCTGGCTCAG", Amp2F = "ACTCCTACGGGAGGCAGC",
             Amp3F = "GAATTGACGGAAGGGCACC", Amp4F = "YGGTGRTGCATGGCCGYT",
             Amp5F = "AAAAACCCCGGGGGGTTTTT", Amp6F = "AGAGTTTGATCCTGCCTCAG")

primerR <- c(Amp1R = "CTGCWGCCNCCCGTAGG", Amp2R = "GACTACHVGGGTATCTAATCC",
             Amp3R = "AAGGGCATCACAGACCTGTTAT", Amp4R = "TCCTTCTGCAGGTTCACCTAC",
             Amp5R = "AAAAACCCCGGGGGGTTTTT", Amp6R = "CCTACGGGTGGCAGATGCAG")

PPS <- PrimerPairsSet(primerF, primerR)

MA <- MultiAmplicon(PPS, PRF)

## ----sortAmplicon, message=FALSE----------------------------------------------
MA1 <- sortAmplicons(MA, filedir=tempfile("dir"))

## ----showRC, results='asis'---------------------------------------------------
knitr::kable(getRawCounts(MA1))

## ----plotAmpliconNumbers, fig.show='hold'-------------------------------------
clusters <- plotAmpliconNumbers(MA1)

## ----subsetCluster, result='asis'---------------------------------------------
two.clusters.row <- cutree(clusters$tree_row, k=2)
two.clusters.col <- cutree(clusters$tree_col, k=2)

knitr::kable(two.clusters.row)
knitr::kable(two.clusters.col)

MA.sub <- MA1[which(two.clusters.row==1), which(two.clusters.col==2)]

knitr::kable(getRawCounts(MA.sub))

## ----errEst, cache=TRUE, message=FALSE----------------------------------------
errF <- learnErrors(unlist(getStratifiedFilesF(MA1)), verbose=0)
errR <- learnErrors(unlist(getStratifiedFilesR(MA1)), verbose=0)

## ----pipeline, cache=TRUE, message=FALSE--------------------------------------
MA3 <- dadaMulti(MA1, Ferr=errF, Rerr=errR,  pool=FALSE)

MA4 <- mergeMulti(MA3, justConcatenate=TRUE)

MA5 <- makeSequenceTableMulti(MA4)

MA6 <- removeChimeraMulti(MA5, mc.cores=1)


## ----altdada, cache=TRUE, message=FALSE---------------------------------------
MA.alt <- dadaMulti(MA3, selfConsist=TRUE, pool=FALSE)

MA.alt <- mergeMulti(MA.alt, justConcatenate=TRUE)

MA.alt <- makeSequenceTableMulti(MA.alt)

MA.alt <- removeChimeraMulti(MA.alt, mc.cores=1)

## ----compdata, cache=TRUE, message=FALSE--------------------------------------
together <- getSequencesFromTable(MA6)
alt <- getSequencesFromTable(MA.alt)

lapply(seq_along(together), function (i){
  union <- union(alt[[i]], together[[i]])
  table("combined Inference"=union%in%together[[i]],
        "seperate amplicon"=union%in%alt[[i]])
})


## ----mixmerge, cache=TRUE-----------------------------------------------------
MA.mixed <- mergeMulti(MA3, justConcatenate=c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE))

## ----mixdada, cache=TRUE------------------------------------------------------
MA.mixed <- dadaMulti(MA3,
                      Ferr=c(NULL, errF,  NULL, errF, NULL, NULL),
                      Rerr=c(NULL, errR,  NULL, errR, NULL, NULL),
               selfConsist=c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE))

## ----taxonAnnot, cache=TRUE---------------------------------------------------
## echo=TRUE for this?
inputFasta <- system.file("extdata", "in.fasta", package = "MultiAmplicon")
outputBlast <- system.file("extdata", "out.blt", package = "MultiAmplicon")

MA7 <- blastTaxAnnot(MA6, db="/SAN/db/blastdb/nt/nt",
                     negative_gilist = "/SAN/db/blastdb/uncultured_gilist.txt",
                     taxonSQL="/SAN/db/taxonomy/taxonomizr.sql", num_threads=20,
                     infasta = inputFasta,
                     outblast = outputBlast)

## ----toPhloseq, cache=TRUE----------------------------------------------------
## echo=TRUE for this?
library(phyloseq)

PHY <- toPhyloseq(MA7, samples=colnames(MA7), multi2Single=TRUE)
PHY.list <- toPhyloseq(MA7, samples=colnames(MA7), multi2Single=FALSE)
PHYgenus <- tax_glom(PHY, "genus")

knitr::kable(table(tax_table(PHYgenus)[, "superkingdom"]))
knitr::kable(table(tax_table(PHYgenus)[, "phylum"]))

