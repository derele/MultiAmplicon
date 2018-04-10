## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----buildup, message=FALSE----------------------------------------------
library(MultiAmplicon)

primerF <- c(Amp1F = "AGAGTTTGATCCTGGCTCAG", Amp2F = "ACTCCTACGGGAGGCAGC",
             Amp3F = "GAATTGACGGAAGGGCACC", Amp4F = "YGGTGRTGCATGGCCGYT",
             Amp5F = "AAAAACCCCGGGGGGTTTTT", Amp6F = "AGAGTTTGATCCTGCCTCAG")

primerR <- c(Amp1R = "CTGCWGCCNCCCGTAGG", Amp2R = "GACTACHVGGGTATCTAATCC",
             Amp3R = "AAGGGCATCACAGACCTGTTAT", Amp4R = "TCCTTCTGCAGGTTCACCTAC",
             Amp5R = "AAAAACCCCGGGGGGTTTTT", Amp6R = "CCTACGGGTGGCAGATGCAG")

PPS <- PrimerPairsSet(primerF, primerR)

fastq.dir <- system.file("extdata", "fastq", package = "MultiAmplicon")
fastq.files <- list.files(fastq.dir, full.names=TRUE)

Ffastq.file <- fastq.files[grepl("F_filt", fastq.files)]
Rfastq.file <- fastq.files[grepl("R_filt", fastq.files)]

PRF <- PairedReadFileSet(Ffastq.file, Rfastq.file)

MA <- MultiAmplicon(PPS, PRF)

## ----sortAmplicon, message=FALSE-----------------------------------------
MA1 <- sortAmplicons(MA)

## ----showRC, results='asis'----------------------------------------------
knitr::kable(rawCounts(MA1))

## ----plotAmpliconNumbers, fig.show='hold'--------------------------------
clusters <- plotAmpliconNumbers(MA1)

## ----subsetCluster, result='asis'----------------------------------------
two.clusters.row <- cutree(clusters$tree_row, k=2)
two.clusters.col <- cutree(clusters$tree_col, k=2)

knitr::kable(two.clusters.row)
knitr::kable(two.clusters.col)

MA.sub <- MA1[which(two.clusters.row==1), which(two.clusters.col==2)]

knitr::kable(rawCounts(MA.sub))

## ----errEst, cache=TRUE, message=FALSE-----------------------------------
errF <- learnErrors(unlist(stratifiedFilesF(MA1)), nread=1e6,
                    verbose=0)
errR <- learnErrors(unlist(stratifiedFilesR(MA1)), nread=1e6,
                    verbose=0)

## ----pipeline, cache=TRUE, message=FALSE---------------------------------
MA2 <- derepMulti(MA.sub, mc.cores=1) 

MA3 <- dadaMulti(MA2, Ferr=errF, Rerr=errR,  pool=FALSE)

MA4 <- mergeMulti(MA3)

MA5 <- sequenceTableMulti(MA4)

MA6 <- noChimeMulti(MA5, mc.cores=1)

## ----altdada, cache=TRUE, message=FALSE----------------------------------
MA.alt <- dadaMulti(MA2, selfConsist=TRUE, pool=FALSE)

MA.alt <- mergeMulti(MA.alt)

MA.alt <- sequenceTableMulti(MA.alt)

MA.alt <- noChimeMulti(MA.alt, mc.cores=1)

## ----mixdata, cache=TRUE-------------------------------------------------
MA.mixed <- mergeMulti(MA3, justConcatenate=c(TRUE, FALSE, FALSE, TRUE))

