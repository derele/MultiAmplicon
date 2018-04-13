## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----buildup, message=FALSE----------------------------------------------
library(MultiAmplicon)

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
errF <- learnErrors(unlist(stratifiedFilesF(MA1)), verbose=0)
errR <- learnErrors(unlist(stratifiedFilesR(MA1)), verbose=0)

## ----pipeline, cache=TRUE, message=FALSE---------------------------------
MA2 <- derepMulti(MA.sub, mc.cores=1) 

MA3 <- dadaMulti(MA2, Ferr=errF, Rerr=errR,  pool=FALSE)

MA4 <- mergeMulti(MA3, justConcatenate=TRUE)

MA5 <- sequenceTableMulti(MA4)

MA6 <- noChimeMulti(MA5, mc.cores=1)

## ----altdada, cache=TRUE, message=FALSE----------------------------------
MA.alt <- dadaMulti(MA2, selfConsist=TRUE, pool=FALSE)

MA.alt <- mergeMulti(MA.alt, justConcatenate=TRUE)

MA.alt <- sequenceTableMulti(MA.alt)

MA.alt <- noChimeMulti(MA.alt, mc.cores=1)

## ----compdata, cache=TRUE, message=FALSE---------------------------------
together <- getSequencesFromTable(MA6)
alt <- getSequencesFromTable(MA.alt)

lapply(seq_along(together), function (i){
    setequal(alt[[i]], together[[i]])
    table(alt[[i]]%in%together[[i]])
})

## ----mixmerge, cache=TRUE------------------------------------------------
MA.mixed <- mergeMulti(MA3, justConcatenate=c(TRUE, FALSE, FALSE, TRUE))

## ----mixdada, cache=TRUE-------------------------------------------------
MA.mixed <- dadaMulti(MA2, Ferr=c(NULL, errF, NULL, errF),
                      Rerr=c(NULL, errF, NULL, errF),
                      selfConsist=c(TRUE, FALSE, TRUE, FALSE))

