## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----filter, message=FALSE----------------------------------------------------
library(MultiAmplicon)
## And we'll also use some dada2 functions directly
library(dada2)
path <- "download_sra" ## change according to where you downloaded

fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE)

fastqF <- grep("_1.fastq.gz", fastqFiles, value = TRUE)
fastqR <- grep("_2.fastq.gz", fastqFiles, value = TRUE)

samples <- gsub("_1.fastq\\.gz", "\\1", basename(fastqF))

filt_path <- "filtered_sra"
if(!file_test("-d", filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples

## some files will be filtered out completely, therefore allowing 50
## files less present and still don't redo filtering
redo <- sum(file.exists(fastqF)) -
    sum(file.exists(filtFs)) > 50

if(redo){
    lapply(seq_along(fastqF),  function (i) {
        filterAndTrim(fastqF[i], filtFs[i], fastqR[i], filtRs[i],
                      truncLen=c(170,170), minLen=c(170,170), 
                      maxN=0, maxEE=2, truncQ=2, 
                      compress=TRUE, verbose=TRUE)
    })
}

names(filtFs) <- names(filtRs) <- samples

files <- PairedReadFileSet(filtFs, filtRs)

## ----prepPrimers, message=FALSE-----------------------------------------------
primer.file <- system.file("extdata", "real_world_primers.csv",
                           package = "MultiAmplicon")

ptable <- read.csv(primer.file, sep=",", header=TRUE, stringsAsFactors=FALSE)

primerF <- ptable[, "TS.SequenceF"]
primerR <- ptable[, "TS.SequenceR"]

names(primerF) <- as.character(ptable[, "corrected.NameF"])
names(primerR) <- as.character(ptable[, "corrected.NameR"])

primers <- PrimerPairsSet(primerF, primerR)

MA <- MultiAmplicon(primers, files)

## ----sortAmps, message=FALSE--------------------------------------------------
filedir <- "stratified_files"
if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
MA <- sortAmplicons(MA, n=1e+05, filedir=filedir)

## ----pipeline, message=FALSE--------------------------------------------------
errF <- learnErrors(unlist(getStratifiedFilesF(MA)), nbase=1e8,
                    verbose=0)
errR <- learnErrors(unlist(getStratifiedFilesR(MA)), nbase=1e8,
                    verbose=0)

MA <- dadaMulti(MA, Ferr=errF, Rerr=errR,  pool=FALSE,
                verbose=0)

## ----merger, message=FALSE----------------------------------------------------
MA <- mergeMulti(MA)

propMerged <- MultiAmplicon::calcPropMerged(MA)

summary(propMerged)
table(propMerged<0.8)

## ----tabulator, mesage=FALSE--------------------------------------------------
MA <- mergeMulti(MA, justConcatenate=propMerged<0.8)
MA <- makeSequenceTableMulti(MA)
MA <- removeChimeraMulti(MA)

## ----tracking, mesage=FALSE, fig.show='hold'----------------------------------
tracking <- getPipelineSummary(MA)
plotPipelineSummary(tracking)

## ----plottrack, mesage=FALSE, fig.show='hold'---------------------------------
library(ggplot2)
plotPipelineSummary(tracking) + scale_y_log10()

## ----handover, mesage=FALSE---------------------------------------------------
library(phyloseq)
PH <- toPhyloseq(MA, samples=colnames(MA))
PHlist <- toPhyloseq(MA, samples=colnames(MA), multi2Single=FALSE)

