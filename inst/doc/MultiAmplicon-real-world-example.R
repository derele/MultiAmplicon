## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----filter, message=FALSE-----------------------------------------------
library(MultiAmplicon)

path <- "~/download_sra" ## change according to where you downloaded

fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE)

fastqF <- grep("_1.fastq.gz", fastqFiles, value = TRUE)
fastqR <- grep("_2.fastq.gz", fastqFiles, value = TRUE)

samples <- gsub("_1.fastq\\.gz", "\\1", basename(fastqF))

filt_path <- "~/filtered_sra/"
if(!file_test("-d", filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples

## some files will be filtered out completely, therefore allowing 50
## files less present and still don't redo filtering
if(sum(file.exists(fastqF)) -
   sum(file.exists(filtFs)) > 50){
    lapply(seq_along(fastqF),  function (i) {
        filterAndTrim(fastqF[i], filtFs[i], fastqR[i], filtRs[i],
                      truncLen=c(170,170), minLen=c(170,170), 
                      maxN=0, maxEE=2, truncQ=2, 
                      compress=TRUE, verbose=TRUE)
    })
}

names(filtFs) <- names(filtRs) <- samples

files <- PairedReadFileSet(filtFs, filtRs)

## ----prepPrimers---------------------------------------------------------
primer.file <- system.file("extdata", "real_world_primers.csv",
                           package = "MultiAmplicon")

ptable <- read.csv(primer.file, sep=",", header=TRUE, stringsAsFactors=FALSE)

primerF <- ptable[, "TS.SequenceF"]
primerR <- ptable[, "TS.SequenceR"]

names(primerF) <- as.character(ptable[, "corrected.NameF"])
names(primerR) <- as.character(ptable[, "corrected.NameR"])

primers <- PrimerPairsSet(primerF, primerR)

MA <- MultiAmplicon(primers, files)

## ----sortAmps------------------------------------------------------------
MA <- sortAmplicons(MA)

## ----pipeline------------------------------------------------------------
errF <- learnErrors(unlist(stratifiedFilesF(MA)), nread=1e6,
                    verbose=0)
errR <- learnErrors(unlist(stratifiedFilesR(MA)), nread=1e6,
                    verbose=0)

MA <- derepMulti(MA, mc.cores=1) 
MA <- dadaMulti(MA, Ferr=errF, Rerr=errR,  pool=FALSE)

## ----merger--------------------------------------------------------------
MA <- mergeMulti(MA)

propMerged <- MultiAmplicon::calcPropMerged(MA)

summary(propMerged)
table(propMerged<0.8)

## ----tabulator-----------------------------------------------------------
MA <- mergeMulti(MA, justConcatenate=propMerged<0.8)
MA <- sequenceTableMulti(MA)
MA <- noChimeMulti(MA, mc.cores=1)

