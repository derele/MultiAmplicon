
## Primers used in the arrays, primer pairs in single processin are part of this
ptable <- read.csv(file = "/SAN/Victors_playground/Metabarcoding/AA_Hyena/primer_list.csv",
                   sep=",", header=TRUE, stringsAsFactors=FALSE)
primerF <- ptable[, "Seq_F"]
primerR <- ptable[, "Seq_R"]
names(primerF) <- as.character(ptable[, "Name_F"])
names(primerR) <- as.character(ptable[, "Name_R"])
primer <- PrimerPairsSet(primerF, primerR)


### working with a small subset of Hyena data
fastq.dir <- "/SAN/MA_tests/fastqHy" 
fastq.files <- list.files(fastq.dir, full.names=TRUE)

fastqF <- grep("_R1_001.fastq.gz", fastq.files, value = TRUE)
fastqR <- grep("_R2_001.fastq.gz", fastq.files, value = TRUE)

samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF))

filt_path <- "/SAN/MA_tests/fastqHy/filtered"

if(!file_test("-d", filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples

filter.track <- lapply(seq_along(fastqF),  function (i) {
    filterAndTrim(fastqF[i], filtFs[i], fastqR[i], filtRs[i],
                  truncLen=c(220,200), minLen=c(220,200), 
                  maxN=0, maxEE=2, truncQ=2, 
                  compress=TRUE, verbose=TRUE,
                  matchIDs=TRUE) ## forward and reverse not matching otherwise 
})

PRF <- PairedReadFileSet(filtFs, filtRs)

MA <- MultiAmplicon(primer, PRF)

MA1 <- sortAmplicons(MA, filedir=tempfile())

MA3 <- dadaMulti(MA1, selfConsist=TRUE, pool=FALSE, 
                 multithread=TRUE)

MA4 <- mergeMulti(MA3, justConcatenate=TRUE)
MA5 <- makeSequenceTableMulti(MA4)
MA6 <- removeChimeraMulti(MA5)

lapply(mapReadsStratTab(MA6), print)

trackReadSorting(MA6)



