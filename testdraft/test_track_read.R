#Primers used in the arrays, primer pairs in single processin are part of this
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

### The UNIQUE BUG doesn't change with 
MA1 <- sortAmplicons(MA, filedir=tempfile())

## MA2 <- derepMulti(MA1, mc.cores=1)
## MA3 <- dadaMulti(MA2, selfConsist=TRUE, pool=FALSE, 
##                  multithread=TRUE)
## MA4 <- mergeMulti(MA3, justConcatenate=TRUE)
## MA5 <- makeSequenceTableMulti(MA4)
## MA6 <- removeChimeraMulti(MA5)

#################### SAMPLE CONFUSION #########################
## BAD2 <- derepMulti(MA1[c(6:4,1L,3L,2L), c(6:4,1L,3L,2L, 7L)], mc.cores=1)
BAD2 <- derepMulti(MA1, mc.cores=1)

### THIS LEADS TO the same SAMPLE CONFUSION
## BAD3 <- dadaMulti(BAD2[c(6:4,1L,3L,2L), c(6:4,1L,3L,2L, 7L)],
##                   selfConsist=TRUE, pool=FALSE, multithread=TRUE)

## resorting samples is enough to cause trouble
## BAD3 <- dadaMulti(BAD2[, c(6:4,1L,3L,2L, 7L)],
##                   selfConsist=TRUE, pool=FALSE, multithread=TRUE)

## resorting samples is enough to cause trouble... easy resort does
## the job
BAD3 <- dadaMulti(BAD2[, 7:1],
                   selfConsist=TRUE, pool=FALSE, multithread=TRUE)

MA3 <- dadaMulti(BAD2,
                   selfConsist=TRUE, pool=FALSE, multithread=TRUE)


## resorting samples is enough to cause trouble... easy resort does
## the job... BUT NOT subsetting without resorting!!!!
## BAD3 <- dadaMulti(BAD2[, 1:7],
##                  selfConsist=TRUE, pool=FALSE, multithread=TRUE)

### this resorting ampicons doesn't matter
## BAD3 <- dadaMulti(MA2[c(6:4,1L,3L,2L), ], selfConsist=TRUE, pool=FALSE, 
##                   multithread=TRUE)

## here it doesn't produce an error
## BAD4 <- mergeMulti(BAD3[c(6:4,1L,3L,2L), c(6:4,1L,3L,2L, 7L)],
##                    justConcatenate=TRUE)
BAD4 <- mergeMulti(BAD3, justConcatenate=TRUE)

MA4 <- mergeMulti(MA3, justConcatenate=TRUE)

## here it doesn't produce an error
## BAD5 <- makeSequenceTableMulti(BAD4[c(6:4,1L,3L,2L), c(6:4,1L,3L,2L, 7L)])
BAD5 <- makeSequenceTableMulti(BAD4)
MA5 <- makeSequenceTableMulti(MA4)

## here it doesn't produce an error
## BAD6 <- removeChimeraMulti(BAD5[c(6:4,1L,3L,2L), c(6:4,1L,3L,2L, 7L)])
BAD6 <- removeChimeraMulti(BAD5)
MA6 <- removeChimeraMulti(MA5)

################# EVALUATE #############
seqtabs <- getSequenceTableNoChime(MA6)
BADtabs <- getSequenceTableNoChime(BAD6)

SamSums <- lapply(seqtabs, rowSums)
BADSums <- lapply(BADtabs, rowSums)

samorder <- colnames(MA)

confusion <- lapply(names(SamSums), function(name) {
    df <- cbind(SamSums[[name]][samorder], BADSums[[name]][samorder])
    rownames(df) <- samorder
    colnames(df) <- c("Correct", "Shuffle")
    df[is.na(df)] <- 0
    df
})


### OKAY THIS IS ONE BUG found! FIX IT!!!!

#### Summary:

## If before the dada step samples are resorted we have a
## mess!!!######


########################################
### but the reads confusion must be in the sorting is 


trackReadSorting <- function (MA) { 
### somehow stupid that this is necessary (inherit form list) to
### improve (small TODO)
    readsFL <- lapply(seq_along(MA@PairedReadFileSet), function(i) {
        rawFiles <- MA@PairedReadFileSet@readsF[[i]]
        rawFiles <- rawFiles[file.exists(rawFiles)]
        readFastq(rawFiles)
    })
    snames <- colnames(MA)
    names(readsFL) <- snames
    sort_track <- lapply(snames, function (sampl) { 
        strat <- lapply(MA@stratifiedFiles, function(x) {
            grep(sampl, x@readsF, value=TRUE)
        })
        readsF_stratified <- readFastq(unlist(strat))
        IDsStrat <- ShortRead::id(readsF_stratified)
        readsF <- readsFL[[sampl]]        
        IDsRaw <- ShortRead::id(readsF)
        cbind(NrawReads = length(readsF),
              NUIDrawReads = length(unique(IDsRaw)),
              NstratReads = length(readsF_stratified),
              NUIDstratReads = length(unique(IDsStrat)),
              Problems = length(IDsStrat[!IDsStrat %in% IDsRaw]))
    })
    sort_track <- as.data.frame(do.call(rbind, sort_track))
    rownames(sort_track) <- snames
    sort_track$lessUraw <- sort_track$NrawReads - sort_track$NUIDrawReads
    sort_track$lessUstrat  <- sort_track$NstratReads - sort_track$NUIDstratReads
    sort_track
}


trackReadSorting(MA1)


### This shows another (small?) problem, some reads are sorted into
### multiple amplicons
lapply(MA1@stratifiedFiles, function (x){
    readsF_stratified <- readFastq(unlist(x@readsF))
    table(duplicated(ShortRead::id(readsF_stratified)))
})
## this shows that reads are still unique within on amplicon




