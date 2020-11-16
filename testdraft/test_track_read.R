## test file for use whith testthat
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

MA1 <- sortAmplicons(MA, filedir=tempfile())
MA2 <- derepMulti(MA1, mc.cores=1)
MA3 <- dadaMulti(MA2, selfConsist=TRUE, pool=FALSE, 
                 multithread=TRUE)
MA4 <- mergeMulti(MA3, justConcatenate=TRUE)
MA5 <- makeSequenceTableMulti(MA4)
MA6 <- removeChimeraMulti(MA5)

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

## here it doesn't produce an error
## BAD5 <- makeSequenceTableMulti(BAD4[c(6:4,1L,3L,2L), c(6:4,1L,3L,2L, 7L)])
BAD5 <- makeSequenceTableMulti(BAD4)

## here it doesn't produce an error
## BAD6 <- removeChimeraMulti(BAD5[c(6:4,1L,3L,2L), c(6:4,1L,3L,2L, 7L)])
BAD6 <- removeChimeraMulti(BAD5)



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


########################################
### but the reads confusion must be in the sorting is 


trackReadSorting <- function (MA) { 
### somehow stupid that this is necessary (inherit form list) to
### improve (small TODO)
    readsFL <- lapply(seq_along(MA@PairedReadFileSet), function(i) {
        readFastq(MA@PairedReadFileSet@readsF[[i]])
    })
    snames <- colnames(MA)
    names(readsFL) <- snames
    sort_track <- lapply(snames, function (sampl) { 
        readsF <- readsFL[[sampl]]
        strat <- lapply(MA@stratifiedFiles, function(x) {
            grep(sampl, x@readsF, value=TRUE)
        })
        readsF_stratified <- readFastq(unlist(strat))
        cbind(length(readsF), 
              length(readsF_stratified),
              length(readsF_stratified[!id(readsF_stratified) %in% id(readsF)])
              )
    })
    sort_track <- as.data.frame(do.call(rbind, sort_track))
    colnames(sort_track) <- c("rawReads", "sortedReads", "Problems")
    rownames(sort_track) <- snames
    sort_track
}


