### Useful functions
fixTableNames <- function (MA, slot="sequenceTableNoChime") {
    sampleNames <- colnames(MA)
    pattern <- paste(sampleNames, collapse="|")
    pattern <- paste0(".*(", pattern, ").*")
    stList <- slot(MA, slot)
    lapply(stList, function(ST) {
        newST <- ST
        rownames(newST) <- gsub(pattern, "\\1", rownames(ST))
        newST
    })
}


mapReadsStratTab <- function(MA) {
    getReadsBySample <- function(MA){
    sreads <- lapply(colnames(MA), function (sampl) { 
        strat <- lapply(MA@stratifiedFiles, function(x) {
            grep(sampl, x@readsF, value=TRUE)
        })
        readFastq(unlist(strat))
    })
    names(sreads) <- colnames(MA)
    sreads
    }
    RbyS <- getReadsBySample(MA)
    RbyS <- lapply(RbyS, ShortRead::sread)
    
    seqtabL <- getSequenceTableNoChime(MA)
    seqtabL <- seqtabL[unlist(lapply(seqtabL, function(x) all(dim(x)>0)))]

    SbyS <- lapply(colnames(MA), function(sampl){
        sbys <- lapply(seqtabL, function(ST) {
            sampl.here <- sampl[sampl%in%rownames(ST)]
            cn <- colnames(ST)[which(ST[sampl.here,]>0)]
            split <- strsplit(cn, "NNNNNNNNNN")
            unlist(lapply(split, "[[", 1))
            ## it might be necessary to also track the counts here
        })
        sbys[!unlist(lapply(sbys, is.null))]
    })
    names(SbyS) <- colnames(MA)

    SbyS <- SbyS[intersect(names(SbyS), names(RbyS))]
    RbyS <- RbyS[intersect(names(SbyS), names(RbyS))]

    sapply(names(SbyS), function (na) {
        unlist(SbyS[[na]]) %in% unique(as.vector(RbyS[[na]]))
    })
}



### Test read confusion in the sortin
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
        cbind(RawReads = length(readsF),
              UniqueRawReads = length(unique(IDsRaw)),
              SortedReads = length(readsF_stratified),
              UniqueSortedReads = length(unique(IDsStrat)),
              Problems = length(IDsStrat[!IDsStrat %in% IDsRaw]))
    })
    sort_track <- as.data.frame(do.call(rbind, sort_track))
    rownames(sort_track) <- snames
    sort_track$DuplicateRaw <- sort_track$RawReads -
        sort_track$UniqueRawReads
    sort_track$DuplicateSorted  <- sort_track$SortedReads -
        sort_track$UniqueSortedReads

    sort_track
}




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

MA6@sequenceTableNoChime <- fixTableNames(MA6)

lapply(mapReadsStratTab(MA6), all)

trackReadSorting(MA6)


### THIS analyses SAMPLE CONFUSION caused by resorting before dada
BAD3 <- dadaMulti(MA1[, c(6:4,1L,3L,2L, 7L, 10:8)],
                  selfConsist=TRUE, pool=FALSE, multithread=TRUE)

## in all other steps resorting does not produce an error
BAD4 <- mergeMulti(BAD3, justConcatenate=TRUE)
BAD5 <- makeSequenceTableMulti(BAD4)
BAD6 <- removeChimeraMulti(BAD5)

## BAD6@sequenceTableNoChime <- fixTableNames(BAD6)

lapply(mapReadsStratTab(BAD6), print)
trackReadSorting(BAD6)

################# EVALUATE sorting #############
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

confusion
### WRITE A TEST FOR IT!!


