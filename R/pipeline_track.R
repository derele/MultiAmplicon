##' Get a sum of uniques from an input object.
##'
##' This function calculates the number of unique sequences, ASVs or
##' merged sequences from the relevant slot in a MultiAmplicon object
##' 
##' @title getCounts
##' @param MaMat A matrix or vector of derep, dada or merger objects
##' @return a matrix, list or vector of counts, depending on the class
##'     of the input
##' @export
##' @author Emanuel Heitlinger
getCounts <- function(MaMat, what="input"){
    if (!what%in%c("uniques", "input")){
        stop("please specify what should be returned: the number of 'uniques' or 'input'")
    }
    if(is.matrix(MaMat) && !is.null(dim(MaMat)) && all(dim(MaMat)>0)) {
        apply(MaMat, 2, function (x) {
            sapply(x, .sumThingsOrNull, what)
        })
    } else if(is.list(MaMat)||(is.vector(MaMat)&&is.character(MaMat))) {
        t(sapply(MaMat, .sumThingsOrNull, what))
    } else {
        stop(paste("nothing calculated, supply a list or a matrix of fastq file names",
                   "derep, dada or merger objects"))
    }
}

### ##  TODO: this doesn't work yet!


##' Get a sum of the uniques vector or matrix or 0 if empty
##'
##' This function uses getUniques from dada2 to extract (potentially)
##' unique entries (ASVs, sequences or merged sequences) from a vector
##' of suitable objects.
##' 
##' @title .sumThingsOrNull
##' @param y a vector of dada, derep, merger objects from dada2
##' @return Numeric vector for the count of entries (uniques or
##'     sequences, ASVs)
##' @importFrom dada2 getUniques
##' @importFrom ShortRead readFastq
##' @author Emanuel Heitlinger
.sumThingsOrNull <- function(y, what){
    calcu <- function(x, what) {
       switch(what,
              input = sum(x),
              uniques = length(x[x>0])) ## >0 or the table objects
    }
    ## derep, dada or merger objects (mergers have no proper class)
    if(class(y) %in% c("derep", "dada")||
       (is.data.frame(y) &&
        all(c("sequence", "abundance") %in% 
            colnames(y)))){ 
        calcu(getUniques(y), what)
    }
    ## or stratified files, if they exist (also no proper class)
    else if(is.character(y) && any(file.exists(y))) {
        ### if they contain data
        if(length(readFastq(y))>0) {
            calcu(getUniques(y), what)
        } else { ## they exist but contain no data
            return(0)
        }
    }
    ## those are sequence tables
    else if (is.matrix(y)){
        apply(y, 1, calcu, what)
    }
    ### if those are NULL objects
    else if (is.null(y)) {
        return(0)
    }
    ###  if all stratified files don't exist 
    else if((is.character(y) && all(!file.exists(y)))){ 
        return(0)
    } else {
        stop("incompatible object to count uniques")
    }
}

### a SHITLOAD of good tests possible with this now!!
## getRawCounts(MA6) == getCounts(getStratifiedFilesF(MA6, dropEmpty=FALSE), "input") ==
## getRawCounts(MA6) == getCounts(getDerepF(MA6, dropEmpty=FALSE), "input")



## Summary function for pipeline ------------------------------------

##' Obtain summary data for amplicons run through the MultiAmplicon pipeline.
##'
##' Get statistics on the number of samples (with read data), the
##' number of unique sequence variants and the number of reads left
##' after processing of amplicons in the MultiAmplicon pipeline. In
##' some steps of the pipeline dada2 performs quality filtering
##' excluding non-credible sequence variants.
##' 
##' @title getPipelineSummary
##' @param MA MultiAmplicon object with all slots filled for tracking.
##' @return a data.frame of sample, unique sequences and sequencing
##'     reads numbers per amplicon.
##' @importFrom plyr revalue
##' @export
##' @author Emanuel Heitlinger
getPipelineSummary <- function(MA){
    slots <- c("stratifiedFiles", "derep", "dada", "mergers",
               "sequenceTable", "sequenceTableNoChime")
    slotFilled <- unlist(sapply(slots, function (x) length(slot(MA, x)))>0)




### helper functions
    track.l <- lapply(seq_along(getDadaF(MA)), function (i) {
        samples <- list(
            sorted=length(getRawCounts(MA[i, ])[getRawCounts(MA[i, ])>0]),
            ##  ## derep is currently defunct 
            ##  derep=length(getDerepF(MA[i, ])),
            denoised = if(slotFilled["dada"]) {
                           length(getDadaF(MA[i,]))
                       } else{0},
            merged = if(slotFilled["mergers"]) {
                         length(getMergers(MA[i,]))
                     } else{0},
            tabulated = if(slotFilled["sequenceTable"]) {
                            nrow(getSequenceTable(MA[i,]))
                        } else{0},
            noChime= if(slotFilled["sequenceTableNoChime"]) {
                         nrow(getSequenceTableNoChime(MA[i,]))
                     } else{0}
        )
        uniques <- list(
            ##  ## derep is currently defunct 
            ##  derep=sum(getU(getDerepF(MA[i, ]))),
            denoised = if(slotFilled["dada"]) {
                           sum(getU(getDadaF(MA[i, ])))
                       } else{0},
            merged = if(slotFilled["mergers"]) {
                         sum(getU(getMergers(MA[i, ])))
                     } else{0},
            tabulated = if(slotFilled["sequenceTable"]) {
                            ncol(getSequenceTable(MA[i, ]))
                        } else{0},
            noChime = if(slotFilled["sequenceTableNoChime"]) {
                          ncol(getSequenceTableNoChime(MA[i, ]))
                      } else{0}
        )
        reads <- list(
            sorted=sum(getRawCounts(MA[i, ])),
            ##  ## derep is currently defunct 
            ##  derep=sum(getN(getDerepF(MA[i, ]))),
            denoised = if(slotFilled["dada"]) {
                           sum(getN(getDadaF(MA[i, ])))
                       } else{0},
            merged = if(slotFilled["mergers"]) {
                         sum(getN(getMergers(MA[i, ])))
                     } else {0}, 
            tabulated = if(slotFilled["sequenceTable"]) {
                            sum(getSequenceTable(MA[i, ]))
                        } else{0},
            noChime= if(slotFilled["sequenceTableNoChime"]) {
                         sum(getSequenceTableNoChime(MA[i, ]))
                     } else {0}
        )
        list(samples, uniques, reads)
    })
    track <- reshape::melt(track.l)
    track$L2 <- revalue(as.factor(track$L2), c("1"="samples",
                                               "2"="uniques", "3"="reads"))
    track$L3 <- factor(track$L3, levels = c("sorted", ## "derep",
                                            "denoised", "merged",
                                            "tabulated", "noChime"))
    names(track) <- c("value", "pipeStep", "what", "primer")
    return(track)
}
