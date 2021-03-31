.concatenatePairedReadFileSets <- function(x, y){
    PairedReadFileSet(readsF=c(x@readsF, y@readsF),
                      readsR=c(x@readsR, y@readsR))
}


.concatenateSampleData <- function(x, y){
    bound <- rbind(x@sampleData, y@sampleData)
    new("sample_data", bound)
}

.concatenateSequenceTable <- function(x, y){
    st <- lapply(seq_along(x), function (i){
        bound <- rbind(x[[i]], y[[i]])
        ## ugly hack to remove dimnames in empty IMPROVE??
        if(all(dim(bound)==0)) {
            attr(bound, "dimnames") <- NULL
        }
        bound
    })
    names(st) <- names(x)
    st
}

.concatenateTaxonTable <- function(x, y){
    tt <- lapply(seq_along(x), function (i){
        bound <- rbind(x[[i]], y[[i]])
        if(!is.null(bound)){
            new("taxonomyTable", bound)
        } else NULL
    })
    names(tt) <- names(x)
    tt
}

##' Concatenate two MultiAmplicon objects
##'
##' Two MultiAmplicon objects are concatenated. Currently only
##' concatenation of different samples with the same amplicons is
##' implemented. 
##' 
##' @title concatenateMultiAmplicon
##' @param MA1 MultiAmplicon object that should be concatenated.
##' @param MA2 Second MultiAmplicon object to be concatenated with the
##'     first.
##' @param what Should either "samples" or "amplicons" be
##'     concatenated? Currently only "samples" are implemented.
##' @return A concatenated MultiAmplicon object
##' @export
##' @author Emanuel Heitlinger
concatenateMultiAmplicon <- function (MA1, MA2, what="samples") {
    if(!all(what%in%c("samples", "amplicons"))){
        stop("please indicate `what` you want to concatenate, `samples` or
`amplicons`?")
    }
    if(what%in%"amplicons"){
        stop("concatenation of amplicons is not yet implemented, plaese open and issue on github and ask me prioritize this feature if you need it!")
    }
    MultiAmplicon(PrimerPairsSet = getPrimerPairsSet(MA1),
                  PairedReadFileSet =.concatenatePairedReadFileSets(
                      getPairedReadFileSet(MA1), getPairedReadFileSet(MA2)),
                  rawCounts = cbind(getRawCounts(MA1), getRawCounts(MA2)),
                  stratifiedFilesF = cbind(
                      getStratifiedFilesF(MA1, dropEmpty=FALSE),
                      getStratifiedFilesF(MA2, dropEmpty=FALSE)),
                  stratifiedFilesR = cbind(
                      getStratifiedFilesR(MA1, dropEmpty=FALSE),
                      getStratifiedFilesR(MA2, dropEmpty=FALSE)),
                  sampleData = .concatenateSampleData(MA1, MA2),
                  derepF = cbind(getDerepF(MA1, dropEmpty=FALSE),
                                 getDerepF(MA2, dropEmpty=FALSE)),
                  derepR = cbind(getDerepR(MA1, dropEmpty=FALSE),
                                 getDerepR(MA2, dropEmpty=FALSE)),
                  dadaF = cbind(getDadaF(MA1, dropEmpty=FALSE),
                                 getDadaF(MA2, dropEmpty=FALSE)),
                  dadaR = cbind(getDadaR(MA1, dropEmpty=FALSE),
                                 getDadaR(MA2, dropEmpty=FALSE)),
                  mergers = cbind(getMergers(MA1, dropEmpty=FALSE),
                                  getMergers(MA2, dropEmpty=FALSE)),
                  .concatenateSequenceTable(getSequenceTable(MA1),
                                            getSequenceTable(MA2)),
                  .concatenateSequenceTable(getSequenceTableNoChime(MA1),
                                            getSequenceTableNoChime(MA2)),
                  .concatenateTaxonTable(getTaxonTable(MA1),
                                         getTaxonTable(MA2))
                  )
}


