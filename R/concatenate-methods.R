.concatenatePairedReadFileSets <- function(x, y){
    PairedReadFileSet(readsF=c(x@readsF, y@readsF),
                      readsR=c(x@readsR, y@readsR))
}

.concatenateStratifiedFiles <- function(x, y){
    ## to avoid but when one of the MA is missing a sample
    along <- 1:max(length(x), length(y)) 
    cf <- lapply(along, function (i){
        .concatenatePairedReadFileSets(x[[i]], y[[i]])
    })
    names(cf) <- names(x)
    cf
}
.concatenateRawCounts <- function(x, y){
    cbind(getRawCounts(x), getRawCounts(y))
}

.concatenateSampleData <- function(x, y){
    bound <- rbind(x@sampleData, y@sampleData)
    new("sample_data", bound)
}

.concatenateDerep <- function(x, y){
    dr <- lapply(seq_along(x@derep), function (i){
        c(x@derep[[i]], y@derep[[i]])
    })
    names(dr) <- names(x@derep)
    dr
}

.concatenateDada <- function(x, y){
    FX <- getDadaF(x)
    FY <- getDadaF(y)
    RX <- getDadaR(x)
    RY <- getDadaR(y)
    PD <-  lapply(seq_along(FX), function (i){
        dF <- c(FX[[i]], FY[[i]])
        dR <- c(RX[[i]], RY[[i]])
        new("PairedDada", dadaF=dF, dadaR=dR)
    })
    names(PD) <- names(x@dada)
    PD
}

.concatenateMergers <- function(x, y){
    mg <- lapply(seq_along(x@mergers), function (i){
        c(x@mergers[[i]], y@mergers[[i]])
    })
    names(mg) <- names(x@mergers)
    mg
}

.concatenateSequenceTable <- function(x, y){
    st <- lapply(seq_along(x), function (i){
        rbind(x[[i]], y[[i]])
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
    MultiAmplicon(MA1@PrimerPairsSet,
                  .concatenatePairedReadFileSets(MA1@PairedReadFileSet,
                                                 MA2@PairedReadFileSet),
                  .concatenateRawCounts(MA1, MA2),
                  .concatenateStratifiedFiles(MA1@stratifiedFiles,
                                              MA2@stratifiedFiles),
                  .concatenateSampleData(MA1, MA2),
                  .concatenateDerep(MA1, MA2),
                  .concatenateDada(MA1, MA2),
                  .concatenateMergers(MA1, MA2),
                  .concatenateSequenceTable(MA1@sequenceTable,
                                            MA2@sequenceTable),
                  .concatenateSequenceTable(MA1@sequenceTableNoChime,
                                            MA2@sequenceTableNoChime),
                  .concatenateTaxonTable(MA1@taxonTable,
                                         MA2@taxonTable)
                  )
}


