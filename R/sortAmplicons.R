##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title sortAmplicons
##' @param MA 
##' @param ... 
##' @return MultiAmplicon
##' @author Emanuel Heitlinger

## setGeneric(name="sortAmplicons",
##            def=function(MA, ...) {
##                standardGeneric("sortAmplicons", ...)
##            })

setMethod("sortAmplicons", "MultiAmplicon", function(MA, n=1e6, ...){
    data <- matrix(0,
                   nrow=length(MA@PrimerPairsSet),
                   ncol=length(MA@PairedReadFileSet))
    ## colnames are sample names taken from file names
    colnames(data) <- names(MA@PairedReadFileSet)
    ## rownames have to come from (matched) primers
    rownames(data) <- names(MA@PrimerPairsSet)
    readsF <- MA@PairedReadFileSet@readsF
    readsR <- MA@PairedReadFileSet@readsR
    for(i in seq_along(readsF)) {
        f <- FastqStreamer(readsF[[i]], n = n)
        r <- FastqStreamer(readsR[[i]], n = n)
        tmpbaseF <- paste0(tempfile(), basename(readsF[[i]]))
        tmpbaseR <- paste0(tempfile(), basename(readsR[[i]]))
        ## request forward and reverse file simultaneously
         while(length(suppressWarnings(Ffq <- yield(f))) &&
               length(suppressWarnings(Rfq <- yield(r)))){
                   fM <- lapply(MA@PrimerPairsSet@.uniqueF, function(x){
                       as.vector(isMatchingStartingAt(x, sread(Ffq),
                                                      fixed=FALSE))
                   })
                   rM <- lapply(MA@PrimerPairsSet@.uniqueR, function(x){
                       as.vector(isMatchingStartingAt(x, sread(Rfq),
                                                      fixed=FALSE))
                   })
                   matches <- sapply(seq_along(MA@PrimerPairsSet), function(i){
                       map.primerF <- MA@PrimerPairsSet@.mapF[[i]]
                       map.primerR <- MA@PrimerPairsSet@.mapR[[i]]
                       select <- fM[[map.primerF]] & rM[[map.primerR]]
                       tmppathF <- paste0(tmpbaseF, names(MA@PrimerPairsSet)[[i]], ".fastq.gz")
                       tmppathR <- paste0(tmpbaseR, names(MA@PrimerPairsSet)[[i]], ".fastq.gz")
                       lengthF <- length(foo@PrimerPairsSet@primerF[[i]])
                       lengthR <- length(foo@PrimerPairsSet@primerR[[i]])
                       F <- narrow(Ffq[select],
                                   lengthF, width(Ffq[select]))
                       R <- narrow(Rfq[select],
                                   lengthR, width(Rfq[select]))
                       writeFastq(F, file=tmppathF, mode="a")
                       writeFastq(R, file=tmppathR, mode="a")
                       number.matches <- length(select[select==TRUE])
                       return(number.matches)
                   })
                   ## need to add over the while loop because of fastq streaming 
                   data[, basename(readsF[[i]])] <-
                       data[, basename(readsF[[i]])] + matches
               }
        close(f)
        close(r)
        cat("\n finished sorting amplicons in \n ", readsF[[i]], " and \n ", readsR[[i]], "\n")
    }
    return(new("MultiAmplicon",
               PrimerPairsSet = MA@PrimerPairsSet,
               PairedReadFileSet = MA@PairedReadFileSet,
               rawCounts=data))
})
