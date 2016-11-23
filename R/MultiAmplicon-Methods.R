################################################################################
##' .. content for description (no empty lines) ..
##'
##' .. content for details
##' 
##' @title sortAmplicons
##' @param MA 
##' @param n 
##' @param ... 
##' @return MultiAmplicon
##' @author Emanuel Heitlinger
##' @importFrom ShortRead FastqStreamer yield narrow
##' @importFrom Biostrings isMatchingStartingAt
setMethod("sortAmplicons", "MultiAmplicon", function(MA, n=1e6, ...){
    ## the data matrix of amplicons x samples stratified counts 
    data <- matrix(0,
                   nrow=length(MA@PrimerPairsSet),
                   ncol=length(MA@PairedReadFileSet))
    ## colnames are sample names taken from file names
    colnames(data) <- names(MA@PairedReadFileSet)
    ## rownames have to come from (matched) primers
    rownames(data) <- names(MA@PrimerPairsSet)
    tmppathF <- matrix("",
                       nrow=length(MA@PrimerPairsSet),
                       ncol=length(MA@PairedReadFileSet))
    tmppathR <- matrix("",
                       nrow=length(MA@PrimerPairsSet),
                       ncol=length(MA@PairedReadFileSet))
    readsF <- MA@PairedReadFileSet@readsF
    readsR <- MA@PairedReadFileSet@readsR
    ## add sample data and metadata in columns
    for(x in seq_along(readsF)) {
        f <- ShortRead::FastqStreamer(readsF[[x]], n = n)
        r <- ShortRead::FastqStreamer(readsR[[x]], n = n)
        tmpbaseF <- paste0(tempfile(), basename(readsF[[x]]))
        tmpbaseR <- paste0(tempfile(), basename(readsR[[x]]))
        ## request forward and reverse file simultaneously
         while(length(suppressWarnings(Ffq <- ShortRead::yield(f))) &&
               length(suppressWarnings(Rfq <- ShortRead::yield(r)))){
                   fM <- lapply(MA@PrimerPairsSet@.uniqueF, function(x){
                       as.vector(Biostrings::isMatchingStartingAt(x, sread(Ffq),
                                                      fixed=FALSE))
                   })
                   rM <- lapply(MA@PrimerPairsSet@.uniqueR, function(x){
                       as.vector(Biostrings::isMatchingStartingAt(x, sread(Rfq),
                                                      fixed=FALSE))
                   })
                   matches <- numeric(length=length(MA@PrimerPairsSet))
                   ## add primer pair data and metadata in rows 
                   for(y in seq_along(MA@PrimerPairsSet)) {
                       map.primerF <- MA@PrimerPairsSet@.mapF[[y]]
                       map.primerR <- MA@PrimerPairsSet@.mapR[[y]]
                       select <- fM[[map.primerF]] & rM[[map.primerR]]
                       tmppathF[y, x] <- paste0(tmpbaseF, names(MA@PrimerPairsSet)[[y]], ".fastq.gz")
                       tmppathR[y, x] <- paste0(tmpbaseR, names(MA@PrimerPairsSet)[[y]], ".fastq.gz")
                       lengthF <- length(foo@PrimerPairsSet@primerF[[y]])
                       lengthR <- length(foo@PrimerPairsSet@primerR[[y]])
                       F <- ShortRead::narrow(Ffq[select],
                                              lengthF, width(Ffq[select]))
                       R <- Rfq[select]
                       R <- ShortRead::narrow(Rfq[select],
                                              lengthR, width(Rfq[select]))
                       ShortRead::writeFastq(F, file=tmppathF[[y, x]], mode="a")
                       ShortRead::writeFastq(R, file=tmppathR[[y, x]], mode="a")
                       matches[y] <- length(select[select==TRUE])
                   }
                   ## need to add over the while loop because of fastq streaming 
                   data[, basename(readsF[[x]])] <-
                       data[, basename(readsF[[x]])] + matches
               }
        close(f)
        close(r)
        cat("\n finished sorting ", sum(data[, basename(readsF[[x]])]), " sequencing reads in",
            "\n ", readsF[[x]], " and \n ", readsR[[x]], "\n",
            " into ", sum(data[, basename(readsF[[x]])]>0), "amplicons \n" )
    }
    return(new("MultiAmplicon",
               PrimerPairsSet = MA@PrimerPairsSet,
               PairedReadFileSet = MA@PairedReadFileSet,
               rawCounts=data,
               FstratifiedFiles=tmppathF,
               RstratifiedFiles=tmppathR
               ))
})

