##' Sort different amplicons into a fully stratified samples x
##' amplicons structure based on primer matches.
##'
##' This functions uses Biostrings::isMatchingStartingAt to match
##' primer sequences at the first postion of forward and reverse
##' sequences. These sequences are written to temporary files to allow
##' processing via standard metabarcoding pipelines (i.e. dada2 using
##' fastq streaming).
##' 
##' @title sortAmplicons
##' @param MA MultiAmplicon-class object containing a set of paired
##'     end files and a primer-pairs set.
##' @param n parameter passed to the yield functions of package
##'     ShortRead. This controls the memory consumption during
##'     streaming. Lower values mean lower memory but longer
##'     processing time.
##' @param ... addtional parameter so be passed to
##'     Biostrings::isMatchingStartingAt. Be careful when using
##'     multiple starting positions or allowing error. This could lead
##'     to read pairs being assigned to multiple amplicons.
##' @return MultiAmplicon
##' @rdname sortAmplicons
##' @author Emanuel Heitlinger
##' @importFrom ShortRead FastqStreamer yield narrow sread
##' @importFrom Biostrings isMatchingStartingAt
##' @export
##' @aliases sortAmplicons, sortAmplicons-Method
setGeneric(name="sortAmplicons",
           def=function(MA, ...) {
               standardGeneric("sortAmplicons", ...)
           })
################################################################################
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
        tmpbaseF <- paste0(tempfile(), basename(readsF[[x]]))
        tmpbaseR <- paste0(tempfile(), basename(readsR[[x]]))
        if(!file.exists(readsF[[x]]) | !file.exists(readsR[[x]])){
            data[, names(readsF)[[x]]] <- 0
            tmppathF[,x] <- paste0(tmpbaseF, names(MA@PrimerPairsSet), ".fastq.gz")
            tmppathR[,x] <- paste0(tmpbaseR, names(MA@PrimerPairsSet), ".fastq.gz")
            break()
        }
        f <- ShortRead::FastqStreamer(readsF[[x]], n = n)
        r <- ShortRead::FastqStreamer(readsR[[x]], n = n)
        ## request forward and reverse file simultaneously
         while(length(suppressWarnings(Ffq <- ShortRead::yield(f))) &&
               length(suppressWarnings(Rfq <- ShortRead::yield(r)))){
                   fM <- lapply(MA@PrimerPairsSet@.uniqueF, function(x){
                       as.vector(Biostrings::isMatchingStartingAt(x,
                                                                  ShortRead::sread(Ffq),
                                                                  fixed=FALSE))
                   })
                   rM <- lapply(MA@PrimerPairsSet@.uniqueR, function(x){
                       as.vector(Biostrings::isMatchingStartingAt(x, ShortRead::sread(Rfq),
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
                       lengthF <- length(MA@PrimerPairsSet@primerF[[y]])
                       lengthR <- length(MA@PrimerPairsSet@primerR[[y]])
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
                   data[, names(readsF)[[x]]] <-
                       data[, names(readsF)[[x]]] + matches
               }
        close(f)
        close(r)
        cat("\n finished sorting", sum(data[, names(readsF)[[x]]]),
            "sequencing reads for", names(readsF)[[x]], "in",
            "\n ", readsF[[x]], " and \n ", readsR[[x]], "\n",
            " into ", sum(data[, names(readsF)[[x]]]>0), "amplicons \n" )
    }
    ## capture the warnings for non-existing files only if zeros are
    ## repored in the report of what was mapped
    stratifiedFiles <- lapply(seq_along(MA@PrimerPairsSet),
                              function(i) PairedReadFileSet(tmppathF[i,],
                                                            tmppathR[i,]))
    return(new("MultiAmplicon",
               PrimerPairsSet = MA@PrimerPairsSet,
               PairedReadFileSet = MA@PairedReadFileSet,
               rawCounts = data,
               stratifiedFiles = stratifiedFiles
               ))
})

