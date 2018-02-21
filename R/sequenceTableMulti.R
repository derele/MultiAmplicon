##' Create a sequence table inside a MultiAmplcion object.
##'
##' This is a wrapper for the \code{\link[dada2]{makeSequenceTable}}
##' function of \code{dada2}. It works on an
##' \code{\link{MultiAmplicon-class}} object with the mergers slot
##' filled. Use \code{\link{dadaMulti}} and \code{\link{derepMulti}}
##' on an amplicon sorted (see \code{\link{sortAmplicons}})
##' \code{\link{MultiAmplicon-class}} object to preprocess your
##' multi-marker data to this point.
##'
##' @title sequenceTableMulti 
##'
##' @param MA \code{\link{MultiAmplicon-class}} object with mergers
##'     slot filled
##' @param ... additional parameters to be passed to the
##'     \code{\link[dada2]{makeSequenceTable}} function of
##'     \code{dada2}
##' @return A MultiAmplicon-class object with the sequenceTable slot
##'     filled
##' @importFrom dada2 makeSequenceTable
##' @export
##' @author Emanuel Heitlinger
sequenceTableMulti <- function(MA, ...){
    sequenceTable <- lapply(seq_along(MA@mergers), function (i){
        makeSequenceTable(MA@mergers[[i]], ...)
    })
    names(sequenceTable) <- MA@PrimerPairsSet@names
    initialize(MA, sequenceTable = sequenceTable)
}
