##' Remove chimeric sequencing read pairs from a MultiAmplcion
##' object.
##'
##' This is a wrapper for the \code{\link[dada2]{removeBimeraDenovo}}
##' function of \code{dada2}. It works on an
##' \code{\link{MultiAmplicon-class}} object with the sequenceTable
##' slot filled. Use \code{\link{dadaMulti}},
##' \code{\link{derepMulti}}, \code{\link{mergeMulti}} and
##' \code{\link{sequenceTableMulti}} on a amplicon sorted (see
##' \code{\link{sortAmplicons}}) \code{\link{MultiAmplicon-class}}
##' object to preprocess your multi-marker data to this point.
##'
##' @title noChimeMulti
##' @param MA A \code{\link{MultiAmplicon-class}} object pre-processed
##'     to have a sequenceTableMulti populated
##' @param mc.cores integer number of cores to use for parallelization
##' @param ... paramter passed through to
##'     \code{\link[dada2]{removeBimeraDenovo}}
##' @return a code{\link{MulitAmplicon-class}} object with the
##'     sequenceTableNoChime filled
##' @importFrom dada2 removeBimeraDenovo
##' @importFrom parallel mclapply
##' @export
##' @author Emanuel Heitlinger
noChimeMulti <- function(MA, mc.cores, ...){
    sequenceTableNoChime <-
        mclapply(MA@sequenceTable, function (x) { 
            removeBimeraDenovo(x, ...)
        },
        mc.cores=mc.cores)
            ## fix me to get the correct rownames on those...
    initialize(MA, sequenceTableNoChime = sequenceTableNoChime)
}
