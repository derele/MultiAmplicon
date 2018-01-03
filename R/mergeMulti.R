##' merge denoised pairs of forward and reverse reads inside an
##' MultiAmplcion object.
##'
##' This is a wrapper for \code{\link[dada2]{mergePairs}} from
##' \code{dada2}. It works on an \code{\link{MultiAmplicon-class}}
##' object with derep and dada slots filled. Use
##' \code{\link{dadaMulti}} and \code{\link{derepMulti}} on a amplicon
##' sorted (see \code{\link{sortAmplicons}})
##' \code{\link{MultiAmplicon-class}} object to preprocess your
##' multi-marker data to this point.
##'
##' @title mergeMulti
##' @param MA \code{\link{MultiAmplicon-class}} object with derep and
##'     dada slots filled.
##' @param ... additional parameters to be passed to the
##'     \code{\link[dada2]{mergePairs}} function of \code{dada2}
##' @return A MultiAmplicon-class object with the mergers slot filled
##' @importFrom dada2 mergePairs
##' @export
##' @author Emanuel Heitlinger
mergeMulti <- function(MA, ...){
    mergers <- lapply(seq_along(MA@PrimerPairsSet), function (i){     
        if(length(MA@dada[[i]]) > 0){
            daF <- unlist(lapply(MA@dada[[i]], "slot", "dadaF"), recursive=FALSE)
            deF <- lapply(MA@derep[[i]], "slot", "derepF")
            daR <- unlist(lapply(MA@dada[[i]], "slot", "dadaR"), recursive=FALSE)
            deR <- lapply(MA@derep[[i]], "slot", "derepR")
            cat("merging sequences from " , length(MA@dada[[i]]),
                "samples for amplicon ",
                MA@PrimerPairsSet@names[[i]], "\n")
            MP <- mergePairs(daF, deF, daR, deR, ...)
            all.samples <- colnames(MA@rawCounts)[
                MA@rawCounts[i, ]>0]
            names(MP) <- all.samples
            cat("DONE\n\n")
            return(MP)} else{
                          cat("skipping empty amplicon (sequences for" ,
                              length(MA@dada[[i]]), "samples)  ",
                              MA@PrimerPairsSet@names[[i]], "\n\n")
                          return(list())}
    })
    initialize(MA, mergers=mergers)
}
