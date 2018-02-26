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
##' @param ... additional arguments to be passed to the
##'     \code{\link[dada2]{mergePairs}} function of \code{dada2}, all
##'     arguments to the function can be given as a vector of the same
##'     length as the number of primer pairs in the MultiAmplicon
##'     object, allowing to specify e.g. justConcatenate to be set
##'     TRUE for only some of the amplicons, or to specify different
##'     minOverlap for each amplicon. If a shorter vector is given it
##'     will be recycled to match the number of ampicons.
##' @return A MultiAmplicon-class object with the mergers slot filled
##' @importFrom dada2 mergePairs
##' @export
##' @author Emanuel Heitlinger
mergeMulti <- function(MA, ...){
    args <- list(...)
    lapply(args, function (x) {
        if(nrow(MA)%%length(x)>0){
            stop("argument of length ", length(x),
                 " can't be recycled to number of amplicons (",
                 nrow(MA), ")\n")
        }
    })
    exp.args <- lapply(args, function(x) {
        rep(x, times=nrow(MA)/length(x))
    })
    mergers <- lapply(seq_along(MA@PrimerPairsSet), function (i){     
        if(length(MA@dada[[i]]) > 0){
            daF <- slot(MA@dada[[i]], "dadaF")
            daR <- slot(MA@dada[[i]], "dadaR")
            deF <- lapply(MA@derep[[i]], "slot", "derepF")
            deR <- lapply(MA@derep[[i]], "slot", "derepR")
            cat("merging sequences from " , length(MA@dada[[i]]),
                "samples for amplicon ",
                MA@PrimerPairsSet@names[[i]], "\n")
            args.here <- lapply(exp.args, "[", i)
            print.args.here <- paste(names(args.here), unlist(args.here),
                                     sep="=")
            cat("calling merge Pairs with parameters",
                print.args.here, "\n")
            MP <- do.call(mergePairs,
                          c(list(dadaF=daF, derepF=deF,
                                 dadaR=daR, derepR=deR), args.here))
            ## correct the case of one sample / amplicon 
            if(class(MP)%in%"data.frame"){MP <- list(MP)}
            cat("DONE\n\n")
            return(MP)} else{
                          cat("skipping empty amplicon (sequences for" ,
                              length(MA@dada[[i]]), "samples)  ",
                              MA@PrimerPairsSet@names[[i]], "\n\n")
                          return(list())}
    })
    names(mergers) <- MA@PrimerPairsSet@names
    initialize(MA, mergers=mergers)
}
