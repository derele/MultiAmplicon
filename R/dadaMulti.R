##' A wrapper around \code{\link[dada2]{dada}} from \code{dada2}
##'
##' The function runs \code{\link[dada2]{dada}} from
##' \code{\link[dada2]{dada}} to performe 'High resolution sample
##' inference from amplicon data' on multiple amplicons stored as
##' dereplicated sequences in a MultiAmplicon-class object
##' 
##' @title dadaMulti
##' @param MA MultiAmplicon-class object
##' @param ... additional parameters to be passed to
##'     \code{\link[dada2]{dada}} from \code{\link[dada2]{dada}}
##' @return MultiAmplicon object with dadaF and dadaR slots filled
##' @importFrom dada2 dada
##' @importFrom methods initialize new slot
##' @export
##' @author Emanuel Heitlinger
dadaMulti <- function(MA, ...){
    ## needs to be computed pair amplicon
    PPdada <- lapply(seq_along(MA@PrimerPairsSet), function (i){
       dF <- lapply(MA@derep[[i]], function (x) slot(x, "derepF"))
       dR <- lapply(MA@derep[[i]], function (x) slot(x,  "derepR"))
       cat("\n\namplicon", MA@PrimerPairsSet@names[i],
           "dada estimation of sequence variants from ",
            length(dF), " of ",
           length(MA@PairedReadFileSet), "possible sample files\n\n")
       if(length(dF)>0 && length(dR)>0){
           ## run functions for reverse and forward
           dadaF <- dada(dF, ...)
           ## make it a list of length 1 in case of only one sample,
           ## otherwise it is simplified and can't be handled
           if (class(dadaF)%in%"dada"){dadaF <- list(dadaF)}
           dadaR <- dada(dR, ...)
           ## make it a list in case of only one sample
           if (class(dadaR)%in%"dada"){dadaR <- list(dadaR)}
           Pdada <- PairedDada(dadaF = dadaF, dadaR = dadaR)
       } else {
           Pdada <- PairedDada()
           cat("skipping empty amplicon")
       }
       return(Pdada)
    })
    initialize(MA, dada=PPdada)
}

