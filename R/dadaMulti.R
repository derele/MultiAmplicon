##' A wrapper around \link{\code{dada}} from \link{\code{dada2}}
##'
##' The function runs \link{\code{dada}} from \link{\code{dada2}} to
##' performe 'High resolution sample inference from amplicon data' on
##' multiple amplicons stored as dereplicated sequences in a
##' MultiAmplicon-class object
##' 
##' @title dadaMulti
##' @param MA MultiAmplicon-class object
##' @param ... additional parameters to be passed to
##'     \code{\link{dada}} from \code{\link{dada}}
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
           ## run both functions
           dadaF <- dada(dF, ...)
           dadaR <- dada(dR, ...)
           Pdada <- lapply(seq_along(dadaF), function (w){
               PairedDada(dadaF = dadaF[w], dadaR = dadaR[w])
           })
       } else {
           Pdada <- PairedDada()
           cat("skipping amplicon")
       }
       return(Pdada)
    })
    initialize(MA, dada=PPdada)
}

