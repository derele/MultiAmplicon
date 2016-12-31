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
##' @export
##' @author Emanuel Heitlinger
dadaMulti <- function(MA, ...){
    PPdada <- lapply(seq_along(MA@PrimerPairsSet), function (i){
        cat("amplicon", rownames(MA)[i], "dada estimation of sequence variants from ",
            length(MA@derep[[i]]@derepF), " of ",
            ncol(MA), "possible sample files\n")
        dadaF <- dada(MA@derep[[i]]@derepF, ...)
        dadaR <- dada(MA@derep[[i]]@derepR, ...)
        Pdada <- lapply(seq_along(dadaF), function (w){
            new("PairedDada",
                dadaF = dadaF[w],
                dadaR = dadaR[w])
        })
        return(Pdada)
    })
    return(PPdada)
}

