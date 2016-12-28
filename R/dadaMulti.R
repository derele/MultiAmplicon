##
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title dadaMulti
##' @param MA MultiAmplicon object
##' @param ... Parameters to be passed to dada2::dada 
##' @return MultiAmplicon object with dadaF and dadaR slots filled
##' @author Emanuel Heitlinger


dadaMulti <- function(MA, ...){
    dada <- lapply(seq_along(MA@PrimerPairsSet), function (i){
        cat("amplicon", rownames(MA)[i], "dada estimateion of sequence variants from ",
            length(MA@dada[[i]]@derepF), " of ",
            ncol(MA), "possible sample files\n")
        dadaF <- dada(MA@dada[i]@derepF, ...)
        dadaR <- dada(MA@dada[i]@derepR, ...)
        new("PairedDada",
            dadaF = dadaF,
            dadaR = dadaR,
            names = rownames(MA)[i])
    })
    new("MultiAmplicon",
        PrimerPairsSet = MA@PrimerPairsSet,
        PairedReadFileSet = MA@PairedReadFileSet,
        stratifiedFiles = MA@stratifiedFiles,
        derep = MA@derep,
        dada = new("PairedDadaSet", dada))
}

