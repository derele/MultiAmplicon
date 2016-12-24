##
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title dadaMulti
##' @param MA MultiAmplicon object
##' @param ... Parameters to be passed to dada2::dada 
##' @return MultiAmplicon object with dadaF and dadaR slots filled
##' @author Emanuel Heitlinger

dadaMulti <- function(MA,  ...) {
    MA@dadaF  <- lapply(MA@derepF, function(x){
        dada(x, err=NULL, selfConsist=selfConsist, ...)
    })
    MA@dadaR  <- lapply(MA@derepR, function(x){
        dada(x, err=NULL, selfConsist=selfConsist, ...)
    })
    return(MA)
}

