## name it when putting the below in the package:
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title derepMulti
##' @param MA MultiAmplicon object to be dereplicated
##' @param ... arguments to be passed to dada2::derepFastq
##' @return MultiAmplicon object with derep slots filled
##' @author Emanuel Heitlinger


derepMulti <- function(MA, ...){
    derepEmpty <- function(x, name, ...){
        file.present <- which(file.info(x)$size>21)
        cat("amplicon", name, "dereplicating for ",
            length(file.present), " of ", length(x), "possible sample files\n")
        derepFastq(x[file.present], ...)
    }
    MA@derepF <- lapply(seq_along(MA@PrimerPairsSet), function (i){
        derepEmpty(MA@FstratifiedFiles[i,], names(MA@PrimerPairsSet)[i])
    })
    MA@derepR <- lapply(seq_along(MA@PrimerPairsSet), function (i){
        derepEmpty(MA@RstratifiedFiles[i,], names(MA@PrimerPairsSet)[i])
    })
    return(MA)
}
