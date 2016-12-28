##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param MA 
##' @param ... 
##' @return 
##' @author Emanuel Heitlinger
sequenceTableMulti <- function(MA, ...){
    sequenceTable <- lapply(seq_along(MA@mergers), function (i){
        ST <- makeSequenceTable(MA@mergers[i], ...)
        ## fix me to get the correct rownames on those...
        rownames(ST) <- rownames(MA[i])
        
    })
    initialize(MA, sequenceTable = sequenceTable)
}
