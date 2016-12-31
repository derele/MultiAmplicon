## @export
sequenceTableMulti <- function(MA, ...){
    sequenceTable <- lapply(seq_along(MA@mergers), function (i){
        ST <- makeSequenceTable(MA@mergers[[i]], ...)
        ## fix me to get the correct rownames on those...
##        rownames(ST) <- rownames(MA)[i]
        
    })
    initialize(MA, sequenceTable = sequenceTable)
}
