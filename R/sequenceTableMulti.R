## @export
sequenceTableMulti <- function(MA, ...){
    sequenceTable <- lapply(seq_along(MA@mergers), function (i){
        makeSequenceTable(MA@mergers[[i]], ...)
    })
    initialize(MA, sequenceTable = sequenceTable)
}
