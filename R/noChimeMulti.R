noChimeMulti <- function(MA, mc.cores, ...){
    sequenceTableNoChime <-
        mclapply(MA@sequenceTable, function (x) { 
            removeBimeraDenovo(x, ...)
        },
        mc.cores=mc.cores)
            ## fix me to get the correct rownames on those...
    initialize(MA, sequenceTableNoChime = sequenceTableNoChime)
}
