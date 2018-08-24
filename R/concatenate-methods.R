.concatenatePairedReadFileSets <- function(x, y){
    PairedReadFileSet(readsF=c(x@readsF, y@readsF),
                      readsR=c(x@readsR, y@readsR))
}

.concatenateStratifiedFiles <- function(x, y) 
lapply(seq_along(x), function (i){
    .concatenatePairedReadFileSets(x[[i]], y[[i]])
})

.concatenateRawCounts <- function(x, y){
    cbind(x@rawCounts, y@rawCounts)
}

.concatenateDerep <- function(x, y){
    lapply(seq_along(x@derep), function (i){
        c(x@derep[[i]], y@derep[[i]])
    })
}

.concatenateDada <- function(x, y){
    FX <- getDadaF(x)
    FY <- getDadaF(y)
    RX <- getDadaR(x)
    RY <- getDadaR(y)
    PD <-  lapply(seq_along(FX), function (i){
        dF <- c(FX[[i]], FY[[i]])
        dR <- c(RX[[i]], RY[[i]])
        new("PairedDada", dadaF=dF, dadaR=dR)
    })
    names(PD) <- names(x@dada)
    PD
}

.concatenateMergers <- function(x, y){
    lapply(seq_along(x@mergers), function (i){
        c(x@mergers[[i]], y@mergers[[i]])
    })
}

.concatenateSequenceTable <- function(x, y){
    lapply(seq_along(x), function (i){
        rbind(x[[i]], y[[i]])
    })
}



concatenateDadaMulti <- function (MA1, MA2, what="samples") {
    if(!all(what%in%c("samples", "amplicons"))){
        stop("please indicate `what` you want to concatenate, `samples` or
`amplicons`?")
    }
    MultiAmplicon(MA1@PrimerPairsSet,
                  .concatenatePairedReadFileSets(MA1@PairedReadFileSet,
                                                 MA2@PairedReadFileSet),
                  .concatenateRawCounts(MA1, MA2),
                  .concatenateStratifiedFiles(MA1@stratifiedFiles,
                                              MA2@stratifiedFiles),
                  .concatenateDerep(MA1, MA2),
                  .concatenateDada(MA1, MA2),
                  .concatenateMergers(MA1, MA2),
                  .concatenateSequenceTable(MA1@sequenceTable,
                                            MA2@sequenceTable),
                  .concatenateSequenceTable(MA1@sequenceTableNoChime,
                                            MA2@sequenceTableNoChime),
                  .concatenateSequenceTable(MA1@sequenceTableFilled,
                                            MA2@sequenceTableFilled)
                  )
}


