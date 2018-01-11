##' Dereplicate sequences in fastq files
##'
##' An interface to \code{\link[dada2]{derepFastq}} which itself uses
##' \code{\link[ShortRead]{FastqStreamer}} for dereplicating amplicon
##' sequences from fastq or compressed fastq files.
##'
##' @title derepMulti
##' @param MA MultiAmplicon object to be dereplicated
##' @param ... arguments to be passed to \code{\link{derepFastq}}
##' @return MultiAmplicon object with derep slots (forward derepF and
##'     reverse derepR) filled
##' @importFrom dada2 derepFastq
##' @export
##' @author Emanuel Heitlinger
derepMulti <- function(MA, ...){
    PPderep <- lapply(seq_along(MA@PrimerPairsSet), function (i){
        cat("amplicon", rownames(MA)[i],
            "dereplicating for ",
            length(MA@stratifiedFiles[[i]]@readsF), " of ",
            length(MA@PairedReadFileSet@readsF), "possible sample files\n")
        derepF <- derepFastq(MA@stratifiedFiles[[i]]@readsF, ...)
        derepR <- derepFastq(MA@stratifiedFiles[[i]]@readsR, ...)
        ## make it a list even if only one sample was dereplicated
        if (class(derepF)%in%"derep") {derepF <- list(derepF)}
        if (class(derepR)%in%"derep") {derepR <- list(derepR)} 
        Pderep <- lapply(seq_along(derepF), function (w){
            new("PairedDerep",
                derepF = derepF[[w]],
                derepR = derepR[[w]])
        })
        return(Pderep)
    })
    initialize(MA, derep=PPderep)
}
