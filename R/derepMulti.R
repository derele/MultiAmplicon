##' Dereplicate sequences in fastq files
##'
##' An interface to \code{\link{derepFastq}} which itself uses
##' \code{\link{FastqStreamer}} for dereplicating amplicon sequences
##' from fastq or compressed fastq files. 
##'
##' @title derepMulti
##' @param MA MultiAmplicon object to be dereplicated
##' @param ... arguments to be passed to \code{\link{derepFastq}}
##' @return MultiAmplicon object with derep slots (forward derepF and
##'     reverse derepR) filled
##' @export
##' @author Emanuel Heitlinger

derepMulti <- function(MA, ...){
    derep <- lapply(seq_along(MA@PrimerPairsSet), function (i){
        cat("amplicon", rownames(MA)[i], "dereplicating for ",
            length(MA@stratifiedFiles[[i]]@readsF), " of ",
            ncol(MA), "possible sample files\n")
        derepF <- derepFastq(MA@stratifiedFiles[[i]]@readsF, ...)
        derepR <- derepFastq(MA@stratifiedFiles[[i]]@readsR, ...)
        new("PairedDerep",
            derepF = derepF,
            derepR = derepR)
    })
    initialize(MA,
               derep = new("PairedDerepSet", derep))
}
