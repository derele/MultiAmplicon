##' Dereplicate sequences in fastq files
##'
##' An interface to \code{\link[dada2]{derepFastq}} which itself uses
##' \code{\link[ShortRead]{FastqStreamer}} for dereplicating amplicon
##' sequences from fastq or compressed fastq files.
##'
##' @title derepMulti
##' @param MA MultiAmplicon object to be dereplicated
##' @param keep.single.singlets logical argument indicationg whether a
##'     derep slot should be filled with a single sequence recovered
##'     in only one sample. Defaults to FALSE meaning that such empty
##'     derep objects are created in such a case. Keeping such derep
##'     objects (setting this to TRUE) might result in downstream
##'     problems (Errors) in dada inference.
##' @param mc.cores number or compute cores for parallel processing
##' @param ... arguments to be passed to \code{\link{derepFastq}}
##' @return MultiAmplicon object with derep slots (forward derepF and
##'     reverse derepR) filled
##' @importFrom dada2 derepFastq
##' @importFrom parallel mclapply
##' @export
##' @author Emanuel Heitlinger
derepMulti <- function(MA, mc.cores=getOption("mc.cores", 2L),
                       keep.single.singlets=FALSE, ...){
    PPderep <- mclapply(seq_along(MA@PrimerPairsSet), function (i){
        cat("amplicon", MA@PrimerPairsSet@names[i],
            "dereplicating for ",
            length(MA@stratifiedFiles[[i]]@readsF), " of ",
            length(MA@PairedReadFileSet@readsF), "possible sample files\n")
        derepF <- derepFastq(MA@stratifiedFiles[[i]]@readsF, ...)
        derepR <- derepFastq(MA@stratifiedFiles[[i]]@readsR, ...)
        ## make it a list even if only one sample was dereplicated
        if (class(derepF)%in%"derep") {
            ## the same must be true for the revers then to keep them
            ## in the same lenghth
            if(length(derepF$map)==1 && keep.single.singlets==FALSE){
                derepF <- list()
                derepR <- list()
                cat("\nproducing empty derep object for amplicon",
                    MA@PrimerPairsSet@names[i],
                    "as only one sequence is reported for one sample ",
                    "set keep.single.singlets to TRUE to change this behaviour,",
                    "but be warned that this may lead to downstream errors\n\n" )
            } else{
                derepF <- list(derepF)
                derepR <- list(derepR)
            }
        }
        Pderep <- lapply(seq_along(derepF), function (w){
            new("PairedDerep",
                derepF = derepF[[w]],
                derepR = derepR[[w]])
        })
        return(Pderep)
    }, mc.cores=mc.cores)
    names(PPderep) <- MA@PrimerPairsSet@names
    initialize(MA, derep=PPderep)
}
