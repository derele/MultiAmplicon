##' get the unmatched sequences for a MultiAmplicon project
##'
##' This function computes the number of sequencing read pairs not
##' matched by primer pairs after running \code{\link{sortAmplicons}}
##' @title getUnmatched
##' @param MA A code{\link{MultiAmplicon}} object on which
##'     \code{link{sortAmplicons}} has been called previously.
##' @param tabulate.only Binary TRUE or FALSE indcating whether only
##'     the number of matching reads should be tabulated if set to
##'     false non-matched sequencing reads are written into outdir.
##' @param outdir output directory into which files for unmatched
##'     reads are written. Files are unmatched_F.fastq and
##'     unmatched_R.fastq for forward and reverse reads, respectively.
##' @return A table indicting the number of read pairs matched (TRUE)
##'     and unmatched (FALSE) by primer pairs.
##' @author Emanuel Heitlinger
##' @importFrom ShortRead readFastq id writeFastq "%in%"
##' @export
getUnmatched <- function(MA, tabulate.only=TRUE, outdir=tempdir()){
    if(!length(MA@stratifiedFiles)>0){
        stop("stratifiedFiles slot is empty, please run sortAmplicons on your MA object before you try to get the unmatched sequences")
    }
    foundFilesF <- unlist(lapply(MA@stratifiedFiles, "slot", "readsF"))
    ## Don't follow R, as F and R must be read pairs
    ## foundFilesR <- unlist(lapply(MA@stratifiedFiles, "slot", "readsR"))
    foundReadsF<- ShortRead::readFastq(foundFilesF)
    ##  foundReadsR<- readFastq(foundFilesR)
    searchedReadsF <- ShortRead::readFastq(MA@PairedReadFileSet@readsF)
   `%here%` <- ShortRead::`%in%` # needed to call binary operator from package
    is.found <- ShortRead::id(searchedReadsF)%here%ShortRead::id(foundReadsF)
    if(!tabulate.only){
        ## needed only for output
        searchedReadsR <- readFastq(MA@PairedReadFileSet@readsR)
        missedReadsF <- searchedReadsF[!is.found]
        missedReadsR <- searchedReadsR[!is.found]
        if(!dir.exists(outdir)) {
        cat("creating directory ", outdir, "for files of unmatched reads\n")
        dir.create(outdir, recursive=TRUE)
        } else {
            cat("using existing directory ", outdir)
            if(file.exists(paste0(outdir, "/unmatched_F.fastq"))){
                stop("file", paste0(outdir, "/unmatched_F.fastq"),
                     " already exists, delete it or provide another output directory")
            }
        }
        ShortRead::writeFastq(missedReadsF, file=paste0(outdir, "/unmatched_F.fastq"))
        ShortRead::writeFastq(missedReadsR, file=paste0(outdir, "/unmatched_R.fastq"))
        cat("witten", paste0(outdir, "/unmatched_F.fastq"), "and",
            paste0(outdir, "/unmatched_R.fastq"), "\n")
    }
    return(table("reads pairs matched primer pairs"=is.found))
}
