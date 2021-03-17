##############################################################################
##' @param object A \code{PairedReadFileSet-class} object.
##' @importFrom methods show
##' @rdname PairedReadFileSet-class
setMethod("show", "PairedReadFileSet", function(object) {
    cat("A", class(object), "instance of length", length(object@readsF), "\n")
    if (length(object@readsF) != 0) {
        cat("\n", length(object@readsF), " Forward read files:",
            substr(paste(object@readsF, collapse=" "), 1, 180),
            "... \n")
        cat(length(object@readsR), " Reverse read files:",
            substr(paste(object@readsR, collapse=" "), 1, 180),
            "... \n")
    }
})

################################################################################

##' @param object A \code{PrimerPairsSet-class} object.
##' @importFrom methods show
##' @rdname PrimerPairsSet-class
setMethod("show", "PrimerPairsSet", function(object){
    cat("A", class(object),
        "instance of length", length(object@primerF), "\n")
    if (length(object@primerF) != 0) {
        cat("Forward:\n")
        cat(Biostrings:::.XStringSet.show_frame(object@primerF))
        cat("Reverse:\n")
        cat(Biostrings:::.XStringSet.show_frame(object@primerR))
    }
})


################################################################################

##' @param object A \code{MultiAmplicon-class} object.
##' @importFrom methods show
##' @rdname MultiAmplicon-class 
setMethod("show", "MultiAmplicon", function(object){
    cat("A", class(object), "instance of dimensions", dim(object)[1], "",
        dim(object)[2],"\n")
    cat("\nContaining slot PrimerPairsSet: \n")
    show(getPrimerPairsSet(object))
    cat("\nContaining slot PairedReadFileSet: \n")
    show(getPairedReadFileSet(object))
    cat("\nContaining slot sampleData:")
    show(getSampleData(object))        
    cat("\nContaining slot rawCounts: \n")
    show(getRawCounts(object))
    cat("\nContaining slots stratifiedFilesF and stratifiedFilesR with dimensions:",
        dim(getStratifiedFilesF(MA1, dropEmpty=FALSE)),
        "\n of those", length(getStratifiedFilesF(object, dropEmpty=FALSE)),
        ",", length(getStratifiedFilesF(object, dropEmpty=TRUE)),
        "have files present\n")
    cat("\nContaining slots derepF and derepR of dimensions:",
        dim(getDerepF(object)), "\n")
    cat("\nContaining slots dadaF and dadaR of dimensions:",
        dim(getDadaF(object, dropEmpty=FALSE)),
    "\n of those", length(getDadaF(object, dropEmpty=FALSE)),
    ",", length(getDadaF(object, dropEmpty=TRUE)),
    "have dada objects present\n")
    cat("\nContaining slot mergers of dimensions:",
        dim(getMergers(object)),"\n")
    cat("\nContaining slot sequenceTable of dimensions:",
        length(object@sequenceTable), "amplicons  x ",
        paste(.replace_inf_range(lapply(object@sequenceTable, nrow)),
              collapse=" to "),
        " samples ==>",
        paste(.replace_inf_range(lapply(object@sequenceTable, ncol)),
              collapse=" to "),
        "ASVs \n")
    cat("\nContaining slot sequenceTableNoChime of dimensions:",
        length(object@sequenceTableNoChime), "amplicons  x ",
        paste(.replace_inf_range(lapply(object@sequenceTableNoChime, nrow)),
              collapse=" to "),
        " samples ==>",
        paste(.replace_inf_range(lapply(object@sequenceTableNoChime, ncol)),
              collapse=" to "),
        "ASVs \n")
    cat("\nContaining slot taxonTable of dimensions:",
        length(object@taxonTable), "amplicons  x ",
        "NA (without) samples ==>", 
        paste(.replace_inf_range(lapply(object@taxonTable, nrow)),
              collapse=" to "),
        " taxonomically annotated ASVs \n")
})

.replace_inf_range<- function(x){
    r <- suppressWarnings(range(x))
    r[!is.finite(r)] <- 0
    r
}
