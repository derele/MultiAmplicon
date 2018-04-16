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
    cat("A", class(object), "instance of dimensions", dim(object)[1], "object",
        dim(object)[2],"\n")
    cat("\nContaining slot PrimerPairsSet: \n")
    show(object@PrimerPairsSet)
    cat("\nContaining slot PairedReadFileSet: \n")
    show(object@PairedReadFileSet)
    cat("\nContaining slot stratifiedFiles of length:",
        length(object@stratifiedFiles), "\n")
    cat("\nContaining slot derep of length:",
        length(object@derep), "\n")
    cat("\nContaining slot dada of length:",
        length(object@dada), "\n")
    cat("\nContaining slot mergers of length:",
        length(object@mergers), "\n")
    cat("\nContaining slot sequenceTable of length:",
        length(object@sequenceTable), "\n")
    cat("\nContaining slot sequenceTableNoChime of length:",
        length(object@sequenceTableNoChime), "\n")
    cat("\nContaining slot sequenceTableFilled of length:",
        length(object@sequenceTableFilled), "\n")
})

