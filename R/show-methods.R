##############################################################################
##' @param object A \code{PairedReadFileSet-class} object.
##' @rdname PairedReadFileSet-class
setMethod("show", "PairedReadFileSet", function(object) {
    cat("  A ", class(object),
        " instance of length ", length(object@readsF), 
        "\n", sep = "")
    if (length(object@readsF) != 0) {
        cat("\n", length(object@readsF), " Forward read files:",
            substr(paste(object@readsF, collapse=" "), 1, 180),
            "... \n")
        cat(length(object@readsR), " Reverse read files:",
            substr(paste(object@readsR, collapse=" "), 1, 180),
            "... \n")
    }
})

################################################################################' @param object A \code{PrimerPairsSet-class} object.
##' @rdname PrimerPairsSet-class
setMethod("show", "PrimerPairsSet", function(object){
    cat("  A ", class(object),
        " instance of length ", length(object@primerF), 
        "\n", sep = "")
    if (length(object@primerF) != 0) {
        cat("Forward:\n")
        cat(Biostrings:::.OBJECTStringSet.show_frame(object@primerF))
        cat("Reverse:\n")
        cat(Biostrings:::.XStringSet.show_frame(object@primerR))
    }
})


################################################################################' @param object A \code{MultiAmplicon-class} object.
##' @rdname MultiAmplicon-class 
setMethod("show", "MultiAmplicon", function(object){
    cat("  A ", class(object),
        " instance of dimensions ", length(object@PrimerPairsSet), " x ",
        length(object@PairedReadFileSet),
        "\n", sep = "")
    if (length(object@PairedReadFileSet) != 0 &&
        length(object@PrimerPairsSet) !=0) {
        cat("Forward Primers:\n")
        cat(Biostrings:::.XStringSet.show_frame(object@PrimerPairsSet@primerF))
        cat("Reverse Primers:\n")
        cat(Biostrings:::.XStringSet.show_frame(object@PrimerPairsSet@primerR))
        cat("\n", length(object@PairedReadFileSet), " Forward read files:",
            substr(paste(object@PairedReadFileSet@readsF, collapse=" "), 1, 180),
            "... \n")
        cat(length(object@PairedReadFileSet), " Reverse read files:",
            substr(paste(object@PairedReadFileSet@readsR, collapse=" "), 1, 180),
            "... \n")
    }
})

