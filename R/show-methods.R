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
    show(object@PrimerPairsSet)
    cat("\nContaining slot PairedReadFileSet: \n")
    show(object@PairedReadFileSet)
    cat("\nContaining slot stratifiedFiles of dimensions:",
        length(object@stratifiedFiles), "amplicons x ",
        paste(.replace_inf_range(lapply(object@stratifiedFiles, length)), collapse=" to "),
        " samples \n")
    cat("\nContaining slot sampleData:")
    show(object@sampleData)        
    cat("\nContaining slot derep of dimensions:",
        length(object@derep), "amplicons x ",
        paste(.replace_inf_range(lapply(object@derep, length)), collapse=" to "),
        " samples \n")
    cat("\nContaining slot dada of dimensions:",
        length(object@dada), "amplicons x ",
        paste(.replace_inf_range(lapply(object@dada, length)), collapse=" to "),
        " samples \n")
    cat("\nContaining slot mergers of dimensions:",
        length(object@mergers), "amplicons x ",
        paste(.replace_inf_range(lapply(object@mergers, length)),
              collapse=" to "),
        " samples ==>", 
        paste(.replace_inf_range(lapply(object@mergers, function (x) {
            sum(unlist(lapply(x, nrow)))
        })),
        collapse=" to "),
        "ASVs \n")
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
