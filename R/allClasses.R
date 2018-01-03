##' A class representing file paths to paired end reads.
##'
##' Two character vectors of the same length specifying file names of
##' paired end reads can be stored in this class. Filenames a checked
##' for their existence. Usually these are the filnames of quality
##' filtered fastq files already stratified into samples (one file
##' pair for each sample). 
##'
##' @slot readF A charcter vector specifying the file paths to files
##'     containing forward (sometimes called R1) sequencing
##'     reads. This vector can be named to store short short-name of
##'     samples.
##'
##' @slot readR A character vector specifiying the file path to a file
##'     containing reverse (sometimes called R2) sequencing
##'     reads. This vector can be named to store short short-name of
##'     samples. Names of readF and readR should be identical in this
##'     case
##'
##' @usage
##' ## Constructors:
##' PairedReadFileSet(readsF, readsR)
##'
##' @param readsF The path and filenames containing forward (R1) reads
##' @param readsR The path and filenames containing reverse (R2) reads
##' 
##' @return PairedReadFileSet
##' @author Emanuel Heitlinger
##' @export PairedReadFileSet
##'
setClass("PairedReadFileSet",
         slots = c(readsF="character", readsR="character", names="character"),
         contains = c(readsF="character", readsR="character", names="character"),
         validity=function(object) {
             all.files <- c(object@readsF, object@readsR)
             missing.file <- !file.exists(all.files)
             if (any(missing.file)){
                 warning(paste0("\nfile ", all.files[missing.file], " does not exist on your system"))
             }
             if (length(object@readsF) != length(object@readsR)){
                 "Same number of forward and reverse reads files needed to constitute files of paired end reads"}
             else{TRUE}
         })

##' @describeIn PairedReadFileSet-class Constructor for
##'     PairedReadFileSet-class
PairedReadFileSet <- function(readsF,
                              readsR, names=character()){
    ## construct from names if they exist
    if(length(names(readsF)) == length(readsF)) {
        na <- names(readsF)
    } ## overwrite if names are given as argument
    if(length(names) == length(readsF)) {
        na <- names
    } else { na <- basename(readsF)} ## otherwise use filenames
    new("PairedReadFileSet",
        readsF = readsF,
        readsR = readsR,
        names = na)
}

## Methods
##' @rdname PairedReadFileSet-class
setMethod("length", "PairedReadFileSet", function(x) length(x@readsF))

##' A class representing sequences of forward and reverse
##' primers.
##'
##' The PrimerPairsSet class is a container for storing primer pairs.
##' This means exactly two \code{\link[Biostrings]{DNAStringSet}}
##' objects of the same length specifying primer-pairs. Primer
##' sequences can be provided as character strings and will be
##' converted to \code{\link[Biostrings]{DNAStringSet}} by the
##' constructor function of the same name. primerF and primerR have to
##' be of the same length to specify primer pairs
##' (primerF[i].primerR[i]; for i in 1:length(PrimerF)).  Warnings are
##' given if primer sequences are of unusal length (<16 or >26 bases).
##'
##' @slot primerF DNAStringSet. Can be named or unnamed.
##' @slot primerR DNAStringSet of the same length. Can be named or
##'     unnnamed.
##' @slot names Either constructed as a concatenation of names of
##'     forward and reverse primers or of their sequences (if primer
##'     sequences are unnamed)
##' @slot .mapF (automatically generated) maps potentially duplicate
##'     entries in FW primers to unique entries.
##' @slot .mapR (auto-generated) maps potentially duplicate entries in
##'     FW primers to unique entries.
##' @slot .uniqueF (auto-generated) unique forward primers as
##'     character strings.
##' @slot .uniqueR (auto-generated) unique reverse primers as
##'     character strings.
##' 
##' @usage
##' ## Constructor:
##'
##' PrimerPairsSet(primerF, primerR)
##' 
##' @param primerF Character vector or DNAStringSet. Can be named or
##'     unnamed.
##' @param primerR Character vector or DNAStringSet of the same
##'     length. Can be named or unnamed.
##' 
##' @seealso \code{\link[Biostrings]{DNAStringSet}}
##' @importFrom Biostrings DNAStringSet
##' @return PrimerPairsSet-class
##' @author Emanuel Heitlinger
##' @export PrimerPairsSet
####  ##' @exportClass PrimerPairsSet-class
setClass("PrimerPairsSet", contains = "DNAStringSet",
         representation(primerF="DNAStringSet", primerR="DNAStringSet",
                        names="character", 
                        .mapF="numeric", .mapR="numeric",
                        .uniqueF="character", .uniqueR="character"),         
         validity=function(object) {
             if (any(width(c(object@primerF, object@primerR))<16)){
                 warning("short primer (<16nt) provided are you sure about your primer sequences?")
             }
             if (any(width(c(object@primerF, object@primerR))>26)){
                 warning("long primer (>26nt) provided are you sure about your primer sequences?")
             }
             if (length(object@primerF) != length(object@primerR)){
                 "Same number of forward and reverse primer sequences needed to constitute Primer-Pairs"}
})

##' @describeIn PrimerPairsSet-class Constructor for
##'     PrimerPairsSet-class
PrimerPairsSet <- function(primerF, primerR){
    ## if names exist construct primer names
    if(length(names(primerR)) == length(primerR) &&
       length(names(primerF))== length(primerF)) {
        na <- paste0(names(primerF), ".", names(primerR))
    } else {na <- paste0(primerF, ".", primerR)} # otherwise use primer sequences
    new("PrimerPairsSet",
        primerF = DNAStringSet(primerF),
        primerR = DNAStringSet(primerR),
        names = na,
        .mapF = as.numeric(factor(as.character(primerF))),
        .mapR = as.numeric(factor(as.character(primerR))),
        .uniqueF = sort(unique(primerF)),
        .uniqueR = sort(unique(primerR)))
}

## Methods
##' @rdname PrimerPairsSet-class
setMethod(length, "PrimerPairsSet", function(x) length(x@primerF))

##' @rdname PrimerPairsSet-class
setMethod(names, "PrimerPairsSet", function (x) x@names)

##' A pair of two derep objects
##'
##' derep-class objects as defined by the package \code{dada2}
##' (\code{\link[dada2]{derepFastq}}] are bundeled as forward and
##' reverse read pairs in this object
##' @title PiredDerep-class
##'
##' @slot derep object containing forward read pairs created by
##'     \code{dada2}'s \code{\link[dada2]{derepFastq}} function
##' @slot derepR derep object containing reverse read pairscreated by
##'     \code{dada2}'s \code{\link[dada2]{derepFastq}} function
##' @return A PairedDerep-class object
##' @author Emanuel Heitlinger
setClass("PairedDerep",
         slots = c(derepF="list", derepR="list"),
         validity=function(object) {
             if (length(object@derepF) != length(object@derepR)){
                 "Same number of forward and reverse derep objects needed to constitute forward and reverse sequence read pairs"
             }})

##' @rdname PairedDerep-class
setMethod("length", "PairedDerep", function(x){
    length(x@derepF)
})


## completely useless list does it all
## setClass("PairedDerepSet",
##          contains="list",
##          representation(),
##          prototype(elementType="PairedDerep"))


##' A pair of two dada objects 
##'
##' dada-class objects as defined by the package \code{dada2}
##' (function \code{\link[dada2]{dada}}) are bundeled as forward and
##' reverse read pairs in this object
##'
##' @title PairedDada-class
##' @author Emanuel Heitlinger
setClass("PairedDada",
         slots = c(dadaF="list", dadaR="list"),
         validity=function(object) {
             if (length(object@dadaF) != length(object@dadaR)){
                 "Same number of forward and reverse dada objects needed to constitute forward and reverse sequence read pairs"
             }})

##' @rdname PairedDada-class
PairedDada <- function(dadaF=list(), dadaR=list()){
    new("PairedDada",
        dadaF = dadaF,
        dadaR = dadaR)
}

##' @rdname PairedDada-class
setMethod("length", "PairedDada", function(x){
    length(x@dadaF)
})


## completely useless: list does it all
## setClass("PairedDadaSet",
##          contains="list",
##          representation(),
##          prototype(elementType="PairedDada"))



##' A class combining sequences of forward and reverse primers (in a
##' \code{\link{PrimerPairsSet-class}}) plus file names of paired end
##' sequencing files (in a \code{\link{PairedReadFileSet-class}}).
##'
##' The MultiAmplicon class is a container for storing primer pairs,
##' read files and processed data in an 'amplicon x samples'
##' format. Slots 
##' 
##' @slot PrimerPairsSet The primer pairs used in your experiment to
##'     specify amplicons stored in a
##'     \code{\link{PrimerPairsSet-class}} object.
##'
##' @slot PairedReadFileSet The (quality filtered) fastq files (one
##'     file pair for each sample) that store your sequencing data.
##'
##' @slot rawCounts
##'
##' @slot stratifiedFiles temporary files as a result of stratifying
##'     into amplicons and samples using the function
##'     \code{\link{sortAmplicons}}. Forward (sometimes called R1) and
##'     reverse (sometimes called R2) files are stored as a (amplicons
##'     x samples) matrix of \code{\link{PairedReadFileSet-class}}
##'     objects.
##'
##' @slot derep A \code{\link{PairedDerep-class}} object containing
##'     pairs of derep-class objects created by \code{dada2}’s
##'     \code{\link[dada2]{derepFastq}} function
##'
##' @slot dadaderep A \code{\link{PairedDada-class}} object containing
##'     pairs of dada-class objects created by \code{dada2}’s
##'     \code{\link[dada2]{dada}} function
##'
##' @slot mergers
##'
##' @slot sequenceTable
##'
##' @slot sequenceTableNoChime
##'
##' 
##' @usage
##' ## Constructors:
##'
##' MultiAmplicon(PrimerPairsSet, PairedReadFileSet)
##' 
##' @param PrimerPairsSet a set of primer pairs specifiying your
##'     amplicons see \code{\link{PrimerPairsSet-class}}
##' 
##' @param PairedReadFileSet a set of paired end sequencing data fiels
##'     \code{\link{PairedReadFileSet-class}}
##'
##' Accessor-like methods:
##' 
##' In the code snippets below, 'x' is an
##' \code{\link{MultiAmplicon-class}} object.
##' 
##' 'nrow(x)' An integer giving the number of primer pairs in the
##' object
##'
##' 'ncol(x)' An integer giving the number of paired read files
##' (usually samples and their replicates) in the object
##'
##' 'rownames(x)' A character vector giving the primer-pair names see
##' \code{\link{PrimerPairsSet}}
##'
##' 'colnames(x)' A character vector giving the names of the paired
##' read filnames (e.g. used to store sample names) see also
##' \code{\link{PairedReadFileSet}}
##'
##' 'rawCounts(x)' A matrix of raw sequencing read counts assigned to
##' a particular sample (input file pair) and amplicon (input primer
##' pair)
##' 
##' @seealso \code{\link[dada2]{derepFastq}},\code{\link[dada2]{dada}}
##' @importFrom dada2 derepFastq dada
##' @author Emanuel Heitlinger
##' @export MultiAmplicon
##' @exportClass MultiAmplicon
setClass("MultiAmplicon",
         representation(PrimerPairsSet="PrimerPairsSet",
                        PairedReadFileSet="PairedReadFileSet",
                        rawCounts="matrix",
                        stratifiedFiles="list",
                        derep="list",
                        dada="list",
                        mergers="list",
                        sequenceTable="list",
                        sequenceTableNoChime="list"))

##' @describeIn MultiAmplicon-class Constructor for
##'     MultiAmplicon-class
MultiAmplicon <- function(PrimerPairsSet = PrimerPairsSet(),
                          PairedReadFileSet = PairedReadFileSet(),
                          rawCounts = matrix(),
                          stratifiedFiles = list(),
                          derep = list(),
                          dada = list(),
                          mergers = list(),
                          sequenceTable = list(),
                          sequenceTableNoChime = list()
                          ){
    new("MultiAmplicon",
        PrimerPairsSet = PrimerPairsSet,
        PairedReadFileSet = PairedReadFileSet        
        )
}

##' @rdname MultiAmplicon-class
setMethod("colnames", "MultiAmplicon", function (x) x@PairedReadFileSet@names)

##' @rdname MultiAmplicon-class
setMethod("rownames", "MultiAmplicon", function (x) x@PrimerPairsSet@names)

##' @rdname MultiAmplicon-class
setMethod("ncol", "MultiAmplicon", function (x) length(x@PairedReadFileSet))


setMethod("nrow", "MultiAmplicon", function (x) length(x@PrimerPairsSet))

##' Access rawCounts of a MultiAmplicon object
##'
##' An accessor for the 'rawCounts' slot of a
##' \code{\link{MultiAmplicon-class}} object
##' 
##' @title rawCounts
##' @param x A \code{\link{MultiAmplicon-class}} object
##' @return A numerical matrix of counts for each amplicon
##'     (primer-pair) matched in each sample (paired read file)
##' @author Emanuel Heitlinger
##' @export
setGeneric("rawCounts", function(x){standardGeneric("rawCounts")})

##' @rdname rawCounts
setMethod("rawCounts", "MultiAmplicon", function(x) x@rawCounts)


