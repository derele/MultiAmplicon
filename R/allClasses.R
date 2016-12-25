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
##' @rdname PairedReadFileSet-class
##' @title PairedReadFileSet-Class
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

PairedReadFileSet <- function(readsF,
                              readsR){
    if(length(names(readsF)) == length(readsF)) {
        na <- names(readsF)
    } else { na <- basename(readsF)}
    new("PairedReadFileSet",
        readsF = readsF,
        readsR = readsR,
        names = na)
}


## Methods
setMethod("length", "PairedReadFileSet", function(x) length(x@readsF))

##' A class representing sequences of forward and reverse
##' primers.
##'
##' The PrimerPairsSet class is a container for storing primer pairs.
##' This means exactly two \code{\link{DNAStrinSet}} objects of the
##' same length specifying primer-pairs. Primer sequences can be
##' provided as character strings and will be converted to
##' \code{\link{DNAStringSet}} by the constructor function of the same
##' name. primerF and primerR have to be of the same length to specify
##' primer pairs (primerF[i]:primerR[i]; for i in 1:length(PrimerF)).
##' Warnings are given if primer sequences are of unusal length (<16
##' or >26 bases).
##'
##' @slot primerF DNAStringSet. Can be named or unnamed.
##' @slot primerR DNAStringSet of the same length. Can be named or
##'     unnnamed.
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
##' ## Constructors:
##'
##' PrimerPairsSet(primerF, primerR)
##' 
##' @param primerF Character vector or DNAStringSet. Can be named or
##'     unnamed.
##' @param primerR Character vector or DNAStringSet of the same
##'     length. Can be named or unnamed.
##' 
##' @seealso \code{\link{DNAStringSet}}
##' @importFrom Biostrings DNAStringSet
##' @rdname PrimerPairsSet-class
##' @title PrimerPairsSet-class
##' @return PrimerPairsSet-class
##' @author Emanuel Heitlinger
##' @export PrimerPairsSet

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

## Constructor
PrimerPairsSet <- function(primerF = DNAStringSet(), primerR = DNAStringSet()){
    new("PrimerPairsSet",
        primerF = DNAStringSet(primerF),
        primerR = DNAStringSet(primerR),
        names = paste0(names(primerF), ":", names(primerR)),
        .mapF = as.numeric(factor(as.character(primerF))),
        .mapR = as.numeric(factor(as.character(primerR))),
        .uniqueF = sort(unique(primerF)),
        .uniqueR = sort(unique(primerR)))
}

## Methods
setMethod(length, "PrimerPairsSet", function(x) length(x@primerF))

setMethod(names, "PrimerPairsSet", function(x){
    if(length(x@primerF) == length(x@names)){
        x@names
    } else {paste0(x@primerF, ":", x@primerR)}
})


##' A class combining sequences of forward and reverse primers (in a
##' \code{\link{PrimerPairsSet-class}}) plus file names of paired end
##' sequencing files (in a \code{\link{PairedReadFileSet-class}}).
##'
##' @slot PrimerPairsSet The primer pairs used in your experiment to
##'     specify amplicons stored in a PrimerPairsSet-class object.
##'
##' @slot PairedReadFileSet The (quality filtered) fastq files (one
##'     file pair for each sample) that store your sequencing data.
##'
##' @slot rawCounts
##'
##' @slot FstratifiedFiles temporary forward (sometimes called R1)
##'     files as a result of stratifying into amplicons and samples
##'     using the function \code{\link{sortAmplicons}}
##'
##' @slot RstratifiedFiles temporary reverse (sometimes called R2)
##'     files as a result of stratifying into amplicons and samples
##'     using the function \code{\link{sortAmplicons}}
##'
##' @slot derep
##'
##' @slot dada
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
##' @description The PrimerPairsSet class is a container for storing
##'     primer pairs.
##' 
##' @seealso \code{\link{derep}},\code{\link{dada}}
##' @importFrom dada2 derepFastq dada
##' @rdname MultiAmplicon-class
##' @title MultiAmplicon-class
##' @return MultiAmplicon-class
##' @author Emanuel Heitlinger
##' @export MultiAmplicon
##' @exportClass MultiAmplicon

setClass("MultiAmplicon",
         representation(PrimerPairsSet="PrimerPairsSet",
                        PairedReadFileSet="PairedReadFileSet",
                        rawCounts="matrix",
                        stratifiedFiles="list",
                        derepF="list",
                        derepR="list",
                        dadaF="list",
                        dadaR="list",
                        mergers="list",
                        sequenceTable="list",
                        sequenceTableNoChime="list"))


MultiAmplicon <- function(PrimerPairsSet = PrimerPairsSet(),
                          PairedReadFileSet = PairedReadFileSet(),
                          rawCounts = matrix(),
                          stratifiedFiles = list(),
                          derepF = list(),
                          derepR = list(),
                          dadaF = list(),
                          dadaR = list(),
                          mergers = list(),
                          sequenceTable = list(),
                          sequenceTableNoChime = list()
                          ){
    new("MultiAmplicon",
        PrimerPairsSet = PrimerPairsSet,
        PairedReadFileSet = PairedReadFileSet        
        )
}

