##' A class representing file paths to paired end reads.
##'
##' Two character vectors of the same length specifying file names of
##' paired end reads can be stored in this class. Filenames a checked
##' for their existence.
##'
##' @slot readF The file path to forward reads
##' @slot readR The file path to reverse reads
##'
##' @usage
##' ## Constructors:
##' PairedReadFileSet(readsF, readsR)
##'
##' @param readsF The path and filenames of forward reads
##' @param readsR The path and filenames of reverse reads
##' 
##' @return PairedReadFileSet
##' @title PairedReadFileSet-Class
##' @author Emanuel Heitlinger
##'
##' @export PairedReadFileSet
setClass("PairedReadFileSet",
         slots = c(readsF="character", readsR="character"),
         contains = c(readsF="character", readsR="character"),
             validity=function(object) {
        all.files <- c(object@readsF, object@readsR)
        missing.file <- !file.exists(all.files)
        if (any(missing.file)){
            warning(paste0("\nfile ", all.files[missing.file], " do(es) not exist on your system"))
        }
        if (length(object@readsF) != length(object@readsR)){
            "Same number of forward and reverse reads files needed to constitute files of paired end reads"}
        else{TRUE}
             }
        )

PairedReadFileSet <- function(readsF = as.character(), readsR = as.character()){
    new("PairedReadFileSet",
        readsF = as.character(readsF),
        readsR = as.character(readsR))
}


## Methods
setMethod("length", "PairedReadFileSet", function(x) length(x@readsF))
setMethod("names", "PairedReadFileSet", function(x) basename(x@readsF))

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
##' @rdname PrimerPairsSet-class
##' @title PrimerPairsSet-class
##' @return PrimerPairsSet-class
##' @author Emanuel Heitlinger
##' @export PrimerPairsSet

setClass("PrimerPairsSet", contains = "DNAStringSet",
         representation(primerF="DNAStringSet", primerR="DNAStringSet",
                        textNames="character", 
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
        textNames = paste0(names(primerF), ":", names(primerR)),
        .mapF=as.numeric(factor(as.character(primerF))),
        .mapR=as.numeric(factor(as.character(primerR))),
        .uniqueF=sort(unique(primerF)),
        .uniqueR=sort(unique(primerR)))
}

## Methods
setMethod(length, "PrimerPairsSet", function(x) length(x@primerF))

setMethod(names, "PrimerPairsSet", function(x){
    paste0(x@primerF, ":", x@primerR)
})

textNames <- function(x){
    if (class(x)!= "PrimerPairsSet"){
        stop("please provide a PrimerPairsSet object")
    } else {
        x@textNames
    }
}

##' A class combining sequences of forward and reverse primers (in a
##' \code{\link{PrimerPairsSet-class}}) plus file names of paired end
##' sequencing files (in a \code{\link{PairedReadFileSet-class}}).
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
##' @seealso \code{\link{DNAStringSet}}
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
                        FstratifiedFiles="matrix",
                        RstratifiedFiles="matrix"))


MultiAmplicon <- function(PrimerPairsSet,
                          PairedReadFileSet, ...){
    new("MultiAmplicon",
        PrimerPairsSet = PrimerPairsSet,
        PairedReadFileSet = PairedReadFileSet        
        )
}

