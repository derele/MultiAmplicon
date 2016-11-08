##' A class representing just file path to paired end reads, meaning
##' basically two character vectors of the same length. 
##'
##' Filenames a checked for existence.  
##' @title PairedReadFileSet-Class
##'
##' @usage
##' ## Constructors:
##' PairedReadFileSet(readsF, readsR)
##'
##' @param readsF The path and filenames of forward reads
##' @param readsR The path and filenames of reverse reads
##' @return Constructors for defined classes 
##' @author Emanuel Heitlinger


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
##' rimers. Exactly two \code{\link{DNAStrinSet}} objects of the same
##' length specifying primer-pairs. Primer sequences can be provided
##' as character strings and will be converted to
##' \code{\link{DNAStringSet}} by the function of the same
##' name. primerF and primerR have to be of the same length to specify
##' primer pairs (primerF[i]:primerR[i]; for i in 1:length(PrimerF)).
##' Warnings are given if primer sequences are of unusal length (<16
##' or >26 bases).
##'
##' @slot primerF DNAStringSet. Can be named or unnamed.
##' @slot primerR DNAStringSet of the same length. Can be named or unnnamed.
##' @slot .mapF Mapps potentially duplicate entries in FW primers to unique entries.
##' @slot .mapR Mapps potentially duplicate entries in FW primers to unique entries.
##' @slot .uniqueF unique forward primers as character strings.
##' @slot .uniqueR unique reverse primers as character strings.
##' 
##' @usage
##' ## Constructors:
##'
##' PrimerPairsSet(primerF, primerR)
##' 
##' @param primerF Character vector or DNAStringSet. Can be named or unnamed.
##' @param primerR Character vector or DNAStringSet of the same length. Can be named or unnamed.
##' 
##' @description The PrimerPairsSet class is a container for storing primer pairs. 
##' 
##' @seealso \code{\link{DNAStringSet}}
##' @rdname PrimerPairsSet-class
##' @title PrimerPairsSet-class
##' @return PrimerPairsSet-class
##' @author Emanuel Heitlinger

setClass("PrimerPairsSet", contains = "DNAStringSet",
         representation(primerF="DNAStringSet", primerR="DNAStringSet",
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


##' A class representing sequences of forward and reverse
##' rimers. Exactly two \code{\link{DNAStrinSet}} objects of the same
##' length specifying primer-pairs. Primer sequences can be provided
##' as character strings and will be converted to
##' \code{\link{DNAStringSet}} by the function of the same
##' name. primerF and primerR have to be of the same length to specify
##' primer pairs (primerF[i]:primerR[i]; for i in 1:length(PrimerF)).
##' Warnings are given if primer sequences are of unusal length (<16
##' or >26 bases).
##' 
##' @usage
##' ## Constructors:
##'
##' PrimerPairsSet(primerF, primerR)
##' 
##' @param primerF Character vector (then converted to) or DNAStringSet. Can be named or unnamed.
##' @param primerR Character vector (then converted to) or DNAStringSet of the same length. Can be named or unnamed.
##' 
##' @description The PrimerPairsSet class is a container for storing primer pairs. 
##' 
##' @seealso \code{\link{DNAStringSet}}
##' @rdname PrimerPairsSet-class
##' @title PrimerPairsSet-class
##' @return PrimerPairsSet-class
##' @author Emanuel Heitlinger
##' @importClassesFrom Biostrings DNAStringSet

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

