##' A class representing file paths to paired end reads.
##'
##' Two character vectors of the same length specifying file names of
##' paired end reads can be stored in this class. Filenames a checked
##' for their existence. Usually these are the filnames of quality
##' filtered fastq files already stratified into samples (one file
##' pair for each sample). 
##'
##' @slot readF A character vector specifying the file paths to files
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
##' @slot names Is either (first choice) from names of the forward
##'     reads, or (if these are empty) constructed from basename of
##'     forward read files (filename without directory).
##' 
##' @return PairedReadFileSet
##' @author Emanuel Heitlinger
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


##' @param readsF The path and filenames containing forward (R1) reads
##' @param readsR The path and filenames containing reverse (R2) reads
##' @export PairedReadFileSet
##' @describeIn PairedReadFileSet-class Constructor for
##'     PairedReadFileSet-class
PairedReadFileSet <- function(readsF=character(), readsR=character()){
    ## construct from names if they exist
    if(length(names(readsF)) == length(readsF)) {
        na <- as.character(names(readsF)) ## as.c to catch empty (NULL)
    } else {na <- basename(readsF)} ## otherwise use filenames
    new("PairedReadFileSet",
        readsF = readsF,
        readsR = readsR,
        names = na)
}

## Methods
##' @rdname PairedReadFileSet-class
##' @export
setMethod("length", "PairedReadFileSet", function(x) length(x@readsF))

## Methods
##' @rdname PairedReadFileSet-class
##' @export
setMethod("names", "PairedReadFileSet", function(x) x@names)



##' A class representing sequences of forward and reverse primers.
##'
##' The PrimerPairsSet class is a container for storing primer pairs.
##' This means exactly two \code{\link[Biostrings]{DNAStringSet}}
##' objects of the same length specifying primer-pairs. Primer
##' sequences can be provided as character strings and will be
##' converted to \code{\link[Biostrings]{DNAStringSet}} by the
##' constructor function of the same name. primerF and primerR have to
##' be of the same length to specify primer pairs. Warnings are
##' given if primer sequences are of unusual length (<16 or >26 bases).
##'
##' @slot primerF DNAStringSet. Can be named or unnamed.
##' @slot primerR DNAStringSet of the same length. Can be named or
##'     unnamed.
##' @slot names Character string. Either constructed as a
##'     concatenation of names of forward and reverse primers or of
##'     their sequences (if primer sequences are unnamed).
##' @slot .mapF (automatically generated) maps potentially duplicate
##'     entries in FW primers to unique entries.
##' @slot .mapR (auto-generated) maps potentially duplicate entries in
##'     FW primers to unique entries.
##' @slot .uniqueF (auto-generated) unique forward primers as
##'     character strings.
##' @slot .uniqueR (auto-generated) unique reverse primers as
##'     character strings.
##' 
##' @seealso \code{\link[Biostrings]{DNAStringSet}}
##' @importFrom Biostrings DNAStringSet
##' @return PrimerPairsSet-class
##' @author Emanuel Heitlinger
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


##' @param primerF Character vector or DNAStringSet. Can be named or
##'     unnamed. 
##' @param primerR Character vector or DNAStringSet of the same
##'     length. Can be named or unnamed.
##' @export PrimerPairsSet
##' @rdname PrimerPairsSet-class
PrimerPairsSet <- function(primerF=character(), primerR=character()){
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
##' Accessor like functions
##' \code{length} gives the number of read pairs in a
##' \code{PrimerPairsSet-class} object
##' @param x A \code{PrimerPairsSet-class} object.
##' @export
##' @rdname PrimerPairsSet-class
setMethod(length, "PrimerPairsSet", function(x) length(x@primerF))

##' \code{names} of primer-pairs (amplicons) in a
##' \code{PrimerPairsSet-class} object
##' @export
##' @rdname PrimerPairsSet-class
setMethod(names, "PrimerPairsSet", function (x) x@names)

##' A pair of two derep objects
##'
##' derep-class objects as defined by the package \code{dada2}
##' (\code{\link[dada2]{derepFastq}}] are bundled as forward and
##' reverse read pairs in this object
##' @title PairedDerep-class
##'
##' @slot derepF derep object containing forward read pairs created by
##'     \code{dada2}'s \code{\link[dada2]{derepFastq}} function
##' @slot derepR derep object containing reverse read pairs created by
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

##' A pair of two dada objects 
##'
##' dada-class objects as defined by the package \code{dada2}
##' (function \code{\link[dada2]{dada}}) are bundled as forward and
##' reverse read pairs in this object
##'
##' @title PairedDada-class
##' @author Emanuel Heitlinger
setClass("PairedDada",
         slots = c(dadaF="list", dadaR="list"),
         validity=function(object) {
             dF <- slot(object, "dadaF")
             dR <- slot(object, "dadaR")
             if (!.isListOf(dF, "dada") || !.isListOf(dR, "dada")){
                 "only lists of dada class objects can form PairedDada-class objects"
             }
             if (length(dF) != length(dR)){"Same number of forward and reverse dada objects needed to constitute forward and reverse sequence read pairs"
             }
             if (any(names(dF) != names(dR))){
                 paste("all forward and reverse dada objects needed to be named exactly the same (being produced from the same samples, names are, Forward:", names(dF), "Reverse:",  names(dR), "\n")
             }
         })

##' @param dadaF a dada object of forward reads
##' @param dadaR a dada object of reverse reads
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


##' The central data structure of the MultiAmplicon package
##'
##' The MultiAmplicon class is a container that stores at least primer
##' pairs, read files and progressively processed data in an 'amplicon
##' x samples' format. The slots in this object are incrementally
##' filled with by running wrappers functions (mostly around functions
##' from the \code{dada2} package). The object is treated (subsetted
##' etc.) like a (pseudo) matrix, colums are samples, rows are
##' different amplicons.
##' 
##' @slot PrimerPairsSet The primer pairs used in your experiment to
##'     specify amplicons stored in a
##'     \code{\link{PrimerPairsSet-class}} object.
##'
##' @slot PairedReadFileSet The (quality filtered) fastq files (one
##'     file pair for each sample) that store your sequencing data.
##'
##' @slot rownames Names of the amplicons: derived from the names of
##'     the primer pairs. Created with the object, should not be
##'     supplied by the user.
##'
##' @slot colnames Names of the samples: these are names of the
##'     origial fastq files. Created with the object, should not be
##'     supplied by the user.
##'
##' @slot sampleData A sample_data object from
##'     \code{\link[phyloseq]{phyloseq}}. The slot is created from
##'     sample names (names of the \code{\link{PrimerPairsSet}} same
##'     as \code{colnames(MA)}) and more data can be added by
##'     \code{\link{addSampleData}}.
##' 
##' @slot rawCounts A numeric matrix of sequencing read counts per
##'     amplicon and sample. Created by the function
##'     \code{\link{sortAmplicons}} in the MultiAmplicon pipeline.
##'
##' @slot stratifiedFiles temporary files as a result of stratifying
##'     into amplicons and samples using the MultiAmplicon pipeline
##'     function \code{\link{sortAmplicons}}. Forward (sometimes
##'     called R1) and reverse (sometimes called R2) files are stored
##'     as a (amplicons x samples) matrix of
##'     \code{\link{PairedReadFileSet-class}} objects.
##'
##' @slot derep A list of \code{\link{PairedDerep-class}} objects
##'     containing pairs of derep-class objects created by
##'     \code{dada2}’s \code{\link[dada2]{derepFastq}} function or
##'     withing the MultiAmplicon pipeline by
##'     \code{\link{derepMulti}}.
##'
##' @slot dada A list of \code{\link{PairedDada-class}} object
##'     containing pairs of dada-class objects created by
##'     \code{dada2}’s \code{\link[dada2]{dada}} function. Within the
##'     MultiAmplicon pipeline this slot is filled by
##'     \code{\link{dadaMulti}}.
##'
##' @slot mergers A list of objects containing merged pairs of forward
##'     and reverse reads as created by by \code{dada2}’s
##'     \code{\link[dada2]{mergePairs}} function. Within the
##'     MultiAmplicon pipeline this slot is filled by
##'     \code{\link{mergeMulti}}.
##'
##' @slot sequenceTable A list of matrix objects created by
##'     \code{dada2}’s \code{\link[dada2]{makeSequenceTable}}.
##'     Samples (in rows) and amplified sequence variants (ASVs) in
##'     columns.  Within the MultiAmplicon pipeline this slot is
##'     filled by \code{\link{makeSequenceTableMulti}}.
##'
##' @slot sequenceTableNoChime A list of matrix objects created by
##'     \code{dada2}’s \code{\link[dada2]{removeBimeraDenovo}}.
##'     Samples (in rows) and ASVs screened for PCR chimeras in
##'     columns. Within the MultiAmplicon pipeline this slot is filled
##'     by \code{\link{removeChimeraMulti}}.
##'
##' @slot taxonTable A list of matrix objects created by a function
##'     for taxonomical annotation (for example
##'     \code{\link{blastTaxAnnot}}.  ASVs are in rows and taxnomical
##'     ranks are in columns. 
##'
##' MultiAmplicon(PrimerPairsSet, PairedReadFileSet)
##' 
##' @param PrimerPairsSet a set of primer pairs specifiying your
##'     amplicons see \code{\link{PrimerPairsSet-class}}
##' 
##' @param PairedReadFileSet a set of paired end sequencing data files
##'     \code{\link{PairedReadFileSet-class}}
##'
##' @param rawCounts Users should not supply this parameter, the slot
##'     is created by \code{\link{sortAmplicons}}.
##'
##' @param sampleData Users should not supply this parameter. It's
##'     filled with a sample_data object from
##'     \code{\link[phyloseq]{phyloseq}}. The slot is created from
##'     sample names (same as \code{colnames(MA)}) and more data can
##'     be added by \code{\link{addSampleData}}.
##'
##' @param stratifiedFiles Users should not supply this parameter, the
##'     slot is created by \code{\link{sortAmplicons}}.
##'
##' @param derep Users should not supply this parameter, the slot is
##'     created by \code{\link{derepMulti}}
##'
##' @param dada Users should not supply this parameter, the slot is
##'     created by \code{\link{dadaMulti}}
##' 
##' @param mergers Users should not supply this parameter, the slot is
##'     created by \code{\link{mergeMulti}}
##'
##' @param sequenceTable Users should not supply this parameter, the
##'     slot is created by \code{\link{makeSequenceTableMulti}}
##'
##' @param sequenceTableNoChime Users should not supply this parameter,
##'     the slot is created by \code{\link{removeChimeraMulti}}
##'
##' @param taxonTable Users should not supply this parameter, the slot
##'     is created by \code{\link{blastTaxAnnot}}. It's filled with a
##'     list of taxonomyTable objects from
##'     \code{\link[phyloseq]{phyloseq}}.
##' 
##' @examples
##'
##' primerF <- c("AGAGTTTGATCCTGGCTCAG", "ACTCCTACGGGAGGCAGC",
##'             "GAATTGACGGAAGGGCACC", "YGGTGRTGCATGGCCGYT")
##' primerR <- c("CTGCWGCCNCCCGTAGG", "GACTACHVGGGTATCTAATCC",
##'              "AAGGGCATCACAGACCTGTTAT", "TCCTTCTGCAGGTTCACCTAC")
##'
##' PPS <- PrimerPairsSet(primerF, primerR)
##' 
##' fastq.dir <- system.file("extdata", "fastq", package = "MultiAmplicon")
##' fastq.files <- list.files(fastq.dir, full.names=TRUE)
##' Ffastq.file <- fastq.files[grepl("F_filt", fastq.files)]
##' Rfastq.file <- fastq.files[grepl("R_filt", fastq.files)]
##'
##' PRF <- PairedReadFileSet(Ffastq.file, Rfastq.file)
##'
##' MA <- MultiAmplicon(PPS, PRF)
##'
##' ## the number of amplicons (primer pairs)
##' nrow(MA)
##'
##' ## the number of samples (sequencing read file pairs)
##' ncol(MA)
##'
##' ## sort into amplicons
##' MA1 <- sortAmplicons(MA, filedir=tempfile(pattern = "dir"))
##'
##' ## dereplication is currently not supported
##' ## MA2 <- derepMulti(MA1)
##'
##' ### use dada directly after sorting
##' MA3 <- dadaMulti(MA1, selfConsist = TRUE)
##'
##' MA4 <- mergeMulti(MA3, justConcatenate=TRUE)
##'
##' MA5 <- makeSequenceTableMulti(MA4)
##'
##' MA6 <- removeChimeraMulti(MA5, mc.cores=1)
##'
##' @seealso \code{\link[dada2]{derepFastq}},\code{\link[dada2]{dada}}
##' @importFrom dada2 derepFastq dada
##' @importClassesFrom phyloseq sample_data
##' @importFrom phyloseq sample_names
##' @author Emanuel Heitlinger
##' @exportClass MultiAmplicon
setClass("MultiAmplicon",
         representation(PrimerPairsSet="PrimerPairsSet",
                        PairedReadFileSet="PairedReadFileSet",
                        rownames = "character",
                        colnames = "character",
                        rawCounts="matrix",
                        stratifiedFiles="list",
                        sampleData="sample_data", 
                        derep="list",
                        dada="list",
                        mergers="list",
                        sequenceTable="list",
                        sequenceTableNoChime="list",
                        taxonTable="list"),
         validity=function(object) {
             ## constructors check for PairedReadFileSet,
             ## PrimerPairsSet, rawCounts and sampleData
             ## directly. Other slots are list can additionally be
             ## checked at deeper levels
             if(!all(rownames(object)%in%sample_names(object@sampleData))){
                 "Sample_names of SampleData must be equeal rownames of  MultiAmplicon object"
             }
             if(!.isListOf(object@dada, "PairedDada")){
                 "PairedDada objects or an empty list must be provided for a valid MultiAmplicon object"
             }
             if(!all(unlist(lapply(object@derep, .isListOf, "PairedDerep")))){
                 "PairedDerep objects or an empty list must be provided for a valid MultiAmplicon object"
             }
             if(!all(unlist(lapply(object@mergers, .isListOf, "data.frame")))){
                 "data.frame objects or an empty list must be provided for valid mergers in MultiAmplicon object"
             }
             if(!.isListOf(object@sequenceTable, "matrix")){
                 "matrix objects or an empty list must be provided for valid sequenceTables in MultiAmplicon object"
             }
             if(!.isListOf(object@sequenceTableNoChime, "matrix")){
                     "matrix objects or an empty list must be provided for valid sequenceTableNoChime in MultiAmplicon object"
             }
             if(!.isListOf(object@taxonTable, "taxonomyTable", nullOk=TRUE)){
                 "taxonomyTable objects or an empty list must be provided in MultiAmplicon object"
             }
         }
         )


##' @export MultiAmplicon
##' @describeIn MultiAmplicon-class Constructor for
##'     MultiAmplicon-class
MultiAmplicon <- function(PrimerPairsSet = PrimerPairsSet(),
                          PairedReadFileSet = PairedReadFileSet(),
                          rawCounts = matrix(ncol=0, nrow=0),
                          stratifiedFiles = list(),
                          sampleData = new("sample_data",
                                           data.frame(row.names=names(PairedReadFileSet),
                                                      readsF=PairedReadFileSet@readsF,
                                                      readsR=PairedReadFileSet@readsF)),
                          derep = list(),
                          dada = list(),
                          mergers = list(),
                          sequenceTable = list(),
                          sequenceTableNoChime = list(),
                          taxonTable = list()
                          ){
    new("MultiAmplicon",
        PrimerPairsSet = PrimerPairsSet,
        PairedReadFileSet = PairedReadFileSet,
        rownames = PrimerPairsSet@names,
        colnames = PairedReadFileSet@names,        
        rawCounts = rawCounts,
        stratifiedFiles = stratifiedFiles,
        sampleData = sampleData,
        derep = derep,
        dada = dada,
        mergers = mergers,
        sequenceTable = sequenceTable,
        sequenceTableNoChime = sequenceTableNoChime,
        taxonTable = taxonTable
        )
}

##' @rdname MultiAmplicon-class
##' @export
setMethod("colnames", "MultiAmplicon", function (x) x@colnames)

##' @rdname MultiAmplicon-class
##' @export
setMethod("rownames", "MultiAmplicon", function (x) x@rownames)

##' @rdname MultiAmplicon-class
##' @export
setMethod("ncol", "MultiAmplicon", function (x) length(x@PairedReadFileSet))

##' @rdname MultiAmplicon-class
##' @export
setMethod("nrow", "MultiAmplicon", function (x) length(x@PrimerPairsSet))

##' @rdname MultiAmplicon-class
##' @export
setMethod("dim", "MultiAmplicon", function (x) {
    as.integer(c(length(x@PrimerPairsSet), length(x@PairedReadFileSet)))
})


##' @rdname MultiAmplicon-class
##' @param MA MultiAmplicon-class object
##' @export
getPrimerPairsSet <- function (MA) slot(MA, "PrimerPairsSet")


##' @rdname MultiAmplicon-class
##' @param MA MultiAmplicon-class object
##' @export
getPairedReadFileSet <- function (MA) slot(MA, "PairedReadFileSet")


##' @rdname MultiAmplicon-class
##' @param MA MultiAmplicon-class object
##' @export
getRawCounts <- function (MA) slot(MA, "rawCounts")

##' @rdname MultiAmplicon-class
##' @param simplify Should a list of objects be simplified to only one
##'     object if it has length one?
##' @export
getStratifiedFilesF <- function(MA, simplify=TRUE) {
    .simpfy(lapply(MA@stratifiedFiles, slot, "readsF"), simplify)
 }

##' @rdname MultiAmplicon-class
##' @export
getStratifiedFilesR <- function(MA, simplify=TRUE) {
    .simpfy(lapply(MA@stratifiedFiles, slot, "readsR"), simplify)
}
    
##' @rdname MultiAmplicon-class
##' @export
getDerepF <-  function(MA, simplify=TRUE) {
    .simpfy(lapply(MA@derep, function (x) lapply(x, slot, "derepF")),
            simplify)
}

##' @rdname MultiAmplicon-class
##' @export
getDerepR <-  function(MA, simplify=TRUE) {
    .simpfy(lapply(MA@derep, function (x) lapply(x, slot, "derepR")),
            simplify)
}

##' @rdname MultiAmplicon-class
##' @export
getDadaF <- function(MA, simplify=TRUE) {
    .simpfy(lapply(MA@dada, slot, "dadaF"), simplify)
}

##' @rdname MultiAmplicon-class
##' @export
getDadaR <- function(MA, simplify=TRUE) {
    .simpfy(lapply(MA@dada, slot, "dadaR"), simplify)
}

##' @rdname MultiAmplicon-class
##' @export
getMergers <- function(MA, simplify=TRUE) .simpfy(MA@mergers, simplify)

##' @rdname MultiAmplicon-class
##' @export
getSequenceTable <- function(MA, simplify=TRUE){
    .simpfy(MA@sequenceTable, simplify)
}

##' @rdname MultiAmplicon-class
##' @export
getSequenceTableNoChime <- function(MA, simplify=TRUE){
    .simpfy(MA@sequenceTableNoChime, simplify)
}

##' @rdname MultiAmplicon-class
##' @export
getTaxonTable <- function(MA, simplify=TRUE){
    .simpfy(MA@taxonTable, simplify)
}




##' @rdname MultiAmplicon-class
##' @export
setGeneric("getSequencesFromTable", function(MA) {standardGeneric("getSequencesFromTable")})

##' @rdname MultiAmplicon-class
##' @importFrom dada2 getSequences
##' @export
setMethod("getSequencesFromTable", "MultiAmplicon",
          function (MA) {
              lapply(getSequenceTableNoChime(MA), function (y) {
                  if(all(dim(y)>1)) getSequences(y) else NULL
              })
          })

.simpfy <- function(x, simplify) {
    if(length(x)==1 && simplify){
        x[[1]]
    } else{x}
}
    

.isListOf <- function (x, what, nullOk=FALSE){
    if(nullOk) {
        what <- c(what, "NULL")
    }
    if(!is.list(x)){
        return(FALSE)
    }
    else{
        classes <- unlist(lapply(x, class))
        all(classes%in%what)
    }
}





## Prevent R CMD check from complaining about the use of pipe expressions
## standard data.table variables
if (getRversion() >= "2.15.1")
    utils::globalVariables(c(".", ".I", ".N", ".SD",
                             "query", "subject", "pident",
                             "length", "mismatch",
                             "gapopen", "qstart", "qend", "sstart", "send",
                             "evalue",  "bitscore", "staxid",
                             "qtaxid", "bitsum", "ampProd",
                             "superkingdom", "phylum", "class", "order", "family",
                             "genus", "species"), utils::packageName(), add=FALSE)

## Make sure data.table knows we know we're using it
.datatable.aware = TRUE

