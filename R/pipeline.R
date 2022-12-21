##' Dereplicate sequences in fastq files.
##'
##' An interface to \code{\link[dada2]{derepFastq}}, which itself uses
##' \code{\link[ShortRead]{FastqStreamer}} for dereplicating amplicon
##' sequences from fastq or compressed fastq
##' files. \code{\link[ShortRead]{FastqStreamer}} uses a streaming
##' interface and can be used on low memory machines (see \code{...}
##' below).
##'
##' @title derepMulti
##' @param MA MultiAmplicon object to be dereplicated.
##' @param keep.single.singlets logical argument indicating whether a
##'     derep slot should be filled with a single sequence recovered
##'     in only one sample. Defaults to FALSE meaning that empty derep
##'     objects are created in such a case. Keeping instead single
##'     sequence derep objects (setting this to \code{TRUE}) might
##'     result in downstream problems (Errors) in dada inference.
##' @param mc.cores number or compute cores for parallel dereplication
##'     of different amplicons (values > 1 are only allowed on
##'     Unix/Linux systems).
##' @param ... arguments to be passed to
##'     \code{\link[dada2]{derepFastq}}. All arguments to the function
##'     can be given as a vector of the same length as the number of
##'     primer pairs in the MultiAmplicon object, allowing to specify
##'     different parameters for each amplicon. If a shorter vector is
##'     given it will be recycled to match the number of amplicons. A
##'     general argument to consider changing is \code{n} (default
##'     1e+06, the maximum number of reads to parse at one
##'     time). Smaller values will decreas memory consumption.
##' @return MultiAmplicon object with paired derep slot (forward
##'     derepF and reverse derepR) filled.
##' @importFrom dada2 derepFastq
##' @importFrom parallel mclapply
##' @export
##' @author Emanuel Heitlinger
derepMulti <- function(MA, mc.cores = getOption("mc.cores", 1L),
                       keep.single.singlets = FALSE, ...){
    .complainWhenAbsent(MA, "stratifiedFilesF")
    exp.args <- .extractEllipsis(list(...), nrow(MA))    
    PPderep <- mclapply(seq_along(rownames(MA)), function (i){
        ## work on possilbe different paramters for this particular amplicon
        args.here <- lapply(exp.args, "[", i)
        .paramMessage("derepFastq", args.here)
        SF <- getStratifiedFilesF(MA[i, ])
        SR <- getStratifiedFilesR(MA[i, ])
        message("amplicon ", rownames(MA)[i],
                    " dereplicating for ",
                    length(SF), " of ",
                    length(getPairedReadFileSet(MA[i, ])), " possible sample files.\n")
        derepF <- do.call(derepFastq,
                              c(list(SF), args.here))
        if (class(derepF)%in%"derep") {
            derepF <- list(derepF)
            ## using the non-zero raw counts for naming, necessary as
            ## derep drops names when only sinlge files
            names(derepF) <- .deriveNames(MA[i, ], derepF)
        }
        if (!.isListOf(derepF, "derep", nullOk=FALSE)){
            stop("not a list of derep objects returned")
        }
        derepR <- do.call(derepFastq,
                          c(list(SR), args.here))
        if (class(derepR)%in%"derep") {
            derepR <- list(derepR)
            ## using the non-zero raw counts for naming, necessary as
            ## derep drops names when only sinlge files
            names(derepR) <- .deriveNames(MA[i, ], derepR) ## reinstating for BAD            
        }
        list(derepF, derepR)
    }, mc.cores = mc.cores)
    names(PPderep) <- rownames(MA)
    derepF_ampXsamples <- lapply(PPderep, "[[", 1)
    derepR_ampXsamples <- lapply(PPderep, "[[", 2)
    derepFmat <- .meltMASlotList(derepF_ampXsamples, MA)
    derepRmat <- .meltMASlotList(derepR_ampXsamples, MA)
    MultiAmplicon(
        .Data = MA@.Data,
        PrimerPairsSet = getPrimerPairsSet(MA),
        PairedReadFileSet = getPairedReadFileSet(MA),
        sampleData = getSampleData(MA),
        stratifiedFilesF = getStratifiedFilesF(MA, dropEmpty=FALSE),
        stratifiedFilesR = getStratifiedFilesR(MA, dropEmpty=FALSE),
        rawCounts = getRawCounts(MA),
        derepF = derepFmat,
        derepR = derepRmat
    )
}


##' A wrapper around \code{\link[dada2]{dada}} of the \code{dada2}
##' packager for multiple amplicons.
##'
##' The function runs \code{\link[dada2]{dada}} from the package
##'  \href{https://benjjneb.github.io/dada2/}{\code{dada2}} to perform
##'  'High resolution sample inference from amplicon data' on multiple
##'  amplicons stored as dereplicated sequences in a
##'  MultiAmplicon-class object.
##' 
##' @title dadaMulti
##' @param MA MultiAmplicon-class object.
##' @param mc.cores mc.cores number or compute cores for parallel
##'     denoising of different amplicons (values > 1 are only allowed
##'     on Unix/Linux systems). Also consider parallelization at the
##'     level of each individual amplicon via
##'     \code{\link[dada2]{dada}}.
##' @param Ferr As in the "err" parameter of dada: the matrix of
##'     estimated rates for each possible nucleotide change. In this
##'     case for forward reads.
##' @param Rerr the same for the reverse reads.
##' @param ... additional parameters to be passed to the
##'     \code{\link[dada2]{dada}} function of \code{dada2}. All
##'     arguments to dadaMulti can be given as a vector of the same
##'     length as the number of primer pairs in the MultiAmplicon
##'     object. This allows to specify different parameters for each
##'     amplicon. If a shorter vector is given it will be recycled to
##'     match the number of amplicons.
##' @return MultiAmplicon object with dadaF and dadaR slots filled.
##' @importFrom dada2 dada
##' @importFrom methods initialize new slot
##' @importFrom parallel mclapply
##' @export
##' @author Emanuel Heitlinger
dadaMulti <- function(MA, mc.cores=getOption("mc.cores", 1L),
                      Ferr=NULL, Rerr=NULL, ...){
    exp.args <- .extractEllipsis(list(...), nrow(MA))
    ## needs to be computed on pairs of amplicons
    PPdada <- mclapply(seq_along(rownames(MA)), function (i){
        ## ## From a derep object 
        dF <- getDerepF(MA[i, ], dropEmpty=TRUE)
        dR <- getDerepR(MA[i, ], dropEmpty=TRUE)
        if(length(dR)>0 && length(dR) >0) {
            message("\ndada estimation from dereplicated read objects")
        } else {
            ## Work directly from stratified files
            .complainWhenAbsent(MA, "stratifiedFilesF")
            message("\ndada estimation directly from stratified files")
            dF <- getStratifiedFilesF(MA[i, ], dropEmpty=TRUE)
            dR <- getStratifiedFilesR(MA[i, ], dropEmpty=TRUE)
        }
        message("\namplicon ", rownames(MA)[i],
           ": dada estimation of sequence variants from ",
            length(dF), " of ",
           length(getPairedReadFileSet(MA[i, ])), " possible samples")
       if(length(dF)>0 && length(dR)>0){
           ## run functions for reverse and forward
           ## work on possilbe different paramters for this particular amplicon
           args.here <- lapply(exp.args, "[", i)
           .paramMessage("dada", args.here)
           dadaF <- do.call(dada, c(list(derep = dF, err=Ferr), args.here))
           ## make it a list of length 1 in case of only one sample,
           ## otherwise it is simplified and can't be handled
           if (class(dadaF)%in%"dada"){
               dadaF <- list(dadaF)
               ## then there should be only one stratified file
               ## using the non-zero raw counts for naming
               names(dadaF) <- .deriveNames(MA[i, ], dadaF)
               if (length(dF) > 1) {
                   stop("dada collapsed but more than one stratified file found")
               }
           }
           if(!.isListOf(dadaF, "dada")) {
               stop(paste("incorrect classes reported for dadaF in amplicon",
                          rownames(MA)[i], "\n"))
           }
           ## using the non-zero raw counts for naming
           dadaR <- do.call(dada, c(list(derep = dR, err=Rerr), args.here))
           ## make it a list in case of only one sample
           if (class(dadaR)%in%"dada"){
               dadaR <- list(dadaR)
               ## using the non-zero raw counts for naming
               names(dadaR) <- .deriveNames(MA[i, ], dadaR)
               
           }
           if(!.isListOf(dadaR, "dada")) {
               stop(paste("incorrect classes reported for dadaR in amplicon",
                          rownames(MA)[i], "\n"))
           }           
       } else {
           dadaF <- list()
           dadaR <- list()
           message("\nskipping empty amplicon")
       }
       return(list(dadaF, dadaR))
    }, mc.cores=mc.cores)
    names(PPdada) <- rownames(MA)
    dadaF_ampXsamples <- lapply(PPdada, "[[", 1)
    dadaR_ampXsamples <- lapply(PPdada, "[[", 2)
    dadaFmat <- .meltMASlotList(dadaF_ampXsamples, MA)
    dadaRmat <- .meltMASlotList(dadaR_ampXsamples, MA)
    MultiAmplicon(
        PrimerPairsSet = getPrimerPairsSet(MA),
        PairedReadFileSet = getPairedReadFileSet(MA),
        sampleData = getSampleData(MA),
        stratifiedFilesF = getStratifiedFilesF(MA, dropEmpty=FALSE),
        stratifiedFilesR = getStratifiedFilesR(MA, dropEmpty=FALSE),
        rawCounts = getRawCounts(MA),
        derepF = getDerepF(MA, dropEmpty=FALSE),
        derepR = getDerepR(MA, dropEmpty=FALSE),
        dadaF = dadaFmat,
        dadaR = dadaRmat
    )
}

##' Merge denoised pairs of forward and reverse reads inside an
##' MultiAmplicon object.
##'
##' This is a wrapper for \code{\link[dada2]{mergePairs}} from
##' \code{dada2}. It works on an \code{\link{MultiAmplicon-class}}
##' object with derep and dada slots filled. Use
##' \code{\link{dadaMulti}} and \code{\link{derepMulti}} on a amplicon
##' sorted (see \code{\link{sortAmplicons}})
##' \code{\link{MultiAmplicon-class}} object to preprocess your
##' multi-marker data to this point.
##'
##' @title mergeMulti
##' @param MA \code{\link{MultiAmplicon-class}} object with derep and
##'     dada slots filled.
##' @param mc.cores number or compute cores for parallel merging of
##'     different amplicons (values > 1 are only allowed on Unix/Linux
##'     systems).
##' @param ... additional arguments to be passed to the
##'     \code{\link[dada2]{mergePairs}} function of \code{dada2}, all
##'     arguments to the function can be given as a vector of the same
##'     length as the number of primer pairs in the MultiAmplicon
##'     object, allowing to specify e.g. justConcatenate to be set
##'     TRUE for only some of the amplicons, or to specify different
##'     minOverlap for each amplicon. If a shorter vector is given it
##'     will be recycled to match the number of amplicons.
##' @return A MultiAmplicon-class object with the mergers slot filled.
##' @importFrom dada2 mergePairs
##' @importFrom parallel mclapply
##' @export
##' @author Emanuel Heitlinger
mergeMulti <- function(MA, mc.cores=getOption("mc.cores", 1L), ...){
    .complainWhenAbsent(MA, "dadaF")
    exp.args <- .extractEllipsis(list(...), nrow(MA))
    mergers <- mclapply(seq_along(MA@PrimerPairsSet), function (i){     
        ## The dada objects
        daF <- getDadaF(MA[i, ], dropEmpty=TRUE)
        daR <- getDadaR(MA[i, ], dropEmpty=TRUE)
        ## From a derep object if present
        deF <- getDerepF(MA[i, ], dropEmpty=TRUE)
        deR <- getDerepR(MA[i, ], dropEmpty=TRUE)
        if(length(deR)>0 && length(deR) >0) {
            message("\nmerging based on dereplicated read objects")
        } else {
            ## Work directly from stratified files
            .complainWhenAbsent(MA, "stratifiedFilesF")
            message("\nmerging directly from stratified files")
            deF <- getStratifiedFilesF(MA[i, ], dropEmpty=TRUE)
            deR <- getStratifiedFilesR(MA[i, ], dropEmpty=TRUE)
        }
        message("\nmerging sequences from " , length(daF),
                " samples for amplicon ",
                rownames(MA)[[i]])
        ## work on possilbe different paramters for this particular amplicon
        args.here <- lapply(exp.args, "[", i)
        .paramMessage("mergePairs", args.here)
        MP <- do.call(mergePairs,
                      c(list(dadaF = daF, derepF = deF,
                             dadaR = daR, derepR = deR), args.here))
        ## correct the case of one sample / amplicon 
        if(class(MP)%in%"data.frame"){
            MP <- list(MP)
            ## again get the (dropped) names from non-empty samples
            names(MP) <- .deriveNames(MA[i, ], MP)
        }
        return(MP)
    }, mc.cores=mc.cores)
    names(mergers) <- rownames(MA)
    mergersmat <- .meltMASlotList(mergers, MA)
    MultiAmplicon(
        PrimerPairsSet = getPrimerPairsSet(MA),
        PairedReadFileSet = getPairedReadFileSet(MA),
        .Data=MA@.Data,
        sampleData = getSampleData(MA),
        stratifiedFilesF = getStratifiedFilesF(MA, dropEmpty=FALSE),
        stratifiedFilesR = getStratifiedFilesR(MA, dropEmpty=FALSE),
        rawCounts = getRawCounts(MA),
        derepF = getDerepF(MA, dropEmpty=FALSE),
        derepR = getDerepR(MA, dropEmpty=FALSE),                        
        dadaF = getDadaF(MA, dropEmpty=FALSE),
        dadaR = getDadaR(MA, dropEmpty=FALSE),                        
        mergers = mergersmat
    )
}

##' Create a sequence table inside a MultiAmplicon object.
##'
##' This is a wrapper for the \code{\link[dada2]{makeSequenceTable}}
##' function of \code{dada2}. It works on an
##' \code{\link{MultiAmplicon-class}} object with the mergers slot
##' filled. Use \code{\link{dadaMulti}} and \code{\link{derepMulti}}
##' on an amplicon sorted (see \code{\link{sortAmplicons}})
##' \code{\link{MultiAmplicon-class}} object to preprocess your
##' multi-marker data to this point.
##'
##' @title makeSequenceTableMulti 
##'
##' @param MA \code{\link{MultiAmplicon-class}} object with mergers
##'     slot filled
##' @param mc.cores number of compute cores to use parallelizing
##'     computations over different amplicons. Only available on *nix
##'     systems.
##' @param ... additional parameters to be passed to the
##'     \code{\link[dada2]{makeSequenceTable}} function of
##'     \code{dada2}. All arguments to the function can be given as a
##'     vector of the same length as the number of primer pairs in the
##'     MultiAmplicon object, allowing to specify different parameters
##'     for each amplicon. If a shorter vector is given it will be
##'     recycled to match the number of amplicons.
##' @return A MultiAmplicon-class object with the sequenceTable slot
##'     filled
##' @importFrom dada2 makeSequenceTable
##' @importFrom parallel mclapply
##' @export
##' @author Emanuel Heitlinger
makeSequenceTableMulti <- function(MA, mc.cores=getOption("mc.cores", 1L), ...){
    .complainWhenAbsent(MA, "mergers")
    exp.args <- .extractEllipsis(list(...), nrow(MA))
    mergers <- apply(MA, 1, getMergers)
    ## ## hack to remove mergers of length 1 which don't get properly
    ## ## named dataframes    
    ## mergers[unlist(lapply(mergers, length ))==1] <- list(list())
    sequenceTable <- mclapply(seq_along(mergers), function (i){
        args.here <- lapply(exp.args, "[", i)
        .paramMessage("makeSequenceTable", args.here)
        do.call(makeSequenceTable, c(list(mergers[[i]]), args.here))
    }, mc.cores=mc.cores)
    names(sequenceTable) <- rownames(MA)
    MultiAmplicon(
        PrimerPairsSet = getPrimerPairsSet(MA),
        PairedReadFileSet = getPairedReadFileSet(MA),
        .Data=MA@.Data,
        sampleData = MA@sampleData,
        stratifiedFilesF = getStratifiedFilesF(MA, dropEmpty=FALSE),
        stratifiedFilesR = getStratifiedFilesR(MA, dropEmpty=FALSE),
        rawCounts = getRawCounts(MA),
        derepF = getDerepF(MA, dropEmpty=FALSE),
        derepR = getDerepR(MA, dropEmpty=FALSE),
        dadaF = getDadaF(MA, dropEmpty=FALSE),
        dadaR = getDadaR(MA, dropEmpty=FALSE),
        mergers = getMergers(MA, dropEmpty=FALSE),
        sequenceTable = sequenceTable
    )
}

##' Remove chimeric sequencing read pairs from a MultiAmplicon
##' object.
##'
##' This is a wrapper for the \code{\link[dada2]{removeBimeraDenovo}}
##' function of \code{dada2}. It works on an
##' \code{\link{MultiAmplicon-class}} object with the sequenceTable
##' slot filled. Use \code{\link{dadaMulti}},
##' \code{\link{derepMulti}}, \code{\link{mergeMulti}} and
##' \code{\link{makeSequenceTableMulti}} on a amplicon sorted (see
##' \code{\link{sortAmplicons}}) \code{\link{MultiAmplicon-class}}
##' object to preprocess your multi-marker data to this point.
##'
##' @title removeChimeraMulti
##' @param MA A \code{\link{MultiAmplicon-class}} object preprocessed
##'     to have the sequenceTable slot filled.
##' @param mc.cores integer number of cores to use for parallelization
##'     accross different amplicons.
##' @param ... passed on to
##'     \code{\link[dada2]{removeBimeraDenovo}}. All arguments to the
##'     function can be given as a vector of the same length as the
##'     number of primer pairs in the MultiAmplicon object, allowing
##'     to specify different parameters for each amplicon. If a
##'     shorter vector is given it will be recycled to match the
##'     number of amplicons.
##' @return a \code{\link{MultiAmplicon-class}} object with the
##'     \code{sequenceTableNoChime} filled
##' @importFrom dada2 removeBimeraDenovo
##' @importFrom parallel mclapply
##' @export
##' @author Emanuel Heitlinger
removeChimeraMulti <- function(MA, mc.cores = getOption("mc.cores", 1L), ...){
    .complainWhenAbsent(MA, "sequenceTable")
    exp.args <- .extractEllipsis(list(...), nrow(MA))
    ST <- getSequenceTable(MA)
    STN <- mclapply(seq_along(ST), function (i) { 
        if (nrow(ST[[i]])>0 && ncol(ST[[i]])>0){
            args.here <- lapply(exp.args, "[", i)
            .paramMessage("removeBimeraDenovo", args.here)
            do.call(removeBimeraDenovo, c(list(ST[[i]]), args.here))
        } else {matrix(nrow=0, ncol=0)}
    }, mc.cores = mc.cores)
    names(STN) <- rownames(MA)
    MultiAmplicon(.Data=MA@.Data,
                  PairedReadFileSet = getPairedReadFileSet(MA),
                  PrimerPairsSet = getPrimerPairsSet(MA),
                  sampleData = getSampleData(MA),
                  stratifiedFilesF = getStratifiedFilesF(MA, dropEmpty=FALSE),
                  stratifiedFilesR = getStratifiedFilesR(MA, dropEmpty=FALSE),
                  rawCounts = getRawCounts(MA),
                  derepF = getDerepF(MA, dropEmpty=FALSE),
                  derepR = getDerepR(MA, dropEmpty=FALSE),
                  dadaF = getDadaF(MA, dropEmpty=FALSE),
                  dadaR = getDadaR(MA, dropEmpty=FALSE),
                  mergers = getMergers(MA, dropEmpty=FALSE),
                  sequenceTable = getSequenceTable(MA),
                  sequenceTableNoChime = STN
                  )
}

##' Calculate the proportion of merged sequences for a MultiAmplicon
##' object.
##'
##' Merging using \code{\link{dadaMulti}} can result in loss of
##' sequences with too little overlap. This function determines what
##' proporton of sequences were merged using \code{\link{mergeMulti}}
##' for each amplicon in a MultiAmplicon object.
##' @title calcPropMerged
##' @param MA MultiAmplicon object with the \code{mergers} slot
##'     filled.
##' @return a vector of proportions of merged reads for each amplicon
##'     (potentially named).
##' @author Emanuel Heitlinger
##' @export
setGeneric("calcPropMerged", function(MA) {standardGeneric("calcPropMerged")})

##' @rdname calcPropMerged
##' @importFrom dada2 getUniques
##' @export
setMethod("calcPropMerged", "MultiAmplicon",
          function(MA){
              .complainWhenAbsent(MA, "mergers")

              input <- getCounts(getDadaF(MA,
                                          dropEmpty=FALSE, name=FALSE),
                                 "input") 
              merged <- getCounts(getMergers(MA,
                                             dropEmpty=FALSE, name=FALSE),
                                  "input")
              nBefore <- rowSums(input)
              nMerged <- rowSums(merged)
              prop <- nMerged/nBefore
              prop[is.nan(prop)] <- 0
              prop
          })


## Util functions for pipeline ------------------------------------
.extractEllipsis <- function(dotargs, n) {
    exp.args <- lapply(dotargs, function (x) {
        if(length(x)){
            if(n%%length(x)>0){
                stop("argument of length ", length(x),
                     " can't be recycled to number of amplicons (",
                     n, ")\n")
            } else{
                rep(x, times=n/length(x))
            }
        }
    })
    exp.args
}

.paramMessage <- function(what, args){
    ### would be better to print the do.call call directly... this
    ### hack to help with it as not implemented
    ## NULL as character
    args[unlist(lapply(args, is.null))] <- "NULL"
    ## other values as character
    args <- lapply(args, as.character)
    print.args.here <- paste(names(args), unlist(args),
                             sep = "=")
    if(!length(print.args.here)) {print.args.here <- "default"}
    message("calling ", what, " with ",
            paste(print.args.here, collapse = " "), " parameters")
}

#' @importFrom methods slotNames
.complainWhenAbsent <- function(MA, slotName){
    if(!slotName %in% slotNames(MA)) {
        stop("no valid slot name for class MultiAmplicon")
    }
    isSlot <- length(slot(MA, name=slotName))
    if(!isSlot){
        stop("Slot ", slotName, " not filled ",
             "for MultiAmplicon object: consult ?MultiAmplicon")
    }
}

## .fixSortedSampleNames <- function (oldNames, sampleNames) {
##     ## most elegant would be one gigantic pattern, but sadly this
##     ## fails with a strange (but known) C-level error. So have to
##     ## break it down.
##     ### pattern <- paste(sampleNames, collapse="|")
##     pattern <- paste0(".*(", sampleNames, ").*")
##     newNames <- oldNames
##     ##     message(paste("replacing", oldNames, "\n"))
##     for(pat in pattern) {
##         newNames <- gsub(pat, "\\1", newNames)
##         ##    message(paste("with", newNames, "\n"))
##     }
##     newNames
## }

.deriveNames <- function(MA, what){
        nnames <- colnames(MA)[getRawCounts(MA)>0]
        if(!length(nnames)==length(what)){
            stop(paste("incompatible names in amplicon",
                       rownames(MA), "for",  what ,"\n"))
        }
        nnames
}

.meltMASlotList <- function(MASlotList, MA){
    sapply(colnames(MA), function(samples) {
        sapply(rownames(MA), function(amp){
            MASlotList[[amp]][samples]
        })
    })
}
