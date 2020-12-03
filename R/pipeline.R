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
    stop(paste("this function is unfunctional at the moment use",
               "stratified files directly in dada/dadaMulti"))
    .complainWhenAbsent(MA, "stratifiedFiles")
    exp.args <- .extractEllipsis(list(...), nrow(MA))    
    PPderep <- mclapply(seq_along(MA@PrimerPairsSet), function (i){
        ## work on possilbe different paramters for this particular amplicon
        args.here <- lapply(exp.args, "[", i)
        .paramMessage("derepFastq", args.here)
        message("amplicon ", names(MA@PrimerPairsSet)[i],
                " dereplicating for ",
                length(MA@stratifiedFiles[[i]]@readsF), " of ",
                length(MA@PairedReadFileSet@readsF), " possible sample files.\n")
        derepF <- do.call(derepFastq,
                          c(list(getStratifiedFilesF(MA[i, ])), args.here))
        derepR <- do.call(derepFastq,
                          c(list(getStratifiedFilesR(MA[i,])), args.here))
        ## make it a list even if only one sample was dereplicated
        if (class(derepF)%in%"derep") {
            ## the same must be true for the revers then to keep them
            ## in the same length
            if(length(derepF$map) == 1 && keep.single.singlets == FALSE){
                derepF <- list()
                derepR <- list()
                message("producing empty derep object for amplicon ",
                    names(MA@PrimerPairsSet)[i],
                    "as only one sequence is reported for one sample ",
                    "set keep.single.singlets to TRUE to change this behaviour, ",
                    "but be warned that this may lead to downstream errors." )
            } else{
                derepF <- list(derepF)
                derepR <- list(derepR)
            }
        }
        Pderep <- lapply(seq_along(derepF), function (w){
            new("PairedDerep",
                derepF = derepF[[w]],
                derepR = derepR[[w]])
        })
        return(Pderep)
    }, mc.cores = mc.cores)
    names(PPderep) <- names(MA@PrimerPairsSet)
    initialize(MA, derep = PPderep)
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
    ## alternatives: 
    ### .complainWhenAbsent(MA, "derep")
    .complainWhenAbsent(MA, "stratifiedFiles")
    exp.args <- .extractEllipsis(list(...), nrow(MA))
    ## needs to be computed on pairs of amplicons
    PPdada <- mclapply(seq_along(MA@PrimerPairsSet), function (i){
        ## ## From a derep object 
        ## dF <- getDerepF(MA[i, ])
        ## dR <- getDerepR(MA[i, ])
        ## ## Or directly from stratified files
        dF <- getStratifiedFilesF(MA[i, ])
        dR <- getStratifiedFilesR(MA[i, ])
        message("\n\namplicon ", names(MA@PrimerPairsSet)[i],
           ": dada estimation of sequence variants from ",
            length(dF), " of ",
           length(MA@PairedReadFileSet), " possible sample files")
       if(length(dF)>0 && length(dR)>0){
           ## run functions for reverse and forward
           ## work on possilbe different paramters for this particular amplicon
           args.here <- lapply(exp.args, "[", i)
           .paramMessage("dada", args.here)
           dadaF <- do.call(dada, c(list(derep = dF, err=Ferr), args.here))
           ## make it a list of length 1 in case of only one sample,
           ## otherwise it is simplified and can't be handled
           if (class(dadaF)%in%"dada"){dadaF <- list(dadaF)}
           ## fix the names of the samples DANGEROUS! 
           names(dadaF) <- .fixSortedSampleNames(names(dadaF),
                                                 MA@PairedReadFileSet@names)
           dadaR <- do.call(dada, c(list(derep = dR, err=Rerr), args.here))
           ## make it a list in case of only one sample
           if (class(dadaR)%in%"dada"){dadaR <- list(dadaR)}
           ## fix the names of the samples DANGEROUS! 
           names(dadaR) <- .fixSortedSampleNames(names(dadaR),
                                                 MA@PairedReadFileSet@names)
           Pdada <- PairedDada(dadaF = dadaF, dadaR = dadaR)
       } else {
           Pdada <- PairedDada()
           message("\nskipping empty amplicon")
       }
       return(Pdada)
    }, mc.cores=mc.cores)
    names(PPdada) <- names(MA@PrimerPairsSet)
    initialize(MA, dada = PPdada)
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
    .complainWhenAbsent(MA, "dada")
    exp.args <- .extractEllipsis(list(...), nrow(MA))
    mergers <- mclapply(seq_along(MA@PrimerPairsSet), function (i){     
        daF <- getDadaF(MA[i, ])
        daR <- getDadaR(MA[i, ])
        ## ## From a derep object 
        ##  deF <- getDerepF(MA[i, ])
        ##  deR <- getDerepR(MA[i, ])
        ## ## Or directly from stratified files
        deF <- getStratifiedFilesF(MA[i, ])
        deR <- getStratifiedFilesR(MA[i, ])
        message("\nmerging sequences from " , length(MA@dada[[i]]),
                " samples for amplicon ",
                MA@PrimerPairsSet@names[[i]])
            ## work on possilbe different paramters for this particular amplicon
            args.here <- lapply(exp.args, "[", i)
            .paramMessage("mergePairs", args.here)
            MP <- do.call(mergePairs,
                          c(list(dadaF = daF, derepF = deF,
                                 dadaR = daR, derepR = deR), args.here))
            ## correct the case of one sample / amplicon 
            if(class(MP)%in%"data.frame"){MP <- list(MP)}
        return(MP)
    }, mc.cores=mc.cores)
    names(mergers) <- names(MA@PrimerPairsSet)
    initialize(MA, mergers = mergers)
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
    mergers <- getMergers(MA)
    ## hack to remove mergers of length 1 which don't get properly
    ## named dataframes    
    mergers[unlist(lapply(mergers, length ))==1] <- list(list())
    sequenceTable <- mclapply(seq_along(mergers), function (i){
        args.here <- lapply(exp.args, "[", i)
        .paramMessage("makeSequenceTable", args.here)
        do.call(makeSequenceTable, c(list(mergers[[i]]), args.here))
    }, mc.cores=mc.cores)
    names(sequenceTable) <- names(MA@PrimerPairsSet)
    initialize(MA, sequenceTable = sequenceTable)
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
    sequenceTableNoChime <-
    mclapply(seq_along(MA@sequenceTable), function (i) { 
            if (nrow(MA@sequenceTable[[i]])>0 && ncol(MA@sequenceTable[[i]])>0){
                args.here <- lapply(exp.args, "[", i)
                .paramMessage("removeBimeraDenovo", args.here)
                do.call(removeBimeraDenovo, c(list(MA@sequenceTable[[i]]), args.here))
            } else {matrix(nrow=0, ncol=0)}
        },
        mc.cores = mc.cores)
    names(sequenceTableNoChime) <- MA@PrimerPairsSet@names
    initialize(MA, sequenceTableNoChime = sequenceTableNoChime)
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
              sgt <- function(x) sum(getUniques(x))
              getN <- function(x) {
                  ## check length for every sample in amplicon for
                  ## mergers or just to have dada otherwise >1 as a
                  if(any(unlist(sapply(x, nrow)) > 0) ||
                     all(sapply(x, class)%in%"dada")){
                      cat("\n we execute\n")
                      cat("\n any(unlist(sapply(x, nrow)) > 0)", 
                          any(unlist(sapply(x, nrow)) > 0), "\n")
                      cat("\n ", "all(sapply(x, class)%in%dada)",
                          all(sapply(x, class)%in%"dada"), "\n")
                      sum(unlist(sapply(x, sgt, simplify=FALSE)))
                  } else {0; cat("\n we skipped\n")}
              }
              nMerged <- sapply(getMergers(MA, simplify=FALSE), getN)
              nBefore <- sapply(getDadaF(MA, simplify=FALSE), getN)
              prop <- nMerged/nBefore
              prop[is.nan(prop)] <- 0
              prop
          })



## Summary function for pipeline ------------------------------------

##' Obtain summary data for amplicons run through the MultiAmplicon pipeline.
##'
##' Get statistics on the number of samples (with read data), the
##' number of unique sequence variants and the number of reads left
##' after processing of amplicons in the MultiAmplicon pipeline. In
##' some steps of the pipeline dada2 performs quality filtering
##' excluding non-credible sequence variants.
##' 
##' @title getPipelineSummary
##' @param MA MultiAmplicon object with all slots filled for tracking.
##' @return a data.frame of sample, unique sequences and sequencing
##'     reads numbers per amplicon.
##' @importFrom plyr revalue
##' @export
##' @author Emanuel Heitlinger
getPipelineSummary <- function(MA){
    slots <- c("rawCounts", "stratifiedFiles", "dada", "mergers",
               "sequenceTable", "sequenceTableNoChime")
    slotFilled <- unlist(sapply(slots, function (x) length(slot(MA, x)))>0)
    ### helper functions
    getN <- function(x) unlist(lapply(x, function (y) sum(getUniques(y))))
    getU <- function(x) unlist(lapply(x, function (y) length(getUniques(y))))
    track.l <- lapply(seq_along(getDadaF(MA)), function (i) {
        samples <- list(
            sorted=length(getRawCounts(MA[i, ])[getRawCounts(MA[i, ])>0]),
            ##  ## derep is currently defunct 
            ##  derep=length(getDerepF(MA[i, ])),
            denoised = if(slotFilled["dada"]) {
                           length(getDadaF(MA[i,]))
                       } else{0},
            merged = if(slotFilled["mergers"]) {
                         length(getMergers(MA[i,]))
                     } else{0},
            tabulated = if(slotFilled["sequenceTable"]) {
                            nrow(getSequenceTable(MA[i,]))
                        } else{0},
            noChime= if(slotFilled["sequenceTableNoChime"]) {
                         nrow(getSequenceTableNoChime(MA[i,]))
                     } else{0}
        )
        uniques <- list(
            ##  ## derep is currently defunct 
            ##  derep=sum(getU(getDerepF(MA[i, ]))),
            denoised = if(slotFilled["dada"]) {
                           sum(getU(getDadaF(MA[i, ])))
                       } else{0},
            merged = if(slotFilled["mergers"]) {
                         sum(getU(getMergers(MA[i, ])))
                     } else{0},
            tabulated = if(slotFilled["sequenceTable"]) {
                            ncol(getSequenceTable(MA[i, ]))
                        } else{0},
            noChime = if(slotFilled["sequenceTableNoChime"]) {
                          ncol(getSequenceTableNoChime(MA[i, ]))
                      } else{0}
        )
        reads <- list(
            sorted=sum(getRawCounts(MA[i, ])),
            ##  ## derep is currently defunct 
            ##  derep=sum(getN(getDerepF(MA[i, ]))),
            denoised = if(slotFilled["dada"]) {
                           sum(getN(getDadaF(MA[i, ])))
                       } else{0},
            merged = if(slotFilled["mergers"]) {
                         sum(getN(getMergers(MA[i, ])))
                     } else {0}, 
            tabulated = if(slotFilled["sequenceTable"]) {
                            sum(getSequenceTable(MA[i, ]))
                        } else{0},
            noChime= if(slotFilled["sequenceTableNoChime"]) {
                         sum(getSequenceTableNoChime(MA[i, ]))
                     } else {0}
        )
        list(samples, uniques, reads)
    })
    track <- reshape::melt(track.l)
    track$L2 <- revalue(as.factor(track$L2), c("1"="samples",
                                               "2"="uniques", "3"="reads"))
    track$L3 <- factor(track$L3, levels = c("sorted", ## "derep",
                                            "denoised", "merged",
                                            "tabulated", "noChime"))
    names(track) <- c("value", "pipeStep", "what", "primer")
    return(track)
}

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

.fixSortedSampleNames <- function (oldNames, sampleNames) {
    ## most elegant would be one gigantic pattern, but sadly this
    ## fails with a strange (but known) C-level error. So have to
    ## break it down.
    ### pattern <- paste(sampleNames, collapse="|")
    pattern <- paste0(".*(", sampleNames, ").*")
    newNames <- oldNames
    for(pat in pattern) {
        newNames <- gsub(pat, "\\1", newNames)
    }
    newNames
}
