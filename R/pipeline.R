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
##' @export
##' @author Emanuel Heitlinger
dadaMulti <- function(MA, Ferr=NULL, Rerr=NULL, ...){
    .complainWhenAbsent(MA, "derep")
    exp.args <- .extractEllipsis(list(...), nrow(MA))
    ## needs to be computed on pairs of amplicons
    PPdada <- lapply(seq_along(MA@PrimerPairsSet), function (i){
       dF <- getDerepF(MA[i, ])
       dR <- getDerepR(MA[i, ])
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
           dadaR <- do.call(dada, c(list(derep = dR, err=Rerr), args.here))
           ## make it a list in case of only one sample
           if (class(dadaR)%in%"dada"){dadaR <- list(dadaR)}
           ## naming the dada objects
           names(dadaF) <- names(dadaR) <-
               names(rawCounts(MA)[i, ])[rawCounts(MA)[i, ]>1]
           Pdada <- PairedDada(dadaF = dadaF, dadaR = dadaR)
##            names(Pdada) <- names(rawCounts(MA)[i, ])[rawCounts(MA)[i, ]>0]
       } else {
           Pdada <- PairedDada()
           message("\nskipping empty amplicon")
       }
       return(Pdada)
    })
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
##' @export
##' @author Emanuel Heitlinger
mergeMulti <- function(MA, ...){
    .complainWhenAbsent(MA, "dada")
    exp.args <- .extractEllipsis(list(...), nrow(MA))
    mergers <- lapply(seq_along(MA@PrimerPairsSet), function (i){     
        daF <- unlist(getDadaF(MA[i, ]), recursive=FALSE)
        daR <- unlist(getDadaR(MA[i, ]), recursive=FALSE)
        if(length(daF)>0 & length(daF)>0){
        deF <- unlist(getDerepF(MA[i, ]), recursive=FALSE)
        deR <- unlist(getDerepR(MA[i, ]), recursive=FALSE)        
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
            return(MP)} else{
                          message("skipping empty amplicon (sequences for" ,
                              length(MA@dada[[i]]), " samples)  ",
                              MA@PrimerPairsSet@names[[i]])
                          return(list())}
    })
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
##' @export
##' @author Emanuel Heitlinger
makeSequenceTableMulti <- function(MA, ...){
    .complainWhenAbsent(MA, "mergers")
    exp.args <- .extractEllipsis(list(...), nrow(MA))
    sequenceTable <- lapply(seq_along(MA@mergers), function (i){
        args.here <- lapply(exp.args, "[", i)
        .paramMessage("makeSequenceTable", args.here)
        do.call(makeSequenceTable, c(list(MA@mergers[[i]]), args.here))
    })
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
removeChimeraMulti <- function(MA, mc.cores = getOption("mc.cores", 2L), ...){
    .complainWhenAbsent(MA, "sequenceTable")
    exp.args <- .extractEllipsis(list(...), nrow(MA))
    sequenceTableNoChime <-
    mclapply(seq_along(MA@sequenceTable), function (i) { 
            if (nrow(MA@sequenceTable[[i]])>0 && ncol(MA@sequenceTable[[i]])>0){
                args.here <- lapply(exp.args, "[", i)
                .paramMessage("removeBimeraDenovo", args.here)
                do.call(removeBimeraDenovo, c(list(MA@sequenceTable[[i]]), args.here))
            } else {matrix()}
        },
        mc.cores = mc.cores)
    names(sequenceTableNoChime) <- MA@PrimerPairsSet@names
    initialize(MA, sequenceTableNoChime = sequenceTableNoChime)
}

##' Fill multiple sequence tables in a MultiAmplicon object to include
##' all or selected samples.
##'
##' In a MultiAmplicon object for some primer pairs some samples might
##' have no amplified sequence variants at all. This function adds a
##' \code{sequenceTableFilled} slot to a MultiAmplicon object. It uses
##' the \code{sequenceTableNoChime} slot and fills tables for all
##' samples appearing in any table or for requested samples.
##' 
##' @title fillSampleTables
##' @param MA MultiAmplicon object
##' @param samples either a character vector of length one giving
##'     "union" to use all samples appearing for any of the amplicons
##'     or a longer character vector giving names of samples to
##'     retain
##' @return MultiAmplicon object with the \code{sequenceTableFilled}
##'     slot filled
##' @export
##' @author Emanuel Heitlinger
fillSampleTables <- function (MA, samples="union"){
    .complainWhenAbsent(MA, "sequenceTableNoChime")
    seqtab <- getSequenceTableNoChime(MA)
    if (length(samples) == 1){
        if(samples %in% "union") {
            all.samples <- unique(unlist(lapply(seqtab, rownames)))
        } else {stop("please specify either \"union\" to use all samples", 
                     "or a list of sample names")
        }
    } else{all.samples <- samples}
    if(any(!all.samples%in%names(MA@PairedReadFileSet))){
        warning("requested samples " ,
                all.samples[!all.samples%in%names(MA@PairedReadFileSet)],
                " not found in original sample names")
    }    
    filledST <- lapply(seqtab, function (ampST){
        missing.samples <- all.samples[!all.samples%in%rownames(ampST)]
        if(length(missing.samples)>0){
            fill <- matrix(0, nrow=length(missing.samples), ncol=ncol(ampST))
            rownames(fill) <- missing.samples
            full <- rbind(ampST, fill)
        } else {full <- ampST}
        full[all.samples, ]
    })
    MA@sequenceTableFilled <- list() ## empty slot for subset operations
    initialize(MA[, which(names(MA@PairedReadFileSet)%in%all.samples)],
               sequenceTableFilled = filledST)
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
                  if(length(x)) {
                      sum(sapply(x, sgt))
                  } else {0}
              }
              nMerged <- sapply(getMergers(MA), getN)
              nBefore <- sapply(getDadaF(MA), getN)
              nMerged/nBefore
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
.complainWhenAbsent <- function(MA, slots=TRUE){
    allSlots <- slotNames(MA)
    slotL <- sapply(allSlots, function (x){
                    length(slot(MA, name=x))})
    fslot <- slotL[slots]
    if(!fslot){
        stop("Slot ", names(fslot), " not filled ",
             "for MultiAmplicon object: consult ?MultiAmplicon")
   }
}
