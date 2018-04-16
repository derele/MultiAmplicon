##' Dereplicate sequences in fastq files.
##'
##' An interface to \code{\link[dada2]{derepFastq}} which itself uses
##' \code{\link[ShortRead]{FastqStreamer}} for dereplicating amplicon
##' sequences from fastq or compressed fastq files.
##'
##' @title derepMulti
##' @param MA MultiAmplicon object to be dereplicated
##' @param keep.single.singlets logical argument indicating whether a
##'     derep slot should be filled with a single sequence recovered
##'     in only one sample. Defaults to FALSE meaning that such empty
##'     derep objects are created in such a case. Keeping such derep
##'     objects (setting this to TRUE) might result in downstream
##'     problems (Errors) in dada inference.
##' @param mc.cores number or compute cores for parallel processing
##' @param ... arguments to be passed to \code{\link{derepFastq}}. All
##'     arguments to the function can be given as a vector of the same
##'     length as the number of primer pairs in the MultiAmplicon
##'     object, allowing to specify different parameters for each
##'     amplicon. If a shorter vector is given it will be recycled to
##'     match the number of amplicons.
##' @return MultiAmplicon object with derep slots (forward derepF and
##'     reverse derepR) filled
##' @importFrom dada2 derepFastq
##' @importFrom parallel mclapply
##' @export
##' @author Emanuel Heitlinger
derepMulti <- function(MA, mc.cores = getOption("mc.cores", 2L),
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
                          c(list(MA@stratifiedFiles[[i]]@readsF), args.here))
        derepR <- do.call(derepFastq,
                          c(list(MA@stratifiedFiles[[i]]@readsR), args.here))
        ## make it a list even if only one sample was dereplicated
        if (class(derepF)%in%"derep") {
            ## the same must be true for the revers then to keep them
            ## in the same lenghth
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


##' A wrapper around \code{\link[dada2]{dada}} from \code{dada2} for
##' multiple amplicons.
##'
##' The function runs \code{\link[dada2]{dada}} from the package
##' \code{dada2} \url{https://benjjneb.github.io/dada2/} to perform
##' 'High resolution sample inference from amplicon data' on multiple
##' amplicons stored as dereplicated sequences in a
##' MultiAmplicon-class object.
##' 
##' @title dadaMulti
##' @param MA MultiAmplicon-class object
##' @param Ferr As in the "err" parameter of dada: the matrix of
##'     estimated rates for each possible nucleotide change. In this
##'     case for forward reads.
##' @param Rerr the same for the reverse reads. 
##' @param ... additional parameters to be passed to
##'     \code{\link[dada2]{dada}} from \code{\link[dada2]{dada}}. All
##'     arguments to the function can be given as a vector of the same
##'     length as the number of primer pairs in the MultiAmplicon
##'     object, allowing to specify different parameters for each
##'     amplicon. If a shorter vector is given it will be recycled to
##'     match the number of amplicons.
##' @return MultiAmplicon object with dadaF and dadaR slots filled
##' @importFrom dada2 dada
##' @importFrom methods initialize new slot
##' @export
##' @author Emanuel Heitlinger
dadaMulti <- function(MA, Ferr=NULL, Rerr=NULL, ...){
    .complainWhenAbsent(MA, "derep")
    exp.args <- .extractEllipsis(list(...), nrow(MA))
    ## needs to be computed on pairs of amplicons
    PPdada <- lapply(seq_along(MA@PrimerPairsSet), function (i){
       dF <- lapply(MA@derep[[i]], function (x) slot(x, "derepF"))
       dR <- lapply(MA@derep[[i]], function (x) slot(x,  "derepR"))
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
               names(rawCounts(MA)[i, ])[rawCounts(MA)[i, ]>0]
           Pdada <- PairedDada(dadaF = dadaF, dadaR = dadaR)
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
##' @return A MultiAmplicon-class object with the mergers slot filled
##' @importFrom dada2 mergePairs
##' @export
##' @author Emanuel Heitlinger
mergeMulti <- function(MA, ...){
    .complainWhenAbsent(MA, "dada")
    exp.args <- .extractEllipsis(list(...), nrow(MA))
    mergers <- lapply(seq_along(MA@PrimerPairsSet), function (i){     
        if(length(MA@dada[[i]])){
            daF <- slot(MA@dada[[i]], "dadaF")
            daR <- slot(MA@dada[[i]], "dadaR")
            deF <- lapply(MA@derep[[i]], "slot", "derepF")
            deR <- lapply(MA@derep[[i]], "slot", "derepR")
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
##' @title sequenceTableMulti 
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
##' \code{\link{sequenceTableMulti}} on a amplicon sorted (see
##' \code{\link{sortAmplicons}}) \code{\link{MultiAmplicon-class}}
##' object to preprocess your multi-marker data to this point.
##'
##' @title removeChimeraMulti
##' @param MA A \code{\link{MultiAmplicon-class}} object preprocessed
##'     to have a sequenceTableMulti populated
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

##' Fill multiple sequence Tables in a MultiAmplicon object for all or
##' selected samples.
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
