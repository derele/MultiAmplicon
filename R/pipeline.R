##' Dereplicate sequences in fastq files
##'
##' An interface to \code{\link[dada2]{derepFastq}} which itself uses
##' \code{\link[ShortRead]{FastqStreamer}} for dereplicating amplicon
##' sequences from fastq or compressed fastq files.
##'
##' @title derepMulti
##' @param MA MultiAmplicon object to be dereplicated
##' @param keep.single.singlets logical argument indicationg whether a
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
##'     match the number of ampicons.  
##' @return MultiAmplicon object with derep slots (forward derepF and
##'     reverse derepR) filled
##' @importFrom dada2 derepFastq
##' @importFrom parallel mclapply
##' @export
##' @author Emanuel Heitlinger
derepMulti <- function(MA, mc.cores = getOption("mc.cores", 2L),
                       keep.single.singlets = FALSE, ...){
    exp.args <- extract.ellipsis(list(...), nrow(MA))    
    PPderep <- mclapply(seq_along(MA@PrimerPairsSet), function (i){
        cat("amplicon", names(MA@PrimerPairsSet)[i],
            "dereplicating for ",
            length(MA@stratifiedFiles[[i]]@readsF), " of ",
            length(MA@PairedReadFileSet@readsF), "possible sample files\n")
        ## work on possilbe different paramters for this particular amplicon
        args.here <- lapply(exp.args, "[", i)
        param.message("derepFastq", args.here)
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
                cat("\nproducing empty derep object for amplicon",
                    names(MA@PrimerPairsSet)[i],
                    "as only one sequence is reported for one sample ",
                    "set keep.single.singlets to TRUE to change this behaviour,",
                    "but be warned that this may lead to downstream errors\n\n" )
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


##' A wrapper around \code{\link[dada2]{dada}} from \code{dada2}
##'
##' The function runs \code{\link[dada2]{dada}} from the package
##' \code{dada2} \url{https://benjjneb.github.io/dada2/} to performe
##' 'High resolution sample inference from amplicon data' on multiple
##' amplicons stored as dereplicated sequences in a
##' MultiAmplicon-class object
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
##'     match the number of ampicons.
##' @return MultiAmplicon object with dadaF and dadaR slots filled
##' @importFrom dada2 dada
##' @importFrom methods initialize new slot
##' @export
##' @author Emanuel Heitlinger
dadaMulti <- function(MA, Ferr=NULL, Rerr=NULL, ...){
    exp.args <- extract.ellipsis(list(...), nrow(MA))
    ## needs to be computed on pairs of amplicons
    PPdada <- lapply(seq_along(MA@PrimerPairsSet), function (i){
       dF <- lapply(MA@derep[[i]], function (x) slot(x, "derepF"))
       dR <- lapply(MA@derep[[i]], function (x) slot(x,  "derepR"))
       cat("\n\namplicon", names(MA@PrimerPairsSet)[i],
           "dada estimation of sequence variants from ",
            length(dF), " of ",
           length(MA@PairedReadFileSet), "possible sample files\n\n")
       if(length(dF)>0 && length(dR)>0){
           ## run functions for reverse and forward
           ## work on possilbe different paramters for this particular amplicon
           args.here <- lapply(exp.args, "[", i)
           param.message("dada", args.here)
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
           cat("skipping empty amplicon\n")
       }
       return(Pdada)
    })
    names(PPdada) <- names(MA@PrimerPairsSet)
    initialize(MA, dada = PPdada)
}

##' merge denoised pairs of forward and reverse reads inside an
##' MultiAmplcion object.
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
    exp.args <- extract.ellipsis(list(...), nrow(MA))
    mergers <- lapply(seq_along(MA@PrimerPairsSet), function (i){     
        if(length(MA@dada[[i]])){
            daF <- slot(MA@dada[[i]], "dadaF")
            daR <- slot(MA@dada[[i]], "dadaR")
            deF <- lapply(MA@derep[[i]], "slot", "derepF")
            deR <- lapply(MA@derep[[i]], "slot", "derepR")
            cat("merging sequences from " , length(MA@dada[[i]]),
                "samples for amplicon ",
                MA@PrimerPairsSet@names[[i]], "\n")

            ## work on possilbe different paramters for this particular amplicon
            args.here <- lapply(exp.args, "[", i)
            param.message("mergePairs", args.here)
            MP <- do.call(mergePairs,
                          c(list(dadaF = daF, derepF = deF,
                                 dadaR = daR, derepR = deR), args.here))
            ## correct the case of one sample / amplicon 
            if(class(MP)%in%"data.frame"){MP <- list(MP)}
            cat("DONE\n\n")
            return(MP)} else{
                          cat("skipping empty amplicon (sequences for" ,
                              length(MA@dada[[i]]), "samples)  ",
                              MA@PrimerPairsSet@names[[i]], "\n\n")
                          return(list())}
    })
    names(mergers) <- names(MA@PrimerPairsSet)
    initialize(MA, mergers = mergers)
}

##' Create a sequence table inside a MultiAmplcion object.
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
##'     recycled to match the number of ampicons.
##' @return A MultiAmplicon-class object with the sequenceTable slot
##'     filled
##' @importFrom dada2 makeSequenceTable
##' @export
##' @author Emanuel Heitlinger
sequenceTableMulti <- function(MA, ...){
    exp.args <- extract.ellipsis(list(...), nrow(MA))
    sequenceTable <- lapply(seq_along(MA@mergers), function (i){
        args.here <- lapply(exp.args, "[", i)
        param.message("makeSequenceTable", args.here)
        do.call(makeSequenceTable, c(list(MA@mergers[[i]]), args.here))
    })
    names(sequenceTable) <- MA@PrimerPairsSet@names
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
##' @title noChimeMulti
##' @param MA A \code{\link{MultiAmplicon-class}} object pre-processed
##'     to have a sequenceTableMulti populated
##' @param mc.cores integer number of cores to use for parallelization
##' @param ... paramter passed through to
##'     \code{\link[dada2]{removeBimeraDenovo}}. All arguments to the
##'     function can be given as a vector of the same length as the
##'     number of primer pairs in the MultiAmplicon object, allowing
##'     to specify different parameters for each amplicon. If a
##'     shorter vector is given it will be recycled to match the
##'     number of ampicons.
##' @return a \code{\link{MultiAmplicon-class}} object with the
##'     sequenceTableNoChime filled
##' @importFrom dada2 removeBimeraDenovo
##' @importFrom parallel mclapply
##' @export
##' @author Emanuel Heitlinger
noChimeMulti <- function(MA, mc.cores = getOption("mc.cores", 2L), ...){
    exp.args <- extract.ellipsis(list(...), nrow(MA))
    sequenceTableNoChime <-
    mclapply(seq_along(MA@sequenceTable), function (i) { 
            if (nrow(MA@sequenceTable[[i]])>0 && ncol(MA@sequenceTable[[i]])>0){
                args.here <- lapply(exp.args, "[", i)
                param.message("removeBimeraDenovo", args.here)
                do.call(removeBimeraDenovo, c(list(MA@sequenceTable[[i]]), args.here))
            } else {matrix()}
        },
        mc.cores = mc.cores)
    names(sequenceTableNoChime) <- MA@PrimerPairsSet@names
    initialize(MA, sequenceTableNoChime = sequenceTableNoChime)
}


extract.ellipsis <- function(dotargs, n) {
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

param.message <- function(what, args){
    ### would be better to print the do.call call directly... this
    ### hack to help with it as not implemented
    ## NULL as character
    args[unlist(lapply(args, is.null))] <- "NULL"
    ## other values as character
    args <- lapply(args, as.character)
    print.args.here <- paste(names(args), unlist(args),
                             sep = "=")
    if(!length(print.args.here)) {print.args.here <- "default"}
    cat("calling", what, "with", print.args.here, "parameters\n")
}

