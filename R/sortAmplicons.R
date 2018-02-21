##' Sort different amplicons into a fully stratified samples x
##' amplicons structure based on primer matches.
##'
##' This functions uses Biostrings::isMatchingStartingAt to match
##' primer sequences at the first postion of forward and reverse
##' sequences. These sequences are written to temporary files to allow
##' processing via standard metabarcoding pipelines (i.e. dada2 using
##' fastq streaming).
##' 
##' @title sortAmplicons
##' @param MA \code{\link{MultiAmplicon-class}} object containing a
##'     set of paired end files and a primer-pairs set.
##' @param n parameter passed to the yield functions of package
##'     ShortRead. This controls the memory consumption during
##'     streaming. Lower values result in lower memory requirements
##'     but might result longer processing time due to more repeated
##'     I/O operations reading the sequence files.
##' @param countOnly logical argument if set TRUE only a matrix of
##'     read counts is returned
##' @param filedir path to an existing or newly to be created folder
##'     on your computer. \code{\link[base]{tempfile}} is used within
##'     this folder to creat unique filnames trying to avoid problems
##'     in case the folder has been used before.
##' @param ... addtional parameter so be passed to
##'     Biostrings::isMatchingStartingAt. Be careful when using
##'     multiple starting positions or allowing error. This could lead
##'     to read pairs being assigned to multiple amplicons.
##' @return MultiAmplicon: By default (countOnly=FALSE) a
##'     \code{\link{MultiAmplicon-class}} object is returned with the
##'     stratifiedFiles slot populated. Stratified file names are
##'     constructed using a unique string created by
##'     \code{\link[base]{tempfile}} and stored in the given filedir
##'     (by default R's \code{\link[base]{tempdir}}). If the countOnly
##'     is set only a numeric matrix of read counts is returned.
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
##' ## sort into amplicons
##' MA1 <- sortAmplicons(MA)
##'
##' @rdname sortAmplicons
##' @author Emanuel Heitlinger
##' @importFrom ShortRead FastqStreamer yield narrow sread width
##' @importFrom Biostrings isMatchingStartingAt
##' @export sortAmplicons
##' @aliases sortAmplicons, sortAmplicons-Method
setGeneric(name="sortAmplicons",
           def=function(MA, n=1e6, countOnly=FALSE, filedir=tempdir(), ...) {
               standardGeneric("sortAmplicons")
           })

##' @rdname sortAmplicons
setMethod("sortAmplicons", "MultiAmplicon", function(MA, n=1e6, countOnly=FALSE,
                                                     filedir=tempdir(), ...){
    ## the data matrix of amplicons x samples stratified counts 
    NR <- length(MA@PrimerPairsSet@primerF)
    NC <- length(MA@PairedReadFileSet@readsF)
    data <- matrix(0, nrow=NR, ncol=NC)
    ## colnames are sample names taken from file names
    colnames(data) <- names(MA@PairedReadFileSet)
    ## rownames have to come from (matched) primers
    rownames(data) <- names(MA@PrimerPairsSet)
    tmppathF <- matrix("", nrow=NR, ncol=NC)
    tmppathR <- matrix("", nrow=NR, ncol=NC)
    readsF <- MA@PairedReadFileSet@readsF
    readsR <- MA@PairedReadFileSet@readsR
    ## test whether filedir exists
    if(!dir.exists(filedir)) {
        cat("creating directory ", filedir,
            "and adding specific prefixes to avoid problems when running code twice \n")
        dir.create(filedir)
    } else {
        cat("using existing directory ", filedir,
            "and adding specific prefixes to avoid problems when running code twice \n")
    }
    ## add sample data and metadata in columns
    for(x in seq_along(readsF)) {
        tmpbaseF <- paste(tempfile(tmpdir=filedir),
                          basename(readsF[[x]]), sep="_")
        tmpbaseR <- paste(tempfile(tmpdir=filedir),
                          basename(readsR[[x]]), sep="_")
        if(!file.exists(readsF[[x]]) | !file.exists(readsR[[x]])){
            data[, colnames(data)[[x]]] <- 0
            tmppathF[,x] <- paste0(tmpbaseF, names(MA@PrimerPairsSet),
                                   ".fastq.gz")
            tmppathR[,x] <- paste0(tmpbaseR, names(MA@PrimerPairsSet),
                                   ".fastq.gz")
            next()
        }
        f <- ShortRead::FastqStreamer(readsF[[x]], n = n)
        r <- ShortRead::FastqStreamer(readsR[[x]], n = n)
        ## request forward and reverse file simultaneously
         while(length(suppressWarnings(Ffq <- ShortRead::yield(f))) &&
               length(suppressWarnings(Rfq <- ShortRead::yield(r)))){
                   fM <- lapply(MA@PrimerPairsSet@.uniqueF, function(x){
                       as.vector(
                           Biostrings::isMatchingStartingAt(x,
                                                            ShortRead::sread(Ffq),
                                                            fixed=FALSE, ...))
                   })
                   rM <- lapply(MA@PrimerPairsSet@.uniqueR, function(x){
                       as.vector(
                           Biostrings::isMatchingStartingAt(x,
                                                            ShortRead::sread(Rfq),
                                                            fixed=FALSE, ...))
                   })
                   matches <- numeric(length=length(MA@PrimerPairsSet))
                   ## add primer pair data and metadata in rows 
                   for(y in seq_along(MA@PrimerPairsSet)) {
                       map.primerF <- MA@PrimerPairsSet@.mapF[[y]]
                       map.primerR <- MA@PrimerPairsSet@.mapR[[y]]
                       # is reverse and forward matched
                       select <- fM[[map.primerF]] & rM[[map.primerR]] 
                       if(!countOnly){ # file operations only if requested
                           tmppathF[y, x] <-
                               paste0(tmpbaseF,
                                      names(MA@PrimerPairsSet)[[y]],
                                      ".fastq.gz")
                           tmppathR[y, x] <-
                               paste0(tmpbaseR,
                                      names(MA@PrimerPairsSet)[[y]],
                                      ".fastq.gz")
                           lengthF <- length(MA@PrimerPairsSet@primerF[[y]])
                           lengthR <- length(MA@PrimerPairsSet@primerR[[y]])
                           F <- ShortRead::narrow(Ffq[select],
                                                  lengthF, width(Ffq[select]))
                           R <- Rfq[select]
                           R <- ShortRead::narrow(Rfq[select],
                                              lengthR, width(Rfq[select]))
                           ShortRead::writeFastq(F, file=tmppathF[[y, x]], mode="a")
                           ShortRead::writeFastq(R, file=tmppathR[[y, x]], mode="a")
                       }
                       matches[y] <- length(select[select==TRUE])
                   }
                   ## need to add over the while loop because of fastq streaming 
                   data[, colnames(data)[[x]]] <-
                       data[, colnames(data)[[x]]] + matches
               }
        close(f)
        close(r)
        doing <- ifelse(countOnly, "counting", "sorting")
        cat("\n finished ", doing, sum(data[, colnames(data)[[x]]]),
            "sequencing reads for", colnames(data)[[x]], "in",
            "\n ", readsF[[x]], " and \n ", readsR[[x]], "\n",
            " into ", sum(data[, colnames(data)[[x]]]>0), "amplicons \n" )
    }
    ## run only on existing files to avoid warnings for non-existing
    ## files. This means don't run on files corresponding to zeros
    ## read counts repored 
    if(!countOnly){
        stratifiedFiles <- lapply(seq_along(MA@PrimerPairsSet),
                                  function(i){
                                      existing <- which(data[i, ]>0)
                                      PairedReadFileSet(tmppathF[i, existing],
                                                        tmppathR[i, existing])
                                  })
        names(stratifiedFiles) <- names(MA@PrimerPairsSet)
        new.MA <- initialize(MA, rawCounts = data,
                             stratifiedFiles = stratifiedFiles)
        return(new.MA)
    }else{
        return(data)
    }
})


