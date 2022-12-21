##' Sort different amplicons into a fully stratified samples x
##' amplicons structure based on primer matches.
##'
##' This function uses \code{\link[Biostrings]{isMatchingStartingAt}}
##' to match primer sequences at the first position of forward and
##' reverse sequences. These primer sequences can be removed. The
##' remaining sequences of interest are written to files to allow
##' processing via standard metabarcoding pipelines.
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
##' @param rmPrimer logical, indicating whether primer sequences
##'     should be removed during sorting
##' @param filedir path to an existing or newly to be created folder
##'     on your computer. If existing it has to be empty.
##' @param ... additional parameter so be passed to
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
           def=function(MA, filedir="stratified_files",
                        n = 1e6, countOnly = FALSE, rmPrimer = TRUE, ...) {
               standardGeneric("sortAmplicons")
           })

##' @rdname sortAmplicons
setMethod(
    "sortAmplicons", "MultiAmplicon",
    function(MA, filedir="stratified_files",
             n = 1e6, countOnly = FALSE, rmPrimer = TRUE, ...){
        ## the data matrix of amplicons x samples stratified counts 
        PPS <- getPrimerPairsSet(MA)
        PRS <- getPairedReadFileSet(MA)

        NR <- length(PPS)
        NC <- length(PRS)
        data <- matrix(0, nrow=NR, ncol=NC)
        ## colnames are sample names taken from file names
        colnames(data) <- names(PRS)
        ## rownames have to come from (matched) primers
        rownames(data) <- names(PPS)
        
        namesGrid <- expand.grid(names(PPS), names(PRS))
        basenames <- apply(namesGrid, 1, paste, collapse="_")          
        
        filebaseF <- paste0(filedir, "/", basenames, "_F.fastq.gz")
        filebaseR <- paste0(filedir, "/", basenames, "_R.fastq.gz")

        stratFilesF <-matrix(filebaseF,
                             nrow=NR, ncol=NC, 
                             dimnames=list(names(PPS), names(PRS)))

        stratFilesR <-matrix(filebaseR,
                             nrow=NR, ncol=NC, 
                             dimnames=list(names(PPS), names(PRS)))
            
        readsF <- slot(PRS, "readsF")
        readsR <- slot(PRS, "readsR")

        ## test whether filedir exists
        if(!dir.exists(filedir)) {
            message("creating directory ", filedir)
            dir.create(filedir)
        } else {
            if (length(list.files("."))==0) {
                message("using existing directory ", filedir)
            } else {
                stop("directory for amplicon sorted files must be empty")
            }
        }
        ## add sample data and metadata in columns
        for(x in seq_along(readsF)) {
            if(!file.exists(readsF[[x]]) | !file.exists(readsR[[x]])){
                data[, colnames(data)[[x]]] <- 0
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
                              lengthF <- length(MA@PrimerPairsSet@primerF[[y]])
                              lengthR <- length(MA@PrimerPairsSet@primerR[[y]])
                              F <- Ffq[select]
                              R <- Rfq[select]
                              if (rmPrimer){
                                  F <- ShortRead::narrow(F,
                                                         lengthF,
                                                         width(Ffq[select]))
                                  R <- ShortRead::narrow(R,
                                                         lengthR,
                                                         width(Rfq[select]))
                              }
                              ShortRead::writeFastq(F,
                                                    file=stratFilesF[[y, x,
                                                                    drop=FALSE]],
                                                    mode="a")
                              ShortRead::writeFastq(R,
                                                    file=stratFilesR[[y, x,
                                                                      drop=FALSE]],
                                                    mode="a")
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
            msg <- paste("finished", doing, sum(data[, colnames(data)[[x]]]),
                         "sequencing reads for", colnames(data)[[x]], "in",
                         "\n", readsF[[x]], "and \n ", readsR[[x]], "\n",
                         "into", sum(data[, colnames(data)[[x]]]>0),
                         " amplicons")
            if(doing%in%"sorting"){
                msg <- paste(msg, "and written into",  normalizePath(filedir))
            }
            message(msg)
        }
        ## run only on existing files to avoid warnings for non-existing
        ## files. This means don't run on files corresponding to zeros
        ## read counts repored 
        if(!countOnly){
            MultiAmplicon(PrimerPairsSet = MA@PrimerPairsSet, 
                          PairedReadFileSet = MA@PairedReadFileSet,
                          .Data=MA@.Data,
                          sampleData=MA@sampleData,
                          stratifiedFilesF = stratFilesF,
                          stratifiedFilesR = stratFilesR,
                          rawCounts = data
                          )
        }else{
            data
        }
    })
