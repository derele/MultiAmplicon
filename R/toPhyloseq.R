##' Populate a phyloseq object with the contents of a MultiAmplicon
##' object
##'
##' Information stored in a \code{MultiAmplicon} object is transferred
##' to a \code{\link[phyloseq]{phyloseq}} object. 
##' 
##' @title toPhyloseq
##' @param MA MultiAmplicon object with the \code{taxonTable} and
##'     \code{sequenceTableNoChime} slots filled.
##' @param samples samples to include in phyloseq object
##' @param multi2Single should data from different amplicon be
##'     combined into one phyloseq object (TRUE) or should a list of
##'     seperate objects, one for each amplicon, be returned?
##' @param ... additional arguments to be passed to
##'     \code{\link[phyloseq]{phyloseq}}
##' @return a \code{\link[phyloseq]{phyloseq}} object
##' @importFrom phyloseq tax_table sample_data otu_table phyloseq
##' @author Emanuel Heitlinger
##' @export
toPhyloseq <- function(MA, samples, multi2Single=TRUE, ...){
    STL <- getSequenceTableNoChime(MA, simplify=FALSE)
    TTL <- getTaxonTable(MA, simplify=FALSE)
    if(length(TTL) == nrow(MA)){
        TAX <- TRUE
    } else if(length(TTL) == 0){
        message("No taxon table provided, so your phyloseq object will lack",
                " taxonomical annotations")
        TAX <- FALSE
    } else {
        stop("\nTaxon tables in provided in Multiamplicon object for are",
             " incongruent with the number of amplicons")
    }
    if(multi2Single){
        ## get sample tables filled with zeros for non-assessed
        ## samples for particular amplicons
        filledST <- lapply(STL, .fillSampleTables, samples=samples)
        allST <- as.matrix(Reduce(cbind, filledST))
        if (!all(dim(allST)>0)) {
            stop(paste("\nempty OTU table\n",
                       "rownames",  unlist(lapply(STL, base::rownames)),
                       "don't match sample names:",  samples, "\n"))
                  }
        ## The same for taxon annotations
        if(TAX){
            all.tax <- as.matrix(Reduce(rbind, TTL))
            ## to avoid problems with duplicated rownames (same sequences
            ## recovered for different amplicons), this can happen after trimming
            base::rownames(all.tax) <- make.unique(base::rownames(all.tax))
        }
        base::colnames(allST) <- make.unique(base::colnames(allST))
        ## wrap it up into one Phyloseq object
        phyloseq(otu_table(allST, taxa_are_rows=FALSE),
                 sample_data(MA@sampleData),
                 if (TAX) tax_table(all.tax),
                 ...)
    } else {
        PS.l <- lapply(seq_along(STL), function (i) {
            ## currently taxa tables are NULL if empty and
            ## sequence Tables have zero dimensions
            seqExists <- all(dim(STL[[i]])>0)
            if(TAX) {
                taxExists <- !is.null(TTL[[i]])
            } else {
                taxExists <- FALSE
            }
            if(seqExists) {                          
                allSampleTable <- .fillSampleTables(STL[[i]], samples=samples)
                phyloseq(otu_table(allSampleTable, taxa_are_rows=FALSE),
                         if(TAX && taxExists) tax_table(TTL[[i]]),
                         sample_data(MA@sampleData[base::rownames(allSampleTable),]),
                         ...)
            } else if(!isTRUE(seqExists) && !isTRUE(taxExists)){
                NULL
            } else  {
                stop(paste("inconsistent taxa data provided for",
                           "sequences in amplicon", i, 
                           rownames(MA)[[i]], "\n")
                     )
            }
        })
        ## somehow can't use rownames(MA)
        names(PS.l) <- rownames(MA)
        PS.l
    }
}
         


##' Add sample data to a MultiAmplicon object
##'
##' The sampleData slot is filled with a
##' \code{\link[phyloseq]{sample_data}} object from phyoseq created
##' merging sample names (colnames) of the MultiAmplicon object with a
##' data frame of sample data. The rownames of the that sampleData
##' data frame must correspond to colnames of the MultiAmplcion
##' object. 
##' 
##' @title addSampleData
##' @param MA A \code{\link{MultiAmplicon}} object
##' @param sampleData A data frame of providing data for samples in
##'     the \code{\link{MultiAmplicon}} object. This has to have the same rownames
##'     as the colnames of the MultiAmplicon object. If NULL, 
##'     (default) sampleData will be added based on names of the
##'     \code{link{PairedReadFileSet}} and the resulting sampleData
##'     will (only) give names of forward and reverse file names for
##'     each sample. 
##' @return A \code{\link{MultiAmplicon}} object with the sampleData
##'     slot filled.
##' @export
##' @author Emanuel Heitlinger
addSampleData <- function (MA, sampleData=NULL) {
    ## needed to provide compatibility with MA objects before
    ## sampleData slot was introduced
    if(is.null(sampleData)){
        MA@sampleData <- new("sample_data",
                             data.frame(row.names=rownames(MA),
                                        readsF=MA@PairedReadFileSet@readsF,
                                        readsR=MA@PairedReadFileSet@readsF))
        return(MA)
    } else {
        ## now merge... discard samples that are not either in the
        ## sequences or in the samples 
        missingSeq <- rownames(sampleData)[!rownames(sampleData)%in%
                                           rownames(MA@sampleData)]
        missingSamples <- rownames(MA@sampleData)[!rownames(MA@sampleData)%in%
                                                  rownames(sampleData)]
        if (length(missingSamples)>0) {
            warning(paste(length(missingSamples), "samples are missing from your",
                          "sampleData but seem to have sequence data reported.",
                          "They will be omitted if you continue"))
        }
        if (length(missingSeq)>0) {
            warning(paste(length(missingSeq), "samples are in your sampleData but",
                          "have no sequence data reported. They will be omitted,",
                          "if you continue"))
        }
        ## Forward and reverse read file-names are now always part of
        ## sampleData in a MA object, they are part of the merge when
        ## subsetting (see subset-methods)
        ## by.cols="row.names"
        ## if (all(c("readsF", "readsR")%in%colnames(sampleData))){
        ##     by.cols <- c(by.cols, "readsF", "readsR")
        ## }
        by.cols <- c("row.names", 
                     intersect(colnames(sampleData),
                               colnames(MA@sampleData)))
        mSampleData <- merge(MA@sampleData, sampleData, by=by.cols)
        rownames(mSampleData) <- mSampleData$Row.names
        mSampleData$Row.names <- NULL
        SData <- new("sample_data", mSampleData)
        MultiAmplicon(
            PrimerPairsSet = getPrimerPairsSet(MA),
            PairedReadFileSet = getPairedReadFileSet(MA),
            .Data=MA@.Data,
            sampleData = SData,
            stratifiedFilesF = getStratifiedFilesF(MA, dropEmpty=FALSE),
            stratifiedFilesR = getStratifiedFilesR(MA, dropEmpty=FALSE),
            derepF = getDerepF(MA, dropEmpty=FALSE),
            derepR = getDerepR(MA, dropEmpty=FALSE),
            dadaF = getDadaF(MA, dropEmpty=FALSE),
            dadaR = getDadaR(MA, dropEmpty=FALSE),
            mergers = getMergers(MA, dropEmpty=FALSE),
            sequenceTable = getSequenceTable(MA),
            sequenceTableNoChime = getSequenceTableNoChime(MA),
            taxonTable = getTaxonTable(MA)
    )
    }
}




##' Fill sequence tables in a MultiAmplicon object to include selected
##' samples for all amplicons.
##'
##' In a MultiAmplicon object for some primer pairs some samples might
##' have no amplified sequence variants at all. This function creates 
##' matrices including consistent samples for all the \code{sequenceTableNoChime}
##' slots of a MultiAmplicon object and returnes filled tables for all requested
##' samplese. Cells are filled with 0 (zeros) for samples originally
##' not recoverd in an amplicon.
##' 
##' @title .fillSampleTables
##' @param ST a sequence tabel from a MultiAmplicon object
##' @param samples a character vector giving names of samples to
##'     retain.
##' @return a sequence table for use in a MultiAmplicon object
##' @export
##' @author Emanuel Heitlinger
.fillSampleTables <- function (ST, samples){
    message("extracting ", length(samples), " requested samples")
    missing.samples <- samples[!samples%in%base::rownames(ST)]
    if(length(missing.samples)>0){
        fill <- matrix(0, nrow=length(missing.samples), ncol=ncol(ST))
        base::rownames(fill) <- missing.samples
        message("filling zeros for amplicon missing ", nrow(fill),
                " samples ")
        full <- rbind(ST, fill)
    } else {full <- ST}
    full[samples, ]
}

