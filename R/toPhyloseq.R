##' Create a phyloseq object from a MultiAmplicon object
##'
##' @title toPhyloseq
##' @param MA MultiAmplicon object with the \code{taxonTable} and
##'     \code{sequenceTableNoChime} slots filled.
##' @param samples samples to include in phyloseq object
##' @param ... additional arguments to be passed to
##'     \code{\link[phyloseq]{phyloseq}}
##' @return a \code{\link[phyloseq]{phyloseq}} object
##' @author Emanuel Heitlinger
##' @export
setGeneric("toPhyloseq", function(MA, samples, ...) {standardGeneric("toPhyloseq")})

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
##' @param ... additional arguments to be passed to
##'     \code{\link[phyloseq]{phyloseq}}
##' @return a \code{\link[phyloseq]{phyloseq}} object
##' @author Emanuel Heitlinger
##' @export
setMethod("toPhyloseq", "MultiAmplicon",
          function(MA, samples, ...){
              .complainWhenAbsent(MA, "taxonTable")
              ## get sample tables filled with zeros for non-assessed
              ## samples for particular amplicons
              filledST <- .fillSampleTables(MA, samples=samples)
              allST <- as.matrix(Reduce(cbind, filledST))
              ## The same for taxon annotations
              all.tax <- as.matrix(Reduce(rbind, getTaxonTable(MA,
                                                               simplify=FALSE)))
              ## wrap it up into one Phyloseq object
              phyloseq(otu_table(allST, taxa_are_rows=FALSE),
                       tax_table(all.tax),
                       sample_data(MA@sampleData),
                       ...)
          })


##' Add sample data to a MultiAmplicon object
##'
##' The sampleData slot is filled with a
##' \code{\link[phylose]{sample_data}} object from phyoseq created
##' merging sample names (colnames) of the MultiAmplicon object with a
##' data frame of sample data. The rownames of the that sampleData
##' data frame must correspond to colnames of the MultiAmplcion
##' object. 
##' 
##' @title addSampleData
##' @param MA A \code{\link{MultiAmplicon}} object
##' @param sampleData A data frame of providing data for samples in
##'     the \code{\link{MultiAmplicon}} object. If set to \code{NULL}
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
                             data.frame(row.names=names(MA@PairedReadFileSet),
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
        initialize(MA, sampleData = SData)
    }
}




##' Fill sequence tables in a MultiAmplicon object to include selected
##' samples for all amplicons.
##'
##' In a MultiAmplicon object for some primer pairs some samples might
##' have no amplified sequence variants at all. This gets matching
##' matrices from the \code{sequenceTableNoChime} slot of a
##' MultiAmplicon object and returnes filled tables for all requested
##' samplese. Tables are filled with 0 (zeros) for samples originally
##' not recoverd in an amplicon.
##' 
##' @title .fillSampleTables
##' @param MA MultiAmplicon object
##' @param samples a character vector giving names of samples to
##'     retain.
##' @return MultiAmplicon object with the \code{sequenceTableFilled}
##'     slot filled
##' @noRd
##' @author Emanuel Heitlinger
.fillSampleTables <- function (MA, samples){
    message("extracting ", length(samples), " requested samples")
    .complainWhenAbsent(MA, "sequenceTableNoChime")
    seqtab <- getSequenceTableNoChime(MA, simplify=FALSE)
    filledST <- lapply(seqtab, function (ampST){
        missing.samples <- samples[!samples%in%rownames(ampST)]
        if(length(missing.samples)>0){
            fill <- matrix(0, nrow=length(missing.samples), ncol=ncol(ampST))
            rownames(fill) <- missing.samples
            full <- rbind(ampST, fill)
            message("filling zeros for amplicon missing ", nrow(fill),
                    " samples ")
        } else {full <- ampST}
        full[samples, ]
    })
    filledST
}

