##' Populate a phyloseq object for 
##'
##' @title toPhyloseq
##' @param MA MultiAmplicon object with the \code{taxonTable} and
##'     \code{sequenceTableNoChime} slots filled.
##' @return a \code{phyloseq} object or a list of such objects
##' @author Emanuel Heitlinger
##' @export
setGeneric("toPhyloseq", function(MA, ...) {standardGeneric("toPhyloseq")})

##' Populate a phyloseq object with the contents of a MultiAmplicon
##' object
##'
##' Information stored in a \code{MultiAmplicon} object is transferred
##' to a \code{\link[phyloseq]{phyloseq}} object. 
##' 
##' @title toPhyloseq
##' @param MA MultiAmplicon object with the \code{taxonTable} and
##'     \code{sequenceTableNoChime} slots filled.
##' @param samples names of samples to be included in the phyloseq object
##' @return a \code{\link[phyloseq]{phyloseq}} object or a list of such objects
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
              all.tax <- as.matrix(Reduce(rbind, getTaxonTable(MA)))
              ## wrap it up into one Phyloseq object
              phyloseq(otu_table(allST, taxa_are_rows=FALSE),
                       tax_table(all.tax),
                       ...)
          })

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
    seqtab <- getSequenceTableNoChime(MA)
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
