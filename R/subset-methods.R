setMethod("[", "PrimerPairsSet", function(x, i){
    newF <- x@primerF[i]
    newR <- x@primerR[i]
    PrimerPairsSet(as.character(newF), as.character(newR))
})

setMethod("[", "PairedReadFileSet", function(x, i){
    newF <- x@readsF[i]
    newR <- x@readsR[i]
    PairedReadFileSet(newF, newR)
})

##' Convenient subsetting for MultiAmplicon objects
##'
##' 
##' @title subset MultiAmplicon
##' @param MultiAmplicon-class object
##' @param i numeric, logical or names vector for subsetting rows (==
##'     amplicons)
##' @param j numeric, logical or names vector for subsetting columns
##'     (== read files =~ samples)
##' @return
##' @author Emanuel Heitlinger
setMethod("[", "MultiAmplicon", function(x, i=TRUE, j=TRUE){
    newPrimer <- x@PrimerPairsSet[i]
    newFiles <- x@PairedReadFileSet[j]
    if(any(dim(x@rawCounts)>0)){
        newRC <- rawCounts(x)[i,j]
    } else{newRC <- matrix()}
    new("MultiAmplicon",
        PrimerPairsSet = newPrimer,
        PairedReadFileSet = newFiles,
        rawCounts = newRC
        )
})

