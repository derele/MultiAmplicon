setMethod("[", c("PrimerPairsSet", "ANY"),  function(x, i){
    newF <- x@primerF[i]
    newR <- x@primerR[i]
    PrimerPairsSet(as.character(newF), as.character(newR))
})

setMethod("[", c("PairedReadFileSet", "ANY"), function(x, i){
    newF <- x@readsF[i]
    newR <- x@readsR[i]
    PairedReadFileSet(newF, newR)
})

##' Convenient subsetting for MultiAmplicon objects
##'
##' Subset a MultiAmplicon object including all potentially filled
##' slots
##' 
##' @title subset MultiAmplicon
##' @param MultiAmplicon-class object
##' @param i numeric, logical or names vector for subsetting rows (==
##'     amplicons)
##' @param j numeric, logical or names vector for subsetting columns
##'     (== read files, corresponding usually to samples)
##' @param drop should 
##' @return a subset of the original MultiAmplicon object
##' @export
##' @author Emanuel Heitlinger
setMethod("[", c("MultiAmplicon", "ANY", "ANY", "logical"),
          function(x, i=TRUE, j=TRUE, drop=FALSE){
              newPrimer <- x@PrimerPairsSet[i, drop=drop]
              newFiles <- x@PairedReadFileSet[j, drop=drop]
              if(any(dim(x@rawCounts)>0)){
                  newRC <- as.matrix(rawCounts(x)[i, j, drop=drop])
##                   newSF <- x@stratifiedFiles[i] ## FIXME wrong
                  ## because j are
                  ## dropped without
                  ## placeholder!!
              } else{newRC <- matrix()}
              new("MultiAmplicon",
                  PrimerPairsSet = newPrimer,
                  PairedReadFileSet = newFiles,
                  rawCounts = newRC #,
##                  stratifiedFiles = newSF
                  )
          })

