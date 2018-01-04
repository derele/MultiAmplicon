### Subsetting methods

##' @param x PrimerPairsSet-class object 
##' @param i numeric to select
##' @rdname PrimerPairsSet-class
setMethod("[", c("PrimerPairsSet", "ANY"),  function(x, i){
    newF <- x@primerF[i]
    newR <- x@primerR[i]
    PrimerPairsSet(primerF=as.character(newF), primerR=as.character(newR))
})

##' @param x PairedReadFileSet-class object
##' @param i numeric to select
##' @rdname PairedReadFileSet-class
setMethod("[", c("PairedReadFileSet", "ANY"), function(x, i){
    newF <- x@readsF[i]
    newR <- x@readsR[i]
    PairedReadFileSet(newF, newR)
})
##' @param x PairedDerep-class
##' @param i numeric to select
##' @rdname PairedDerep-class
setMethod("[", c("PairedDerep", "ANY"), function(x, i){
    newF <- x@derepF[i]
    newR <- x@derepR[i]
    new("PairedDerep",
        derepF=newF, derepR=newR)
})

##' @param x PairedDada-class object
##' @param i numeric to select
##' @rdname PairedDada-class
setMethod("[", c("PairedDada", "ANY"), function(x, i){
    newF <- x@dadaF[i]
    newR <- x@dadaR[i]
    new("PairedDada",
        dadaF=newF, dadaR=newR)
})

##' Convenient subsetting for MultiAmplicon objects
##'
##' Subset a MultiAmplicon object including all potentially filled
##' slots
##' 
##' @title subset MultiAmplicon
##' @param x MultiAmplicon-class object
##' @param i numeric, logical or names vector for subsetting rows (==
##'     amplicons)
##' @param j numeric, logical or names vector for subsetting columns
##'     (== read files, corresponding usually to samples)
##' @param drop should not be used
##' @rdname MultiAmplicon-class
setMethod("[", "MultiAmplicon",
          function(x, i=TRUE, j=TRUE, drop="missing"){
              newPrimer <- x@PrimerPairsSet[i]
              newFiles <- x@PairedReadFileSet[j]
              newRC <- matrix()
              newSF <- list()
              newderep <- list()
              newdada <- list()
              newmergers <- list()
              newST <- list()
              newSTnC <- list()
              if(all(dim(x@rawCounts)>0)){
                  newRC <- as.matrix(rawCounts(x)[i, j])
                  ## we drop empty files (i.e. also empty "samples" in
                  ## derep, dada, etc objects) without a placeholder,
                  for(ii in i){
                      ## therefore we need to find new, fixed j indices
                      zero.indices <- which(x@rawCounts[ii, ]>0)
                      newJ <- which(zero.indices%in%j)
                      ## create updated objects for existing slots
                      if(length(x@stratifiedFiles)>0){
                          newSF[[length(newSF)+1]] <-
                              lapply(x@stratifiedFiles[ii], "[", newJ)
                      }
                      if(length(x@derep)>0){
                          newderep[[length(newderep)+1]] <-
                              lapply(x@derep[ii], "[", newJ)
                      }
                      if(length(x@dada)>0){
                          newdada[[length(newdada)+1]] <-
                              lapply(x@dada[ii], "[", newJ)
                      }
                      if(length(x@mergers)>0){
                          newmergers[[length(newST)+1]] <-
                              lapply(x@mergers[ii], "[", newJ)
                      }
                      if(length(x@sequenceTable)>0){
                          newST[[length(newST)+1]] <-
                              lapply(x@sequenceTable[ii], function (y){
                                  y[newJ, ]
                              })
                      }
                      if(length(x@sequenceTableNoChime)>0){
                          newSTnC[[length(newSTnC)+1]] <-
                              lapply(x@sequenceTableNoChime[ii], function(y){
                                  y[newJ, ]
                              })
                      }
                  }
              }
          initialize(x,
              PrimerPairsSet = newPrimer,
              PairedReadFileSet = newFiles,
              rawCounts = newRC,
              stratifiedFiles = newSF,
              derep = newderep,
              dada = newdada,
              mergers = newmergers,
              sequenceTable = newST,
              sequenceTableNoChime = newSTnC
              )
})

