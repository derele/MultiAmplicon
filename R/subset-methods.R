### Subsetting methods

## already documented in accessors @param x PrimerPairsSet-class
## object
##' @param i integer, or logical value which primer pairs to select or
##'     character string giving the name of the primer pair
##' @param j not used
##' @param ... not used
##' @param drop not used
##' @importClassesFrom Matrix index
##' @rdname PrimerPairsSet-class
setMethod("[", c("PrimerPairsSet", "index", "missing", "ANY"),
          function(x, i, j, ..., drop=NA){
              newF <- x@primerF[i]
              newR <- x@primerR[i]
              na <- names(x)[i]
              PrimerPairsSet(primerF=as.character(newF),
                             primerR=as.character(newR),
                             names=na)
          })

##' @param x PairedReadFileSet-class object
##' @param i numeric to select
##' @param j not used
##' @param ... not used
##' @param drop not used
##' @rdname PairedReadFileSet-class
setMethod("[", c("PairedReadFileSet", "index", "missing", "ANY"),
          function(x, i, j, ..., drop=TRUE){
              newF <- x@readsF[i]
              newR <- x@readsR[i]
              new("PairedReadFileSet",
                  readsF = newF,
                  readsR = newR,
                  names = x@names[i])
          })

## ##' @param x stratifiedFilesMatrix-class
## ##' @param i index, name or logical indicating which sequencing read file
## ##'     pair to select or character string giving its name
## ##' @param j index, name or logical indicating which sequencing read file
## ##'     pair to select or character string giving its name
## ##' @param ... not used
## ##' @param drop not used
## ##' @rdname stratifiedFilesMatrix-class
## ##' @export
## setMethod("[", "stratifiedFilesMatrix", function(x, i, j, ..., drop=TRUE) {
##     m <- callNextMethod()
##     stratifiedFilesMatrix(readsF=x@readsF[m], readsR=x@readsR[m],
##                           ncol=ncol(m), nrow=nrow(m),
##                           dimnames=list(rownames(m), colnames(m)))
## })


##' Subsetting for MultiAmplicon objects should conveniently subset
##' all (potentially) filled slots
##' 
##' @title subset MultiAmplicon
##' @param x MultiAmplicon-class object
##' @param i numeric, logical or names vector for subsetting rows (==
##'     amplicons)
##' @param j numeric, logical or names vector for subsetting columns
##'     (== read files, corresponding usually to samples)
##' @param ... not used
##' @param drop should not be used
##' @rdname MultiAmplicon-class

setMethod("[", "MultiAmplicon",
          function(x, i, j, ..., drop=TRUE){

              if(all(dim(getRawCounts(x)>0))){
                  newRC <- getRawCounts(x)[i, j, drop=FALSE]
                  newSF <- getStratifiedFilesF(x, dropEmpty=FALSE)[i, j, drop=FALSE]
                  newSR <- getStratifiedFilesR(x, dropEmpty=FALSE)[i, j, drop=FALSE]
                  ## we drop empty files from statified files
                  ## therefore we have to find new indices j. These
                  ## later have to be used also for the columns of
                  ## sequence tables.                  
                  new.j <- lapply(seq_along(i), function (ii) {
                      zero.i <- which(getRawCounts(x)[i[[ii]], ]>0) # >1 singl seq rm
                      which(zero.i%in%j)
                  })
              } else {newRC <- newSF <- newSR <- matrix(nrow=0, ncol=0)}
              if(length(getDerepF(x))>0){
                  newderepF <- getDerepF(x, dropEmpty=FALSE)[i, j, drop=FALSE]
                  newderepR <- getDerepR(x, dropEmpty=FALSE)[i, j, drop=FALSE]
              } else {newderepF <- newderepR <- matrix(nrow=0, ncol=0)}
              if(length(getDadaF(x, dropEmpty=FALSE))>0){
                  newdadaF <- getDadaF(x, dropEmpty=FALSE)[i, j, drop=FALSE]
                  newdadaR <- getDadaR(x, dropEmpty=FALSE)[i, j, drop=FALSE]                  
              } else {newdadaF <- newdadaR <- matrix(nrow=0, ncol=0)}
              if(length(getMergers(x))>0){
                  newmergers <- getMergers(x)[i, j, drop=FALSE]
              } else {newmergers <- matrix(nrow=0, ncol=0)}
              if(length(getSequenceTable(x))>0){
                  newST <- lapply(seq_along(i), function (ii){
                      ST <- x@sequenceTable[[i[[ii]]]]
                      if(nrow(ST)>=length(new.j[[ii]])){
                          ST[new.j[[ii]], , drop=FALSE]
                      } else {matrix(ncol=0, nrow=0)}
                  })
                  names(newST) <- names(x@sequenceTable[i])
              } else{newST <- list()}
              if(length(x@sequenceTableNoChime)>0){
                  newSTnC <- lapply(seq_along(i), function (ii){
                      ST <- x@sequenceTableNoChime[[i[[ii]]]]
                      if(nrow(ST)>=length(new.j[[ii]])){
                          ST[new.j[[ii]], , drop=FALSE]
                      } else {matrix(ncol=0, nrow=0)}
                  })
                  names(newSTnC) <- names(x@sequenceTableNoChime[i])         
              } else{newSTnC <- list()}
              if(length(x@taxonTable)>0){
                  newTT <- x@taxonTable[i]
                  names(newTT) <- names(x@taxonTable[i])         
              } else{newTT <- list()}
              ## avoid warnings from ValidityCheck??
              suppressWarnings(
                  MA.out <- MultiAmplicon(
                      ## .Data = m, not this! as
                      ## new indices needed
                      ## components to subset with one index
                      PrimerPairsSet = getPrimerPairsSet(x)[i],
                      PairedReadFileSet = getPairedReadFileSet(x)[j],
                      sampleData = getSampleData(x)[j, , drop=FALSE],

                      ## components that needed  to be tested for presence
                      rawCounts = newRC,
                      stratifiedFilesF = newSF,
                      stratifiedFilesR = newSR,
                      derepF = newderepF,
                      derepR = newderepR,
                      dadaF = newdadaF,
                      dadaR = newdadaR,
                      mergers = newmergers#,

                      ## ## more complex  to subset components
                      ##                      sequenceTable = newST,
                      ##                      sequenceTableNoChime = newSTnC,
                      ##                       taxonTable = newTT
                      ))
              MA.out
          })

## empty column
##' @rdname MultiAmplicon-class
setMethod("[", c("MultiAmplicon", "index", "missing", "ANY"),
          function(x, i, j, ..., drop=FALSE){
          x[i, seq(length(x@PairedReadFileSet))]    
          }
)

## empty row
##' @rdname MultiAmplicon-class
setMethod("[", c("MultiAmplicon", "missing", "index", "ANY"),
          function(x, i, j, ..., drop=FALSE){
          x[seq(length(x@PrimerPairsSet)), j]    
          }
)




