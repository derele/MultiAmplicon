### Subsetting methods

## already documented in accessors @param x PrimerPairsSet-class
## object
##' @param i integer, or logical value which primer pairs to select or
##'     character string giving the name of the primer pair
##' @param j not used
##' @param ... not used
##' @param drop not used
## ##' @importClassesFrom Matrix index 
##' @rdname PrimerPairsSet-class
setMethod("[", c("PrimerPairsSet", "index", "missing", "ANY"),
          function(x, i, j, ..., drop=NA){
              if(class(i)=="logical"){
                  i <- logical.to.numeric(i, length(x))
              }
              if(class(i)=="character"){
                  i <- name.to.numeric(i, names(x))
              }
              newF <- x@primerF[i]
              newR <- x@primerR[i]
              PrimerPairsSet(primerF=as.character(newF),
                             primerR=as.character(newR))
          })

##' @param x PairedReadFileSet-class object
##' @param i numeric to select
##' @param j not used
##' @param ... not used
##' @param drop not used
##' @rdname PairedReadFileSet-class
setMethod("[", c("PairedReadFileSet", "index", "missing", "ANY"),
          function(x, i, j, ..., drop=TRUE){
              if(class(i)=="logical"){
                  i <- logical.to.numeric(i, length(x))
              }
              if(class(i)=="character"){
                  i <- name.to.numeric(i, names(x))
              }
              newF <- x@readsF[i]
              newR <- x@readsR[i]
              PairedReadFileSet(newF, newR)
          })

##' @param x PairedDerep-class
##' @param i index or logical indicating which sequencing read file
##'     pair to select or character string giving its name
##' @param j not used
##' @param ... not used
##' @param drop not used
##' @rdname PairedDerep-class
setMethod("[", c("PairedDerep", "index", "missing", "ANY"),
          function(x, i, j, ..., drop=TRUE){
              if(class(i)=="logical"){
                  i <- logical.to.numeric(i, length(x))
              }
              if(class(i)=="character"){
                  i <- name.to.numeric(i, names(x))
              }
              newF <- slot(x, "derepF")[i]
              newR <- slot(x, "derepR")[i]
              new("PairedDerep",
                  derepF=newF, derepR=newR)
          })

##' @param x PairedDada-class object
##' @param i numeric to select
##' @param j not used
##' @param ... not used
##' @param drop not used
##' @rdname PairedDada-class
setMethod("[", c("PairedDada", "index", "missing", "ANY"),
          function(x, i, j, ..., drop=TRUE){
              if(class(i)=="logical"){
                  i <- logical.to.numeric(i, length(x))
              }
              if(class(i)=="character"){
                  i <- name.to.numeric(i, names(x))
              }
              newF <- slot(x, "dadaF")[i]
              newR <- slot(x, "dadaR")[i]
              ## if we get a single object, but  need a list of them
              if(class(newF)%in%"dada"){
                  newF <- list(newF)
                  newR <- list(newR)
              }
              new("PairedDada",
                  dadaF=newF, dadaR=newR)
          })

## ##' @param x PairedDada-class object
## ##' @param i integer to select
## ##' @param j integer to select
## ##' @param ... not used
## ##' @param drop not used
## ##' @rdname PairedDada-class
## setMethod("[", c("PairedDada", "integer", "integer", "ANY"),
##           function(x, i, j, ..., drop=TRUE){
##     newF <- lapply(x@dadaF[i], "[", j)
##     newR <- lapply(x@dadaR[i], "[", j)
##     new("PairedDada",
##         dadaF=newF, dadaR=newR)
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
##' @importFrom assertthat not_empty
##' @rdname MultiAmplicon-class
setMethod("[", c("MultiAmplicon", "index", "index", "ANY"),
          function(x, i, j, ..., drop=FALSE){
              newPrimer <- x@PrimerPairsSet[i]
              suppressWarnings( ## to avoid validity messages
                  newFiles <- x@PairedReadFileSet[j]
              )
              newRC <- matrix()
              newSF <- list()
              newderep <- list()
              newdada <- list()
              newmergers <- list()
              newST <- list()
              newSTnC <- list()
              if(class(i)!=class(j)){
                  stop("both indices should be off the same class")
              }
              if(class(i)=="logical"){
                  i <- logical.to.numeric(i, nrow(x))
                  j <- logical.to.numeric(j, ncol(x))
              }
              if(class(i)=="character"){
                  i <- name.to.numeric(i, rownames(x))
                  j <- name.to.numeric(j, colnames(x))
              }
              if(not_empty(x@rawCounts)){
                  newRC <- as.matrix(getRawCounts(x)[i, j, drop=FALSE])
                  ## we drop empty files from statified files
                  ## therefore we have to find new indices j. These
                  ## later have to be used also for the columns of
                  ## sequence tables.
                  new.j <- lapply(seq_along(i), function (ii) {
                      zero.i <- which(getRawCounts(x)[i[[ii]], ]>1) # >1 singl seq rm
                      which(zero.i%in%j)
                  })
                  newSF <-  lapply(seq_along(i), function (ii) {
                      x@stratifiedFiles[[i[[ii]]]][new.j[[ii]]]
                  })
                  names(newSF) <- names(x@stratifiedFiles[i])
              }
              if(length(x@derep)>0){
                  newderep <- lapply(seq_along(i), function (ii){
                      x@derep[[i[[ii]]]][new.j[[ii]]]
                  })
                  names(newderep) <- names(x@derep[i])
              }
              if(length(x@dada)>0){
                  newdada <- lapply(seq_along(i), function (ii){
                      x@dada[[i[[ii]]]][new.j[[ii]]]
                  })
                  names(newdada) <- names(x@dada[i])
              }
              if(length(x@mergers)>0){
                  newmergers <- lapply(seq_along(i), function (ii){
                      x@mergers[[i[[ii]]]][new.j[[ii]]]
                  })
                  names(newmergers) <- names(x@mergers[i])
              }
              if(length(x@sequenceTable)>0){
                  newST <- lapply(seq_along(i), function (ii){
                      x@sequenceTable[[i[[ii]]]][new.j[[ii]], , drop=FALSE]
                  })
                  names(newST) <- names(x@sequenceTable[i])
              }
              if(length(x@sequenceTableNoChime)>0){
                  newSTnC <- lapply(seq_along(i), function (ii){
                      x@sequenceTableNoChime[[i[[ii]]]][new.j[[ii]], , drop=FALSE]
                  })
                  names(newSTnC) <- names(x@sequenceTableNoChime[i])         
              }              
              ## avoid warnings from ValidityCheck
              suppressWarnings(
                  MA.out <- initialize(x,
                                       PrimerPairsSet = newPrimer,
                                       PairedReadFileSet = newFiles,
                                       rawCounts = newRC,
                                       stratifiedFiles = newSF,
                                       derep = newderep,
                                       dada = newdada,
                                       mergers = newmergers,
                                       sequenceTable = newST,
                                       sequenceTableNoChime = newSTnC)
              )
              if(length(x@sequenceTableFilled)>0){           
                  fillSampleTables(MA.out, samples="union")
              } else {MA.out}
          })

## empty column
##' @rdname MultiAmplicon-class
setMethod("[", c("MultiAmplicon", "index", "missing", "ANY"),
          function(x, i, j, ..., drop=FALSE){
          x[i, 1:ncol(x)]    
          }
)

## empty row
##' @rdname MultiAmplicon-class
setMethod("[", c("MultiAmplicon", "missing", "index", "ANY"),
          function(x, i, j, ..., drop=FALSE){
          x[1:nrow(x), j]    
          }
)

## all empty
##' @rdname MultiAmplicon-class
setMethod("[", c("MultiAmplicon", "missing", "missing", "ANY"),
          function(x, i, j, ..., drop=FALSE){
          x
          }
)


logical.to.numeric <- function(x, n){
    stop("only numeric/ingeger indices are currently supported for indexing")
         ## if (length(x)!=n){
         ##     warning("logical subscript is recycled to match number of ",
         ##             "slots in object")}
         ## if (n%%length(x)!=0) {
         ##     stop("number of logical indices (", 
         ##          length(i),
         ##          ") is not multiple of slots in object (",
         ##          n, ")")
         ## }
         ## index_i <- rep_len(i, length.out=n)
         ## return(index_i)
    }

name.to.numeric <- function(x, names){
        stop("only numeric/ingeger indices are currently supported for indexing")
##     which(names%in%x)
}
