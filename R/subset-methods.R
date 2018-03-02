### Subsetting methods

## already documented in accessors @param x PrimerPairsSet-class
## object
##' @param i integer, or logical value which primer pairs to select or
##'     character string giving the name of the primer pair
##' @param j not used
##' @param ... not used
##' @param drop not used
##' @rdname PrimerPairsSet-class
setMethod("[", c("PrimerPairsSet", "integer", "missing", "ANY"),
          function(x, i, j, ..., drop=NA){
              newF <- x@primerF[as.integer(i)]
              newR <- x@primerR[as.integer(i)]
              PrimerPairsSet(primerF=as.character(newF),
                             primerR=as.character(newR))
})

##' @rdname PrimerPairsSet-class
setMethod("[", c("PrimerPairsSet", "logical", "missing", "ANY"),
          function(x, i, j, ..., drop=NA){
              if (length(i)!=length(x)){
                  warning("logical subscript is recycled to match length of",
                          "PrimerPairsSet")}
              if (length(x)%%length(i)!=0) {
                  stop("length of PrimerPairsSet (", 
                       length(x),
                       ") is not multiple of length of logical subscript (",
                       length(i), ")")
              }
              index <- rep_len(i, length.out=length(x))
              x[which(index)]
          })

##' @rdname PrimerPairsSet-class
setMethod("[", c("PrimerPairsSet", "character", "missing", "ANY"),
          function(x, i, j, ..., drop=NA){
              x[which(names(x)%in%i)]
          })

##' @param x PairedReadFileSet-class object
##' @param i numeric to select
##' @param j not used
##' @param ... not used
##' @param drop not used
##' @rdname PairedReadFileSet-class
setMethod("[", c("PairedReadFileSet", "integer", "missing", "ANY"),
          function(x, i, j, ..., drop=TRUE){
              newF <- x@readsF[i]
              newR <- x@readsR[i]
              PairedReadFileSet(newF, newR)
})

##' @rdname PairedReadFileSet-class
setMethod("[", c("PairedReadFileSet", "logical", "missing", "ANY"),
          function(x, i, j, ..., drop=NA){
              if (length(i)!=length(x)){
                  warning("logical subscript is recycled to match length of",
                          "PairedReadFileSet")}
              if (length(x)%%length(i)!=0) {
                  stop("length of PairedReadFileSet (", 
                       length(x),
                       ") is not multiple of length of logical subscript (",
                       length(i), ")")
              }
              index <- rep_len(i, length.out=length(x))
              x[which(index)]
          })

##' @rdname PairedReadFileSet-class
setMethod("[", c("PairedReadFileSet", "character", "missing", "ANY"),
          function(x, i, j, ..., drop=NA){
              x[which(names(x)%in%i)]
          })

##' @param x PairedDerep-class
##' @param i integer or logical indicating which sequencing read file
##'     pair to select or character string giving its name
##' @param j not used
##' @param ... not used
##' @param drop not used
##' @rdname PairedDerep-class
setMethod("[", c("PairedDerep", "integer", "missing", "ANY"),
          function(x, i, j, ..., drop=TRUE){
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
setMethod("[", c("PairedDada", "integer", "missing", "ANY"),
          function(x, i, j, ..., drop=TRUE){
    newF <- slot(x, "dadaF")[i]
    newR <- slot(x, "dadaR")[i]
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
##' @param ... not used
##' @param drop should not be used
##' @rdname MultiAmplicon-class
setMethod("[", c("MultiAmplicon", "integer", "integer", "ANY"),
          function(x, i, j, ..., drop=FALSE){
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
                  newRC <- as.matrix(rawCounts(x)[i, j, drop=FALSE])
                  ## we drop empty files from statified files
                  ## therefore we have to find new indices j. These
                  ## later have to be used also for the columns of
                  ## sequence tables.
                  new.j <- lapply(seq_along(i), function (ii) {
                      zero.i <- which(x@rawCounts[ii, ]>1) # >1 singl seq rm
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
          }
)

##' @rdname MultiAmplicon-class
setMethod("[", c("MultiAmplicon", "logical", "logical", "ANY"),
          function(x, i, j, ..., drop=NA){
              if (length(i)!=nrow(x)){
                  warning("logical subscript is recycled to match number of ",
                          "amplicons in object")}
              if (length(j)!=ncol(x)){
                  warning("logical subscript is recycled to match number of ",
                          "samples in object")}
              if (nrow(x)%%length(i)!=0) {
                  stop("number of Amplicons (", 
                       nrow(x),
                       ") is not multiple of length of logical subscript (",
                       length(i), ")")
              }
              if (ncol(x)%%length(j)!=0) {
                  stop("number of samples (", 
                       ncol(x),
                       ") is not multiple of length of logical subscript (",
                       length(j), ")")
              }
              index_i <- rep_len(i, length.out=nrow(x))
              index_j <- rep_len(j, length.out=ncol(x))
              x[which(index_i), which(index_j)]
          })

##' @rdname MultiAmplicon-class
setMethod("[", c("MultiAmplicon", "character", "character", "ANY"),
          function(x, i, j, ..., drop=NA){
              x[which(rownames(x)%in%i), which(colnames(x)%in%j)]
          })
