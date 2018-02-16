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
    newF <- x@derepF[i]
    newR <- x@derepR[i]
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
##' @param ... not used
##' @param drop should not be used
##' @rdname MultiAmplicon-class
setMethod("[", c("MultiAmplicon", "integer", "integer", "ANY"),
          function(x, i, j, ..., drop=TRUE){
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
                  ## we drop empty files from statified files
                  ## therefore we have to find new indices j. These
                  ## later have to be used also for the columns of
                  ## sequence tables.
                  new.j <- lapply(seq_along(i), function (ii) {
                      zero.i <- which(x@rawCounts[ii, ]>0)
                      which(zero.i%in%j)
                  })
                  newSF <-  lapply(seq_along(i), function (ii) {
                      x@stratifiedFiles[[ii]][new.j[[ii]]]
                  })
              }
              if(length(x@derep)>0){
                  newderep <- lapply(seq_along(i), function (ii){
                      x@derep[[ii]][new.j[[ii]]]
                  })
              }
              if(length(x@dada)>0){
                  newdada <- lapply(seq_along(i), function (ii){
                      x@dada[[ii]][new.j[[ii]]]
                  })
              }
              if(length(x@mergers)>0){
                  newmergers <- lapply(seq_along(i), function (ii){
                      x@mergers[[ii]][new.j[[ii]]]
                  })
              }
              if(length(x@sequenceTable)>0){
                  newST <- lapply(seq_along(i), function (ii){
                      x@sequenceTable[[ii]][new.j[[ii]], ]
                  })
              }
              if(length(x@sequenceTableNoChime)>0){
                  newSTnC <- lapply(seq_along(i), function (ii){
                      x@sequenceTableNoChime[[ii]][new.j[[ii]], ]
                  })
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
