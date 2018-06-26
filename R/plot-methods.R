##' Plot the raw counts for primer pairs matched in different sample
##' files as a heatmap
##'
##' The function uses pheatmap to plot the primer matching statistics
##' stored in the rawCounts slot of a MultiAmplicons-class object
##' 
##' @title plotAmpliconNumbers
##' 
##' @param MA The primer matching (amplicon) matrix of a
##'     MultiAmplicon object or the object itself
##' 
##' @param transf transformation to be applied, recommended default is
##'     log10(x+1)
##' 
##' @importFrom pheatmap pheatmap
##' @param ... additional parameter to be passed to pheatmap function
##' @return just like the original \code{\link{pheatmap}} function, A
##'     heatmap is drawn to the graphics output device. Returned are
##'     (invisibly) a list of components: 1."tree_row" the clustering
##'     of rows as hclust object 2. "tree_col" the clustering of
##'     columns as hclust object 3. "kmeans" the kmeans clustering of
##'     rows if parameter kmeans_k was specified.
##' @export
##' @author Emanuel Heitlinger
setGeneric(name="plotAmpliconNumbers",
           def=function(MA, transf=function(x) log10(x+1), ...) {
               standardGeneric("plotAmpliconNumbers")
           })
##' @rdname plotAmpliconNumbers
setMethod("plotAmpliconNumbers", c("matrix", "ANY"),
          function (MA, transf=function(x) log10(x+1), ...){
              ## get the function name for display on the plot
              transf_function <- deparse(transf)[length(deparse(transf))]
              if (is.primitive(transf)){
                  transf_function <- gsub('\\.Primitive\\(\\"(.*)\\"\\)', "\\1",
                                          transf_function)
              }
              pheatmap::pheatmap(transf(MA),
                                 main = paste(transf_function,
                                              "transformed number of",
                                              "sequencing reads"),
                                 ...)
          }
)

##' @rdname plotAmpliconNumbers
setMethod("plotAmpliconNumbers", c("MultiAmplicon", "ANY"),
          function (MA, transf=function(x) log10(x+1), ...){
              MA <- getRawCounts(MA)
              if (nrow(MA) < 2 || ncol(MA) < 2) {
                  stop("No rawCounts found in MultiAmplicon object:
                  Run sortAmplicons to produce a MultiAmplicon-object
                  with at least wo files and two samples")}
              else {plotAmpliconNumbers(MA, transf)}
          })

##' Plot summary data for amplicons run through the MultiAmplicon pipeline.
##'
##' Plot statistics on the number of samples (with read data), the
##' number of unique sequence variants and the number of reads left
##' after processing of amplicons in the MultiAmplicon pipeline. In
##' some steps of the pipeline dada2 performs quality filtering
##' excluding non-credible sequence variants.
##' 
##' @title plotPipelineSummary
##' @param MA MultiAmplicon object with all slots filled for tracking.
##' @return a ggplot object
##' @import ggplot2
## ##' @importFrom assertthat assert_that
##' @export
##' @author Emanuel Heitlinger
setGeneric(name="plotPipelineSummary",
           def=function(MA) {
               standardGeneric("plotPipelineSummary")
           })

##' @rdname plotAmpliconNumbers
setMethod("plotPipelineSummary", "MultiAmplicon",
          function(MA){
              track <- getPipelineSummary(MA)
              plotPipelineSummary(track)
})


##' @rdname plotAmpliconNumbers
setMethod("plotPipelineSummary", "data.frame",
          function(MA){
              if(!all(c("pipeStep", "value", "primer")%in%colnames(MA))){
                  stop("please provide either a MultiAmplicon object or ",
                       "a data.frame produced by `getPipelineSummary`")
              }
              ggplot(MA, aes(pipeStep, value, group=primer))+
                  facet_wrap(~what, scales="free_y")+
                  geom_line()+
                  geom_point()+
                  scale_x_discrete("steps in pipeline")+
                  scale_y_continuous("number of (samples with) reads")+
                  theme_bw() +
                  theme(axis.text.x=element_text(angle = -45, hjust = 0))
          })


