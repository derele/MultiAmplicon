##' Plot the raw counts for primer pairs matched in different sample
##' files as a heatmap
##'
##' The function uses pheatmap to plot the primer matching statistics
##' stored in the rawCounts slot of a MultiAmplicons-class object
##' 
##' @title plot_Amplicon_numbers
##' 
##' @param MAmatrix The primer matching (amplicon) matrix of a
##'     MultiAmplicon object.
##' 
##' @param transf transformation to be applied, recommended default is
##'     log10
##' 
##' @param add offset added to allow e.g. log10 transformation of zero
##'     values
##'
##' @importFrom pheatmap pheatmap
##' @param ... addtional parameter to be passed to pheatmap funciton
##' @return just like the original \code{\link{pheatmap}} function, A
##'     heatmap is drawn to the graphics output device. Returned are
##'     (invisibly) a list of components: 1."tree_row" the clustering
##'     of rows as hclust object 2. "tree_col" the clustering of
##'     columns as hclust object 3. "kmeans" the kmeans clustering of
##'     rows if parameter kmeans_k was specified.
##' @export
##' @author Emanuel Heitlinger

plot_Amplicon_numbers <- function (MAmatrix, transf=function(x) log10(x+1), ...){
    if (nrow(MAmatrix) < 2 || ncol(MAmatrix) < 2) {
        stop(cat("No rawCounts found in MultiAmplicon object:
                  Run sortAmplicons to produce a MultiAmplicon with at least
                  two files and two samples"))
    } else {
        ## get the function name for display on the plot
        transf_function <- deparse(transf)[length(deparse(transf))]
        if (is.primitive(transf)){
            transf_function <- gsub('\\.Primitive\\(\\"(.*)\\"\\)', "\\1",
                                    transf_function)
        }
        pheatmap::pheatmap(transf(MAmatrix),
                           main = paste(transf_function,
                                        "transformed",
                                        "read number"),
                           ...)
    }
}
