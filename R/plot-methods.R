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
##'     columns as hclust object 3. kmeans the kmeans clustering of
##'     rows if parameter kmeans_k was specified.
##' @author Emanuel Heitlinger

plot_Amplicon_numbers <- function (MAmatrix, transf=log10, add=0.1, ...){
    if (nrow(MA@rawCounts) < 2 || ncol(MA@rawCounts) < 2) {
        stop(cat("No rawCounts found in MultiAmplicon object:
                  Run sortAmplicons to produce a MultiAmplicon with at least
                  two files and two samples"))
    } else {
        pheatmap(transf(MAmatrix + add), ...)
    }
}
