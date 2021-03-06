% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-methods.R
\name{plotAmpliconNumbers}
\alias{plotAmpliconNumbers}
\alias{plotAmpliconNumbers,matrix-method}
\alias{plotAmpliconNumbers,MultiAmplicon-method}
\alias{plotPipelineSummary,MultiAmplicon-method}
\alias{plotPipelineSummary,data.frame-method}
\title{plotAmpliconNumbers}
\usage{
plotAmpliconNumbers(MA, transf = function(x) log10(x + 1), ...)

\S4method{plotAmpliconNumbers}{matrix}(MA, transf = function(x) log10(x + 1), ...)

\S4method{plotAmpliconNumbers}{MultiAmplicon}(MA, transf = function(x) log10(x + 1), ...)

\S4method{plotPipelineSummary}{MultiAmplicon}(MA)

\S4method{plotPipelineSummary}{data.frame}(MA)
}
\arguments{
\item{MA}{The primer matching (amplicon) matrix of a
MultiAmplicon object or the object itself}

\item{transf}{transformation to be applied, recommended default is
log10(x+1)}

\item{...}{additional parameter to be passed to pheatmap function}
}
\value{
just like the original \code{\link{pheatmap}} function, A
    heatmap is drawn to the graphics output device. Returned are
    (invisibly) a list of components: 1."tree_row" the clustering
    of rows as hclust object 2. "tree_col" the clustering of
    columns as hclust object 3. "kmeans" the kmeans clustering of
    rows if parameter kmeans_k was specified.
}
\description{
Plot the raw counts for primer pairs matched in different sample
files as a heatmap
}
\details{
The function uses pheatmap to plot the primer matching statistics
stored in the rawCounts slot of a MultiAmplicons-class object
}
\author{
Emanuel Heitlinger
}
