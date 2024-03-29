% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline_track.R
\name{getPipelineSummary}
\alias{getPipelineSummary}
\title{getPipelineSummary}
\usage{
getPipelineSummary(MA)
}
\arguments{
\item{MA}{MultiAmplicon object with all slots filled for tracking.}
}
\value{
a data.frame of sample, unique sequences and sequencing
    reads numbers per amplicon.
}
\description{
Obtain summary data for amplicons run through the MultiAmplicon pipeline.
}
\details{
Get statistics on the number of samples (with read data), the
number of unique sequence variants and the number of reads left
after processing of amplicons in the MultiAmplicon pipeline. In
some steps of the pipeline dada2 performs quality filtering
excluding non-credible sequence variants.
}
\author{
Emanuel Heitlinger
}
