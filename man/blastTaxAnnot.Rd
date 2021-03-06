% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotateTaxa.R
\name{blastTaxAnnot}
\alias{blastTaxAnnot}
\title{blastTaxAnnot}
\usage{
blastTaxAnnot(
  MA,
  db = "nt/nt",
  num_threads = getOption("mc.cores", 1L),
  negative_gilist = system.file("extdata", "uncultured.gi", package = "MultiAmplicon"),
  infasta = paste0(tempfile(), ".fasta"),
  outblast = paste0(tempfile(), ".blt"),
  taxonSQL,
  ...
)
}
\arguments{
\item{MA}{A MultiAmplicon for which taxa shoudl be annotated. This
should contain a \code{sequenceTableNoChime} slot, meaning
that the MultiAmplicon pipeline needs to be followed to that
point.}

\item{db}{The blast database. Either a full path or the path
relative to the you data base directory set in an
environmental variable.}

\item{num_threads}{The number of threads used for the blast
search.}

\item{negative_gilist}{A file containing NCBI GI numbers to
exclude from blast searches.}

\item{infasta}{A fasta file generated for the blast searche, a
temporary file in the respective R temporary folder by
default.}

\item{outblast}{A blast tabular output file generated by the blast
searche, a temporary file in the respective R temporary folder
by default.}

\item{taxonSQL}{An SQL file generated by the package
\code{\link[taxonomizr]{taxonomizr}}.}

\item{...}{String (of other options) passed to blastn (see blastn
-help in the terminal)}
}
\value{
A MultiAmplicon object with the taxonTable slot filled
}
\description{
Annotate amplified sequence variants (ASVs) with taxa labels
derived from a BLAST search.
}
\details{
Based on a BLAST search taxonomic labels are assigned to
ASVs. Currently supported taxonomy levels are c("superkingdom",
"phylum", "class", "order", "family", "genus", "species"). For
each taxonomic level unique best taxa (based on bitscores) are
reported. If no unique best taxon exists at a particular level NA
is returned for this. The function combines evidence for taxnomic
annotations in case of seperate HSPs for (concatenated) non-merged
sequences.
}
\author{
Emanuel Heitlinger
}
