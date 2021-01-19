#' Download sequencing run files for a carnivore biome dataset
#'
#' For use in the [real-world data
#' vignette](https://derele.github.io/MultiAmplicon/articles/MultiAmplicon-real-world-example.html)),
#' we download a dataset containing details of the SRA dataset for
#' BioProject PRJNA386767
#' "Intestinal biome sequencing of carnivores". We download the sample
#' data from SRA, re-format it for our later analysis. And use the SRA
#' run ids to download the raw sequencing files for precessing in the
#' vignette.
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/bioproject/PRJNA386767}
#'
#' @examples
#' ## This data was generated in R using the following code:
#' \dontrun{
#' ## create an SRA database connection
#' library(SRAdb)
#' sqlfile <- "SRAmetadb.sqlite"
#' if(!file.exists(sqlfile)){
#'     sqlfile <- getSRAdbFile()
#'
#' sra_con <- dbConnect(SQLite(),sqlfile)
#'
#' carnivoreSeqRuns <- getSRA(search_terms = '"Intestinal biome sequencing of carnivores"',
#'                            sra_con=sra_con, acc_only=FALSE)
#'
#' carnivoreSeqRuns <- carnivoreSeqRuns[, c('run','study','sample',
#'                                        'experiment', 'sample_attribute')]
#'
#' sample.vars <- strsplit(carnivoreSeqRuns$sample_attribute, " \\|\\| ")
#'
#' sample.variable <- lapply(sample.vars, function(x){
#'     what <- strsplit(x, ": ")
#'     variable <- lapply(what, "[[", 2)
#'     names(variable) <- lapply(what, "[[", 1)
#'     variable
#' })
#'
#' carnivoreSeqRuns <- cbind(carnivoreSeqRuns, do.call(rbind, sample.variable))
#' 
#' carnivoreSeqRuns$sample_attribute <- NULL
#'
#' runs <- carnivoreSeqRuns$run
#' 
#' destDir <- "download_sra"
#' 
#' if(!file_test("-d", destDir)) dir.create(destDir)
#' 
#' fastqFiles <- list.files(destDir, pattern=".fastq.gz", full.names=TRUE)
#' 
#' if(!length(fastqFiles)){
#'     getSRAfile(runs, sra_con, fileType = 'fastq' , srcType = "ftp",
#'                destDir=destDir)
#'     fastqFiles <- list.files(destDir, pattern=".fastq.gz", full.names=TRUE)
#'   }
#' }
#' }
"carnivoreSeqRuns"
