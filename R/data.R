#' Sequencing run files for the carnivore biome dataset
#'
#' A dataset containing SRA submission, details of the SRA dataset for
#' BioProject PRJNA386767
#' "Intestinal biome sequencing of carnivores". 
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/bioproject/PRJNA386767}
#'
#' @examples
#' ## This data was generated in R using the following code:
#' \dontrun{
#' ## create an SRA database connection
#' sqlfile <- "SRAmetadb.sqlite"
#' if(!file.exists(sqlfile)){
#'     sqlfile <- getSRAdbFile()
#' }
#'
#' sra_con <- dbConnect(SQLite(),sqlfile)
#'
#' ## download the samle data from SRA
#' carnivoreSeqRuns <- getSRA(search_terms = '"Intestinal biome sequencing of carnivores"',
#'                            sra_con=sra_con, acc_only=TRUE)
#'
#' ## add sample data stored in sample_attribute
#' carnivoreSeqRuns <- carnivoreSeqRuns[, c('run','study','sample',
#'                                         'experiment', 'sample_attribute')]
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
#' ## this creates the object distributed with the package
#' carnivoreSeqRuns <- cbind(carnivoreSeqRuns, do.call(rbind, sample.variable))
#' 
#' ## The data linked at
#' ##\url{https://svalbard.biologie.hu-berlin.de:443/d/c291b870178f4abc8c59/}
#' ##for download in the vignette #
#' ##\url{https://derele.github.io/MultiAmplicon/articles/MultiAmplicon-real-world-example.html}
#' ##can then be downloadd as follows
#' carnivoreSeqRuns$sample_attribute <- NULL
#'
#' runs <- carnivoreSeqRuns$run
#' 
#' destDir <- "~/download"
#' 
#' if(!file_test("-d", destDir)) dir.create(destDir)
#' fastqFiles <- list.files(destDir, pattern=".fastq.gz", full.names=TRUE)
#' if(!length(fastqFiles)){
#'     getSRAfile(runs, sra_con, fileType = 'fastq' , srcType = "ftp",
#'                destDir=destDir)
#'     fastqFiles <- list.files(destDir, pattern=".fastq.gz", full.names=TRUE)
#' }
#' }
#' 
"carnivoreSeqRuns"
