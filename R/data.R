#' Sequencing run files for the carnivore biome dataset
#'
#' A dataset containing SRA submission, details of the SRA dataset for
#' BioProject PRJNA386767
#' "Intestinal biome sequencing of carnivores". This data was
#' generated in R using the following code:
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/bioproject/PRJNA386767}
#' 
#' @examples
#' ## showing how this was generated.
#' \dontrun{sqlfile <- "SRAmetadb.sqlite"
#' if(!file.exists(sqlfile)){
#'     sqlfile <- getSRAdbFile()
#' }
#' sra_con <- dbConnect(SQLite(),sqlfile)
#' carnivoreSeqRuns <- getSRA(search_terms = '"Intestinal biome sequencing of carnivores"',
#'                           sra_con=sra_con, acc_only=TRUE)
#' }
#' 
"carnivoreSeqRuns"
