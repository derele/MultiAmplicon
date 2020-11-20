MA1.alt <- sortAmplicons(MA1[2:5, 2:6], filedir=tempfile())
MA3.alt <- dadaMulti(MA1.alt, selfConsist=TRUE, pool=FALSE, 
                     multithread=TRUE)

test_that("sorting a subsetted object same as subsetting a sorted object", {
   expect_equal(MA3.alt@rawCounts, MA3[2:5, 2:6]@rawCounts)
   ##  expect_equal(MA3.alt, MA3[2:5, 2:6])
   ##  expect_equal(MA3.alt@stratifiedFiles, MA3[2:5, 2:6]@stratifiedFiles)
})

## Thought I had a bug in stratified file subsetting before, but this
## can't work as subsetting creates different stratified file names by
## default

## This shouldn't be the case and needs FIXING!!!

## expect_equal(getStratifiedFilesF(MA3.alt), getStratifiedFilesF(MA3[2:5, 2:6]))
## expect_equal(getStratifiedFilesR(MA2.alt), getStratifiedFilesR(MA2[2:5, 2:6]))


## context("Taxonomy annotation works")

## MA7 <- blastTaxAnnot(MA6,
##                         infasta=system.file("extdata", "in.fasta", package = "MultiAmplicon"),
##                         outblast=system.file("extdata", "out.blt", package = "MultiAmplicon"))

## test_that("sequence and taxon annotation are in same order", {
##     expect_equal(unlist(lapply(MA7@taxonTable, rownames)),
##                  unlist(lapply(MA7@sequenceTableNoChime, colnames)))
##  })

## ## ## This is  more stringent then even only the number
## ## STnoC <- getSequenceTableNoChime(MA7)
## ## nAnnot <- lapply(getTaxonTable(MA7), nrow)
## #### now are they in sync??
## ## cbind(nAnnot, unlist(lapply(STnoC, ncol)))


## ## ## this might be a suitable test for something...
## ## all(unique(unlist(lapply(getSequenceTableNoChime(MA), rownames)))%in%colnames(MA))
## ## all(colnames(MA)%in%unique(unlist(lapply(getSequenceTableNoChime(MA), rownames))))

## context("Subsetting and concatenating MultiAmplicon objects")

## foo  <- concatenateMultiAmplicon((MA7[, 1:3]), MA7[, 4:7])

## ## Again this fails, despite the fact that there is no new
## ## stratification involved!  Hove to FIX THIS!!!
## ## test_that("subsetting and concatenation go hand in hand", {
## ##     expect_equal(MA6@stratifiedFiles, foo@stratifiedFiles)
## ## })

## test_that("subsetting and concatenation go hand in hand", {
##     expect_equal(MA7@rawCounts, foo@rawCounts)
## })

## test_that("subsetting and concatenation go hand in hand", {
##     expect_equal(MA7@dada, foo@dada)
## })

## test_that("subsetting and concatenation go hand in hand", {
##     expect_equal(MA7@derep, foo@derep)
## })

## Again this fails, the problem being likely the one empty amplicon
## mess...HAVE TO FIX THIS!!!

## test_that("subsetting and concatenation
## go hand in hand", { expect_equal(MA7@sequenceTable,
## foo@sequenceTable) })

## test_that("subsetting and concatenation
## go hand in hand", { expect_equal(MA7@sequenceTableNoChime,
## foo@sequenceTableNoChime) })

## test_that("subsetting and concatenation
## go hand in hand", { expect_equal(MA7@taxonTable,
## foo@taxonTable) })

## getPipelineSummary(MA7)

## failing  from here TODO!!!

## MA3.alt <- dadaMulti(MA2.alt, selfConsist=TRUE, pool=FALSE, multithread=TRUE)

## test_that("dada same after subsetting", {
##     expect_equal(MA3.alt@dada, MA3[2:5, 2:6]@dada)
## })

## MA4.alt <- mergeMulti(MA3.alt, justConcatenate=TRUE)

## test_that("mergers same after subsetting", {
##     expect_equal(MA4.alt@mergers, MA4[2:5, 2:6]@mergers)
## })

## MA5.alt <- makeSequenceTableMulti(MA4.alt)

## test_that("sequence table same after subsetting", {
##     expect_equal(MA5.alt@sequenceTable, MA5[2:5, 2:6]@sequenceTable)
## })

## MA6.alt <- removeChimeraMulti(MA5.alt)

## test_that("noChime same after subsetting", {
##     expect_equal(MA6.alt@sequenceTableNoChime, MA5[2:5, 2:6]@sequenceTableNoChime)
## })


## ## logical and name indexing does not work yet

## ## this works because class(1:5) "integer" and class(1:4*2)
## "numeric" by my own error message

## MA6[1:5, 1:4*2]@rawCounts

## ## This works not really ... maybe define index within subset function
## MA6[c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), c(FALSE, TRUE)]@rawCounts

## ## This shouldn't be empty
## MA6[c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), c(FALSE, TRUE)]@dada

## ## This is still buggin
## MA6[c(FALSE, TRUE), TRUE]

## ## this is good
## expect_identical(MA6[TRUE, TRUE], MA6[1:6, 1:7])

## ## this has a names problem 
## ## expect_identical(MA6, MA6[1:6, 1:7])

## ## this works
## MA6[c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
##     TRUE]

## ## this is wrong... without an error!
## foo <- MA6[c("Amp1F.Amp1R","Amp2F.Amp2R"),

## ## lost the sample names in rownames
## lapply(MA6[1:2, 1:2]@sequenceTable, rownames)

## lapply(MA6[1:6, c(1L, 4L)]@sequenceTableNoChime, dim)

## ## here I have them back
## lapply(MA6[1:6, c(1L, 4L, 5L)]@sequenceTableNoChime, rownames)

## ## seems to work
## MA6[1:6, 5L]@derep


## ## colnames(MA6[c(TRUE, FALSE), c(FALSE, TRUE)])




