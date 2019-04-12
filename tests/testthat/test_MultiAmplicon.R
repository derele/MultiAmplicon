## test file for use whith testthat
primerF <- c(Amp1F = "AGAGTTTGATCCTGGCTCAG", Amp2F = "ACTCCTACGGGAGGCAGC",
             Amp3F = "GAATTGACGGAAGGGCACC", Amp4F = "YGGTGRTGCATGGCCGYT",
             Amp5F = "AAAAACCCCGGGGGGTTTTT", Amp6F = "AGAGTTTGATCCTGCCTCAG")
primerR <- c(Amp1R = "CTGCWGCCNCCCGTAGG", Amp2R = "GACTACHVGGGTATCTAATCC",
             Amp3R = "AAGGGCATCACAGACCTGTTAT", Amp4R = "TCCTTCTGCAGGTTCACCTAC",
             Amp5R = "AAAAACCCCGGGGGGTTTTT", Amp6R = "CCTACGGGTGGCAGATGCAG")
PPS <- PrimerPairsSet(primerF, primerR)
fastq.dir <- system.file("extdata", "fastq", package = "MultiAmplicon")
fastq.files <- list.files(fastq.dir, full.names=TRUE)
Ffastq.file <- fastq.files[grepl("F_filt", fastq.files)]
Rfastq.file <- fastq.files[grepl("R_filt", fastq.files)]
PRF <- PairedReadFileSet(Ffastq.file, Rfastq.file)

MA <- MultiAmplicon(PPS, PRF)
SA <- MultiAmplicon(PrimerPairsSet(primerF[1], primerR[1]), PRF)

context("sortAmplicons resut doesn't change results over executions")

MA1 <- sortAmplicons(MA, filedir=tempfile())

test_that("sortAmplicons resut doesn't change over executions", {
    ## For multi amplicon objects
    tmp <- tempfile()
    expect_known_output(MA1, tmp, print=TRUE)
})

SA1 <- sortAmplicons(SA, filedir=tempfile())

## this need major debuggin!!!
## getUnmatched(MA1)

context("Do empty files produce empty data?")
test_that("rowCounts is zero for empty file", {
    ## For multi amplicon objects
    expect_equal(colSums(getRawCounts(MA1))[["S00_F_filt.fastq.gz"]], 0)
    ## For single amplicon objects
    expect_equal(colSums(getRawCounts(SA1))[["S00_F_filt.fastq.gz"]], 0)
})

context("Do nonsensical primers result in zero matches?")
test_that("rowCounts is zero for nonsensical primer", {
    ## For multi amplicon objects
    expect_equal(rowSums(getRawCounts(MA1))[["Amp5F.Amp5R"]], 0)
    ## For single amplicon objects ... no empty amplicon used
    ##    expect_equal(rowSums(getRawCounts(SA1))[["Amp5F.Amp5R"]], 0)
})

context("SortAmplcion produced two files for each non-empty sample and amplicon?")
                   
## get only non empty samples raw counts
test_that("number of files written equals non-zero samples in rawCounts", {
    ## For multi amplicon objects
    F.files <- unlist(lapply(MA1@stratifiedFiles, function (x) x@readsF))
    R.files <- unlist(lapply(MA1@stratifiedFiles, function (x) x@readsF))
    expect_equal(sum(getRawCounts(MA1)>0), length(F.files))
    expect_equal(sum(getRawCounts(MA1)>0), length(R.files))
    ## For single amplicon objects
    F.files <- unlist(lapply(SA1@stratifiedFiles, function (x) x@readsF))
    R.files <- unlist(lapply(SA1@stratifiedFiles, function (x) x@readsF))
    expect_equal(sum(getRawCounts(SA1)>0), length(F.files))
    expect_equal(sum(getRawCounts(SA1)>0), length(R.files))
})

test_that("N statified files is N of non-zero samples x amplicons",{
    ## For multi amplicon objects
    expect_equal(rowSums(getRawCounts(MA1)>0), 
                 unlist(lapply(MA1@stratifiedFiles, length)))
    ## For single amplicon objects
    expect_equal(rowSums(getRawCounts(SA1)>0), 
                 unlist(lapply(SA1@stratifiedFiles, length)))
})

test_that("files for each amplicon contain the number of reads reported", {
    ## For multi amplicon objects
    expect_equivalent(
        unlist(lapply(MA1@stratifiedFiles,
                      function (x) length(readFastq(x@readsR)))),
        rowSums(getRawCounts(MA1)))
    ## For single amplicon objects
    expect_equivalent(
        unlist(lapply(SA1@stratifiedFiles,
                      function (x) length(readFastq(x@readsR)))),
        rowSums(getRawCounts(SA1)))
})

context("SortAmplcion can be made less stringent?")
test_that("less stringent sorting results in more reads accepted", {
    expect_true(all(sortAmplicons(MA, filedir=tempfile(), countOnly=TRUE,
                                  starting.at=1:2) >=
                    getRawCounts(MA1)))
    expect_true(all(sortAmplicons(MA, filedir=tempfile(), countOnly=TRUE,
                                  max.mismatch=1) >=
                   getRawCounts(MA1)))
})


MA2 <- derepMulti(MA1, mc.cores=1)

context("Dereplication works?")
test_that("dereplication produces a list of derep objects ", {
    expect_equal(length(MA2@derep), nrow(MA2))
    expect_equal(length(MA2@derep), nrow(getRawCounts(MA2)))
    expect_equal(unlist(lapply(MA2@derep, length)),
                 rowSums(getRawCounts(MA2)>1)) # >1 singl seq rm
})

up1.counts <- t(getRawCounts(MA2))[t(getRawCounts(MA2))>1] # >1 singl seq rm

up1.dereps <- unname(unlist(lapply(MA2@derep, function (x){
    lapply(x, function (y) sum(slot(y, "derepF")$uniques))
})))

test_that("all sequences are dereplicated ", {
expect_equal(up1.counts, up1.dereps)
})

## set OMEGA_C to avoid removing sequences
MA3 <- dadaMulti(MA2, selfConsist=TRUE, pool=FALSE, OMEGA_C=0, 
                 multithread=TRUE)

MA3.R <- dadaMulti(MA2, selfConsist=TRUE, pool=FALSE, 
                   multithread=TRUE)

context("Denoising works?")
test_that("dada2 denoising produces a list of derep objects ", {
    expect_equal(length(MA3@dada), nrow(MA3))
    expect_equal(length(MA3@dada), nrow(getRawCounts(MA3)))
    expect_equal(unlist(lapply(MA3@dada, length)),
                 rowSums(getRawCounts(MA3)>1)) # >1 singl seq rm
})

up1.dadas <- unname(unlist(lapply(MA3@dada, function (x)
    lapply(slot(x, "dadaF"), function (y) sum(getUniques(y))))))


test_that("all sequences are dereplicated ", {
    expect_equal(up1.counts, up1.dadas)
})

MA4 <- mergeMulti(MA3, justConcatenate=c(TRUE, TRUE),
                  verbose=FALSE, maxMismatch = c(15, 20, 18))

context("Merging works?")
test_that("merging produces a list of derep objects ", {
    expect_equal(length(MA4@mergers), nrow(MA3))
    expect_equal(length(MA4@mergers), nrow(getRawCounts(MA3)))
    expect_equal(unlist(lapply(MA4@mergers, length)),
                 rowSums(getRawCounts(MA3)>1)) # >1 singl seq rm
})

context("Merging works?")
test_that("proportion of merged is between zero and one ", {
    expect_true(
    all((calcPropMerged(MA4) >= 0 & calcPropMerged(MA4) <= 1))
    )
})

up1.merge <- unname(unlist(lapply(MA4@mergers, function (x)
    lapply(x, function (y) sum(getUniques(y))))))

## this is only expected when concantenating during merge, otherwise
## non-merging sequences would be removed
test_that("all sequences are dereplicated ", {
    expect_equal(up1.counts, up1.merge)
})

MA5 <- makeSequenceTableMulti(MA4)

context("Sequence table has entries corresponding to stratified files ")
test_that("stratified files result in the number of columns of sequence tables ",
{
    expect_true(all(unlist(lapply(MA5@stratifiedFiles, length)) ==
                    unlist(lapply(MA5@sequenceTable, nrow)) |
                    unlist(lapply(MA5@stratifiedFiles, length)) ==
                    unlist(lapply(MA5@sequenceTable, nrow))+1))
    ## last case for if a single sequence was dropped derep object
})

MA6 <- removeChimeraMulti(MA5)

context("Merging works?")
test_that("merging produces a list of derep objects ", {
    expect_equal(length(MA4@mergers), nrow(MA3))
    expect_equal(length(MA4@mergers), nrow(getRawCounts(MA3)))
    expect_equal(unlist(lapply(MA4@mergers, length)),
                 rowSums(getRawCounts(MA3)>1)) # >1 singl seq rm
})

context("Subsetting MultiAmplicon objects")

test_that("subsetting leaves rawCounts intact", {
    expect_equal(getRawCounts(MA6[2, 6]), getRawCounts(MA6)[2, 6, drop=FALSE])
    expect_equal(getRawCounts(MA6[3:4, 2:5]), getRawCounts(MA6)[3:4, 2:5, drop=FALSE])
    })

MA1.alt <- sortAmplicons(MA1[2:5, 2:6], filedir=tempfile())
MA2.alt <- derepMulti(MA1.alt)

test_that("sorting a subsetted object same as subsetting a sorted object", {
   expect_equal(MA2.alt@derep, MA2[2:5, 2:6]@derep)
   expect_equal(MA2.alt@rawCounts, MA2[2:5, 2:6]@rawCounts)
   ## expect_equal(MA2.alt, MA2[2:5, 2:6])
   ## expect_equal(MA2.alt@stratifiedFiles, MA2[2:5, 2:6]@stratifiedFiles)
})

## Thought I had a bug in stratified file subsetting before, but this
## can't work as subsetting creates different stratified file names by
## default
## expect_equal(getStratifiedFilesF(MA2.alt), getStratifiedFilesF(MA2[2:5, 2:6]))
## expect_equal(getStratifiedFilesR(MA2.alt), getStratifiedFilesR(MA2[2:5, 2:6]))


context("Subsetting and concatenating MultiAmplicon objects")

foo  <- concatenateMultiAmplicon((MA6[, 1:3]), MA6[, 4:7])

## Again this fails, despite the fact that there is no new
## stratification involved!  Hove to FIX THIS!!!
## test_that("subsetting and concatenation go hand in hand", {
##     expect_equal(MA6@stratifiedFiles, foo@stratifiedFiles)
## })

Test_that("subsetting and concatenation go hand in hand", {
    expect_equal(MA6@rawCounts, foo@rawCounts)
})

test_that("subsetting and concatenation go hand in hand", {
    expect_equal(MA6@dada, foo@dada)
})

test_that("subsetting and concatenation go hand in hand", {
    expect_equal(MA6@derep, foo@derep)
})

## Again this fails, the problem being likely the one empty amplicon
## mess...HAVE TO FIX THIS!!!

## test_that("subsetting and concatenation
## go hand in hand", { expect_equal(MA6@sequenceTable,
## foo@sequenceTable) })

## test_that("subsetting and concatenation
## go hand in hand", { expect_equal(MA6@sequenceTableNoChime,
## foo@sequenceTableNoChime) })


MA7 <- getBlastTaxAnnot(MA6)

test_that("subsetting and concatenation
go hand in hand", { expect_equal(MA6@taxonTable,
foo@taxonTable) })


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




