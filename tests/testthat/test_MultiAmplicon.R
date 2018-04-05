
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

## MA1 <- sortAmplicons(MA, max.mismatch=3)

MA1 <- sortAmplicons(MA)

test_that("sortAmplicons resut doesn't change over executions", {
    ## For multi amplicon objects
    expect_known_output(MA1, system.file("testdata", "MA_sorted.rda",
                                         package = "MultiAmplicon"))
})

SA1 <- sortAmplicons(SA)

## this need major debuggin!!!
## getUnmatched(MA1)

context("Do empty files produce empty data?")
test_that("rowCounts is zero for empty file", {
    ## For multi amplicon objects
    expect_equal(colSums(rawCounts(MA1))[["S00_F_filt.fastq.gz"]], 0)
    ## For single amplicon objects
    expect_equal(colSums(rawCounts(SA1))[["S00_F_filt.fastq.gz"]], 0)
})

context("Do nonsensical primers result in zero matches?")
test_that("rowCounts is zero for nonsensical primer", {
    ## For multi amplicon objects
    expect_equal(rowSums(rawCounts(MA1))[["Amp5F.Amp5R"]], 0)
    ## For single amplicon objects ... no empty amplicon used
    ##    expect_equal(rowSums(rawCounts(SA1))[["Amp5F.Amp5R"]], 0)
})

context("SortAmplcion produced two files for each non-empty sample and amplicon?")
                   
## get only non empty samples raw counts
test_that("number of files written equals non-zero samples in rawCounts", {
    ## For multi amplicon objects
    F.files <- unlist(lapply(MA1@stratifiedFiles, function (x) x@readsF))
    R.files <- unlist(lapply(MA1@stratifiedFiles, function (x) x@readsF))
    expect_equal(sum(rawCounts(MA1)>0), length(F.files))
    expect_equal(sum(rawCounts(MA1)>0), length(R.files))
    ## For single amplicon objects
    F.files <- unlist(lapply(SA1@stratifiedFiles, function (x) x@readsF))
    R.files <- unlist(lapply(SA1@stratifiedFiles, function (x) x@readsF))
    expect_equal(sum(rawCounts(SA1)>0), length(F.files))
    expect_equal(sum(rawCounts(SA1)>0), length(R.files))
})

test_that("N statified files is N of non-zero samples x amplicons",{
    ## For multi amplicon objects
    expect_equal(rowSums(rawCounts(MA1)>0), 
                 unlist(lapply(MA1@stratifiedFiles, length)))
    ## For single amplicon objects
    expect_equal(rowSums(rawCounts(SA1)>0), 
                 unlist(lapply(SA1@stratifiedFiles, length)))
})

test_that("files for each amplicon contain the number of reads reported", {
    ## For multi amplicon objects
    expect_equivalent(
        unlist(lapply(MA1@stratifiedFiles,
                      function (x) length(readFastq(x@readsR)))),
        rowSums(rawCounts(MA1)))
    ## For single amplicon objects
    expect_equivalent(
        unlist(lapply(SA1@stratifiedFiles,
                      function (x) length(readFastq(x@readsR)))),
        rowSums(rawCounts(SA1)))
})

context("SortAmplcion can be made less stringent?")
test_that("less stringent sorting results in more reads accepted", {
    expect_true(all(sortAmplicons(MA, countOnly=TRUE, starting.at=1:2) >=
                    rawCounts(MA1)))
    expect_true(all(sortAmplicons(MA, countOnly=TRUE, max.mismatch=1) >=
                   rawCounts(MA1)))
})


MA2 <- derepMulti(MA1, mc.cores=1)

context("Dereplication works?")
test_that("dereplication produces a list of derep objects ", {
    expect_equal(length(MA2@derep), nrow(MA2))
    expect_equal(length(MA2@derep), nrow(rawCounts(MA2)))
    expect_equal(unlist(lapply(MA2@derep, length)),
                 rowSums(rawCounts(MA2)>1)) # >1 singl seq rm
})

up1.counts <- t(rawCounts(MA2))[t(rawCounts(MA2))>1] # >1 singl seq rm

up1.dereps <- unname(unlist(lapply(MA2@derep, function (x){
    lapply(x, function (y) sum(slot(y, "derepF")$uniques))
})))

test_that("all sequences are dereplicated ", {
expect_equal(up1.counts, up1.dereps)
})


MA3 <- dadaMulti(MA2, err=NULL, selfConsist=TRUE, pool=FALSE, 
                 multithread=TRUE)

context("Denoising works?")
test_that("dada2 denoising produces a list of derep objects ", {
    expect_equal(length(MA3@derep), nrow(MA3))
    expect_equal(length(MA3@derep), nrow(rawCounts(MA3)))
    expect_equal(unlist(lapply(MA3@derep, length)),
                 rowSums(rawCounts(MA3)>1)) # >1 singl seq rm
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
    expect_equal(length(MA4@mergers), nrow(rawCounts(MA3)))
    expect_equal(unlist(lapply(MA4@mergers, length)),
                 rowSums(rawCounts(MA3)>1)) # >1 singl seq rm
})

up1.merge <- unname(unlist(lapply(MA4@mergers, function (x)
    lapply(x, function (y) sum(getUniques(y))))))

## this is only expected when concantenating during merge, otherwise
## non-merging sequences would be removed
test_that("all sequences are dereplicated ", {
    expect_equal(up1.counts, up1.merge)
})

MA5 <- sequenceTableMulti(MA4)

context("Sequence table has entries corresponding to stratified files ")
test_that("stratified files result in the number of columns of sequence tables ",
{
    expect_true(all(unlist(lapply(MA5@stratifiedFiles, length)) ==
                    unlist(lapply(MA5@sequenceTable, nrow)) |
                    unlist(lapply(MA5@stratifiedFiles, length)) ==
                    unlist(lapply(MA5@sequenceTable, nrow))+1))
    ## last case for if a single sequence was dropped derep object
})

MA6 <- noChimeMulti(MA5)

context("Merging works?")
test_that("merging produces a list of derep objects ", {
    expect_equal(length(MA4@mergers), nrow(MA3))
    expect_equal(length(MA4@mergers), nrow(rawCounts(MA3)))
    expect_equal(unlist(lapply(MA4@mergers, length)),
                 rowSums(rawCounts(MA3)>1)) # >1 singl seq rm
})


## context("Subsetting MultiAmplicon objects")
## ToDO: proper tests from the below... 

MA6[1:2, 1:2]

MA6[2, 4]@rawCounts

MA6[2, 6]@rawCounts

MA6[2, 6]@derep

MA6[2, 6]@dada

MA6[2, 6]@stratifiedFiles

MA6[c(4, 6),  c(4, 6)]@rawCounts

MA6[c(4, 6),  c(4, 6)]@derep

## killed the bug!
MA6[1:6, 1:7]


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




