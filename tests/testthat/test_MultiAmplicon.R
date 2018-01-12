
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
MA1 <- sortAmplicons(MA)
test_that("sortAmplicons resut doesn't change over executions", {
    ## For multi amplicon objects
    expect_known_output(MA, system.file("testdata", "MA_sorted.rda",
                                        package = "MultiAmplicon"))
})


SA1 <- sortAmplicons(SA)


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
    expect_equal(unname(rowSums(rawCounts(MA1)>0)), 
                 unlist(lapply(MA1@stratifiedFiles, length)))
    ## For single amplicon objects
    expect_equal(unname(rowSums(rawCounts(SA1)>0)), 
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


MA2 <- derepMulti(MA1)

MA3 <- dadaMulti(MA2, err=NULL, selfConsist=TRUE,
                 multithread=TRUE)

## bugging here...
## MA4 <- mergeMulti(MA3, justConcatenate=TRUE)

