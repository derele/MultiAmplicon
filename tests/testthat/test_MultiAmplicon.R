
primerF <- c(Amp1F = "AGAGTTTGATCCTGGCTCAG", Amp2F = "ACTCCTACGGGAGGCAGC",
             Amp3F = "GAATTGACGGAAGGGCACC", Amp4F = "YGGTGRTGCATGGCCGYT",
             Amp5F = "AAAAACCCCGGGGGGTTTTT")
primerR <- c(Amp1R = "CTGCWGCCNCCCGTAGG", Amp2R = "GACTACHVGGGTATCTAATCC",
             Amp3R = "AAGGGCATCACAGACCTGTTAT", Amp4R = "TCCTTCTGCAGGTTCACCTAC",
             Amp5R = "AAAAACCCCGGGGGGTTTTT")
PPS <- PrimerPairsSet(primerF, primerR)
fastq.dir <- system.file("extdata", "fastq", package = "MultiAmplicon")
fastq.files <- list.files(fastq.dir, full.names=TRUE)
Ffastq.file <- fastq.files[grepl("F_filt", fastq.files)]
Rfastq.file <- fastq.files[grepl("R_filt", fastq.files)]
PRF <- PairedReadFileSet(Ffastq.file, Rfastq.file)

MA <- MultiAmplicon(PPS, PRF)
MA1 <- sortAmplicons(MA)

context("Do empty files produce empty data?")
test_that("rowCounts is zero for empty file", {
    expect_equal(colSums(rawCounts(MA1))[["S00_F_filt.fastq.gz"]], 0)
})

context("Do nonsensical primers result in zero matches?")
test_that("rowCounts is zero for nonsensical primer", {
    expect_equal(rowSums(rawCounts(MA1))[["Amp5F.Amp5R"]], 0)
})

context("SortAmplcion produced two files for each non-empty sample and amplicon?")
F.files <- unlist(lapply(MA1@stratifiedFiles, function (x) x@readsF))
R.files <- unlist(lapply(MA1@stratifiedFiles, function (x) x@readsF))
                   
## get only non empty samples raw counts
test_that("number of files written equals non-zero samples in rawCounts", {
    expect_equal(sum(rawCounts(MA1)>0), length(F.files))
    expect_equal(sum(rawCounts(MA1)>0), length(R.files))
})

test_that("N statified files is N of non-zero samples x amplicons",{
    expect_equal(unname(rowSums(rawCounts(MA1)>0)), 
                 unlist(lapply(MA1@stratifiedFiles, length)))
})


test_that("files for each amplicon contain the number of reads reported", {
    expect_equivalent(
        unlist(lapply(MA1@stratifiedFiles,
                      function (x) length(readFastq(x@readsR)))),
      rowSums(rawCounts(MA1)))
})


