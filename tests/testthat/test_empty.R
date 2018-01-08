context("Do empty files produce empty data")

## running the whole example MA class example first to get MA1 data
## loaded
test_example("../../man/sortAmplicons.Rd")

test_that("rowCounts is zero", {
    expect_equal(colSums(rawCounts(MA1))[["S00_F_filt.fastq.gz"]], 0)
})
