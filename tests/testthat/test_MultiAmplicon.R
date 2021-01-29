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

SA1 <- sortAmplicons(SA, filedir=tempfile())

MA1 <- sortAmplicons(MA, filedir=tempfile())

context("Does read sorting work?")

test_that("no reads were sorted into different samples" , {
    ## Test read confusion in the read sorting. Think about making this
    ## available to the user with a warning that it's compute intesive
    trackReadSorting <- function (MA) { 
        readsFL <- lapply(seq_along(MA@PairedReadFileSet), function(i) {
            rawFiles <- MA@PairedReadFileSet@readsF[[i]]
            rawFiles <- rawFiles[which(file.exists(rawFiles))]
            ShortRead::readFastq(rawFiles)
        })
        snames <- colnames(MA)
        names(readsFL) <- snames
        sort_track <- lapply(snames, function (sampl) { 
            strat <- lapply(getStratifiedFilesF(MA), function(x) {
                grep(sampl, x, value=TRUE)
            })
            readsF_stratified <- ShortRead::readFastq(unlist(strat))
            IDsStrat <- ShortRead::id(readsF_stratified)
            readsF <- readsFL[[sampl]]        
            IDsRaw <- ShortRead::id(readsF)
            cbind(RawReads = length(readsF),
                  UniqueRawReads = length(unique(IDsRaw)),
                  SortedReads = length(readsF_stratified),
                  UniqueSortedReads = length(unique(IDsStrat)),
                  Problems = length(IDsStrat[which(!as.vector(IDsStrat) %in%
                                                   as.vector(IDsRaw))])
                  )
        })
        sort_track <- as.data.frame(do.call(rbind, sort_track))
        base::rownames(sort_track) <- snames
        sort_track$DuplicateRaw <- sort_track$RawReads -
            sort_track$UniqueRawReads
        sort_track$DuplicateSorted  <- sort_track$SortedReads -
            sort_track$UniqueSortedReads
        sort_track
    }

    sortingStats <- trackReadSorting(MA1)
    cat("\n\n Read sorting statistics\n")
    print(sortingStats)
    cat("\n\n")
    ## has the table been generated with the proper column
    expect_gt(length(sortingStats$Problems), 0)
    ## are there no sorting problems
    expect_equal(sum(sortingStats$Problems), 0)
    ## and are thery the sam for replicated samples
    expect_equal(unname(t(sortingStats["S05_F_filt.fastq.gz", ])),
                 unname(t(sortingStats["S05D_F_filt.fastq.gz", ])))
})


test_that("rowCounts is zero for empty file", {
    ## For multi amplicon objects
    expect_equal(colSums(getRawCounts(MA1))[["S00_F_filt.fastq.gz"]], 0)
    ## For single amplicon objects
    expect_equal(colSums(getRawCounts(SA1))[["S00_F_filt.fastq.gz"]], 0)
})

test_that("rowCounts is zero for nonsensical primer", {
    ## For multi amplicon objects
    expect_equal(rowSums(getRawCounts(MA1))[["Amp5F.Amp5R"]], 0)
    ## For single amplicon objects ... no empty amplicon used
    ##    expect_equal(rowSums(getRawCounts(SA1))[["Amp5F.Amp5R"]], 0)
})

## get only non empty samples raw counts
test_that("number of files written equals non-zero samples in rawCounts", {
    ## For multi amplicon objects
    F.files <- unlist(getStratifiedFilesF(MA1))
    R.files <- unlist(getStratifiedFilesR(MA1))
    expect_equal(sum(getRawCounts(MA1)>0), length(F.files))
    expect_equal(sum(getRawCounts(MA1)>0), length(R.files))
    ## For single amplicon objects
    F.files <- unlist(getStratifiedFilesF(SA1))
    R.files <- unlist(getStratifiedFilesR(SA1))
    expect_equal(sum(getRawCounts(SA1)>0), length(F.files))
    expect_equal(sum(getRawCounts(SA1)>0), length(R.files))
})

test_that("files for each amplicon contain the number of reads reported", {
    ## For multi amplicon objects
    expect_equivalent(
        lapply(seq(nrow(MA1)), function (i){
            length(ShortRead::readFastq(getStratifiedFilesF(MA1[i, ])))
        }), 
        rowSums(getRawCounts(MA1)))
})

test_that("files for each sample contain the number of reads reported", {
    ## For multi amplicon objects
    expect_equivalent(
        lapply(seq(ncol(MA1)), function (i){
            length(ShortRead::readFastq(getStratifiedFilesF(MA1[, i])))
        }), 
        colSums(getRawCounts(MA1)))
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


MA3 <- dadaMulti(MA1, selfConsist=TRUE, pool=FALSE, 
                 multithread=TRUE)

context("Denoising works and doesn't confuse samples?")

test_that("dada2 denoising produces a list of dada objects ", {
    expect_equal(length(MA3@dada), nrow(MA3))
    expect_equal(length(MA3@dada), nrow(getRawCounts(MA3)))
    expect_equal(unlist(lapply(MA3@dada, length)),
                 rowSums(getRawCounts(MA3)>0)) # >1 singl seq rm
})


test_that("dada2 denoising produces identical results for replicate samplesd", {
    expect_equal(lapply(getDadaF(MA3[, "S05_F_filt.fastq.gz"]),
                        unname),
                 ## have to unname the amplicon naming
                 lapply(getDadaF(MA3[, "S05D_F_filt.fastq.gz"]),
                        unname))
})


up1.dadas <- unname(unlist(lapply(MA3@dada, function (x)
    lapply(slot(x, "dadaF"), function (y) sum(getUniques(y))))))


## test_that("all sequences are dereplicated ", {
##     expect_equal(up1.counts, up1.dadas)
## })

MA4 <- mergeMulti(MA3, justConcatenate=c(TRUE, TRUE),
                  verbose=FALSE, maxMismatch = c(15, 20, 18))

context("Merging works?")
test_that("merging produces a list of derep objects ", {
    expect_equal(length(MA4@mergers), nrow(MA3))
    expect_equal(length(MA4@mergers), nrow(getRawCounts(MA3)))
    expect_equal(unlist(lapply(MA4@mergers, length)),
                 rowSums(getRawCounts(MA3)>0)) # >1 singl seq rm
})

## for some weird reason this fails (only on TravisCI) NO IDEA WHY
context("Merging works?")
test_that("proportion of merged is between zero and one ", {
    expect_true(
    all((calcPropMerged(MA4) >= 0 & calcPropMerged(MA4) <= 1))
    )
})



up1.merge <- unname(unlist(lapply(MA4@mergers, function (x)
    lapply(x, function (y) sum(getUniques(y))))))


MA5 <- makeSequenceTableMulti(MA4)

test_that("Identical files produce identical sequence tables ", {
    lapply(MA5@sequenceTable, function(x) {
        if("S05_F_filt.fastq.gz" %in% base::rownames(x)){
            expect_equal(
                unname(t(x["S05_F_filt.fastq.gz", ])),
                unname(t(x["S05D_F_filt.fastq.gz", ])))
        }
    })
})
              

MA6 <- removeChimeraMulti(MA5)


test_that("Identical files produce identical NoChime sequence tables ", {
    lapply(MA6@sequenceTableNoChime, function(x) {
        if("S05_F_filt.fastq.gz" %in% base::rownames(x)){
            expect_equal(
                unname(t(x["S05_F_filt.fastq.gz", ])),
                unname(t(x["S05D_F_filt.fastq.gz", ])))
        }
    })
})


test_that("Reads in sequence tables map to stratified files", {
    ### these funcitons are candicates for making them available in
    ### the package itself
    mapReadsStratTab <- function(MA) {
        getReadsBySample <- function(MA){
            sreads <- lapply(colnames(MA), function (sampl) { 
                strat <- lapply(getStratifiedFilesF(MA), function(x) {
                    grep(sampl, x, value=TRUE)
                })
                ShortRead::readFastq(unlist(strat))
            })
            names(sreads) <- colnames(MA)
            sreads
        }
        RbyS <- getReadsBySample(MA)
        RbyS <- lapply(RbyS, ShortRead::sread)
        
        seqtabL <- getSequenceTableNoChime(MA)
        seqtabL <- seqtabL[unlist(lapply(seqtabL, function(x) all(dim(x)>0)))]

        SbyS <- lapply(colnames(MA), function(sampl){
            sbys <- lapply(seqtabL, function(ST) {
                sampl.here <- sampl[sampl%in%base::rownames(ST)]
                cn <- base::colnames(ST)[which(ST[sampl.here,]>0)]
                split <- strsplit(cn, "NNNNNNNNNN")
                unlist(lapply(split, "[[", 1))
                ## it might be necessary to also track the counts here
            })
            sbys[!unlist(lapply(sbys, is.null))]
        })
        names(SbyS) <- colnames(MA)

        SbyS <- SbyS[intersect(names(SbyS), names(RbyS))]
        RbyS <- RbyS[intersect(names(SbyS), names(RbyS))]

        sapply(names(SbyS), function (na) {
            unlist(SbyS[[na]]) %in% unique(as.vector(RbyS[[na]]))
        })
    }
    ## Map the reads
    map <- mapReadsStratTab(MA6)
    ## and execute the testing
    expect_true(all(unlist(lapply(map, all))))
    expect_true(any(unlist(lapply(map, length))>0))
})


context("Subsetting MultiAmplicon objects")

test_that("subsetting leaves rawCounts intact", {
    expect_equal(getRawCounts(MA6[2, 6]), getRawCounts(MA6)[2, 6, drop=FALSE])
    expect_equal(getRawCounts(MA6[3:4, 2:5]), getRawCounts(MA6)[3:4, 2:5, drop=FALSE])
    })




### THIS analyses SAMPLE CONFUSION caused by resorting before dada
resortedMA3 <- dadaMulti(MA1[, c(6:4,1L,3L,2L, 8L, 7L)],
                  selfConsist=TRUE, pool=FALSE, multithread=TRUE)

## in all other steps resorting does not produce an error
resortedMA4 <- mergeMulti(resortedMA3, justConcatenate=TRUE)
resortedMA5 <- makeSequenceTableMulti(resortedMA4)
resortedMA6 <- removeChimeraMulti(resortedMA5)

## trackReadSorting(resortedMA6)

################# EVALUATE sorting #############

test_that("Resorting produces identical output over samples", {
    seqtabs <- getSequenceTableNoChime(MA6)
    Sorttabs <- getSequenceTableNoChime(resortedMA6)
    SamSums <- lapply(seqtabs, base::rowSums)
    SortSums <- lapply(Sorttabs, base::rowSums)
    samorder <- colnames(MA6)
    ## over samples
    confusion <- lapply(base::colnames(SamSums), function(name) {
        df <- cbind(SamSums[[name]][samorder], SortSums[[name]][samorder])
        base::rownames(df) <- samorder
        base::colnames(df) <- c("Correct", "Shuffle")
        df[is.na(df)] <- 0
        df
    })
    ## evaluate
    lapply(confusion, function(x) {
        expect_equal(x[, "Correct"], x[, "Shuffle"])
    })
})

## context("concatenating MA objects")

## ## this doesn't work because of this
## ##  .concatenateStratifiedFiles(MA[, 1:4]@stratifiedFiles, MA[,
## ##  5:8]@stratifiedFiles) Error in x[[i]] (from allClasses.R#48) :
## ##  subscript out of bounds

## ## a solution would be a final push to make statified files a
## ## matrix probably!!!?

## test_that("Concatenating over samples works", {
##     expect_equal(concatenateMultiAmplicon(MA[, 1:4], MA[, 5:8]), 
##    MA)
## })

context("adding sample information to MA objects")

additionalSD <- data.frame(whatever=letters[seq(ncol(MA6))],
                           row.names=rownames(MA6@sampleData))

test_that("sample data with exact matches of rowname works", {
    expect_s4_class(addSampleData(MA, additionalSD), "MultiAmplicon")
})

additionalSD <- data.frame(whatever=letters[seq(ncol(MA6)+2)],
                           row.names=c(rownames(MA6@sampleData), "foo", "bar"))

test_that("sample data warning when more samples than sequence data", {
    expect_warning(addSampleData(MA, additionalSD),
                   "sampleData but have no sequence data reported\\. They will be omitted")
})

additionalSD <- data.frame(whatever=letters[1:4],
                           row.names=rownames(MA6@sampleData)[1:4])

test_that("sample data warning when more sequence data than sample data", {
    expect_warning(addSampleData(MA, additionalSD),
                   "missing from your sampleData but seem to have sequence data reported")
})




context("Handing over to Phyloseq")

test_that("toPhyloseq multi2Single TRUE/FALSE work and give same resultsw ", {
    PHYWO <- toPhyloseq(MA6, samples=colnames(MA6))
    PHYWO.list <- toPhyloseq(MA6, samples=colnames(MA6), multi2Single=FALSE)
    expect_s4_class(PHYWO, "phyloseq")
    expect_true(all(unlist(lapply(PHYWO.list, class))%in%c("phyloseq", "NULL")))
    PHYWO.list <- PHYWO.list[!unlist(lapply(PHYWO.list, is.null))]
    all.samplesL <- lapply(PHYWO.list, function (x) rownames(otu_table(x)))
    all.samples <- rownames(otu_table(PHYWO))
    lapply(all.samplesL, function (x) expect_equal(x, all.samples))
    all.seqL <- unname(unlist(lapply(PHYWO.list, function (x) colnames(otu_table(x)))))
    all.seq <- colnames(otu_table(PHYWO))
    expect_equal(all.seqL, all.seq)
})

context("Get the pipeline summary")
test_that("pipelin  Summary is a data.frame", {
    expect_s3_class(getPipelineSummary(MA6), "data.frame")
    expect_s3_class(getPipelineSummary(MA3), "data.frame")
})


