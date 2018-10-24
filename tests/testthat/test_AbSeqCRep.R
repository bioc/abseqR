context("AbSeqCRep S4 class")

## Wed Oct 24 12:08:01 AEDT 2018
## This test file will test
##      1. abseqReport
##      2. The instantiation of AbSeqCRep using the "+" operator

# Load example data -- ex -------------------------------------------------
tmpdir <- tempdir()
exdata <- system.file("extdata", "ex", package = "abseqR")
file.copy(exdata, tmpdir, recursive = TRUE)
ex <- file.path(tmpdir, "ex")


# Test S4 class -----------------------------------------------------------
test_that("AbSeqCRep is correctly instantiated", {
    samples <- abseqReport(ex, report = 0)
    expect_length(samples, 3)
    expect_named(samples,
                 c("PCR1", "PCR2", "PCR3"),
                 ignore.order = FALSE,
                 ignore.case = FALSE)
    csamples <- Reduce("+", samples)
    expect_true(isS4(csamples))
    expect_equal(class(csamples)[[1]], "AbSeqCRep")
})

test_that("AbSeqCRep is correctly instantiated when 3 samples are combined", {
    samples <- abseqReport(ex, report = 0)
    expect_length(samples, 3)
    expect_named(samples,
                 c("PCR1", "PCR2", "PCR3"),
                 ignore.order = FALSE,
                 ignore.case = FALSE)
    csamples <- Reduce("+", samples)
    expect_length(.asRepertoireList(csamples), 3)
})

test_that("AbSeqCRep is correctly instantiated during partial addition", {
    samples <- abseqReport(ex, report = 0)
    expect_length(samples, 3)
    pcr1_pcr2 <- samples[["PCR1"]] + samples[["PCR3"]]
    expect_length(.asRepertoireList(pcr1_pcr2), 2)
    expect_equal(
        unlist(lapply(.asRepertoireList(pcr1_pcr2), .asRepertoireName)),
        c("PCR1", "PCR3"),
        ignore.order = FALSE,
        ignore.case = FALSE)
})
