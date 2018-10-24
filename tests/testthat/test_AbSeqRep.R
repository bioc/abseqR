context("AbSeqRep S4 class")

## Wed Oct 24 12:08:01 AEDT 2018
## This test file will test
##      1. (indirectly) accessors
##      2. abseqReport
##      3. The instantiation of AbSeqRep

# Load example data -- ex -------------------------------------------------
tmpdir <- tempdir()
exdata <- system.file("extdata", "ex", package = "abseqR")
file.copy(exdata, tmpdir, recursive = TRUE)
ex <- file.path(tmpdir, "ex")


# Test S4 class -----------------------------------------------------------
test_that("AbSeqRep is correctly instantiated with all 3 samples", {
    samples <- abseqReport(ex, report = 0)
    expect_length(samples, 3)
    expect_named(samples,
                 c("PCR1", "PCR2", "PCR3"),
                 ignore.order = FALSE,
                 ignore.case = FALSE)
})


test_that("AbSeqRep object was loaded correctly using its accessors", {
    samples <- abseqReport(ex, report = 0)
    expectedResults <- list(
        .asRepertoireName = c("PCR1", "PCR2", "PCR3"),
        .asRepertoireDir = rep("ex", 3),
        .asRepertoireChain = rep("hv", 3),
        .asRepertoireBitscore = rep("300 - Inf", 3),
        .asRepertoireQueryStart = rep("1 - Inf", 3),
        .asRepertoireAlignLen = rep("250 - Inf", 3),
        .asRepertoireSubjectStart = rep("1 - 3", 3),
        .asRepertoirePrimer5 = rep("None", 3),
        .asRepertoirePrimer3 = rep("None", 3),
        .asRepertoireUpstream = rep("None", 3)
    )
    lapply(seq_along(expectedResults), function(i) {
        functor <- names(expectedResults)[[i]]
        obs <- unlist(lapply(samples, function(x) {
            basename(do.call(functor, list(x)))
        }))
        exp <- expectedResults[[i]]
        expect_equivalent(obs, exp)
    })
})
