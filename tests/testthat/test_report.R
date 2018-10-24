context("'report' function on AbSeqCRep and AbSeqRep objects")

## Wed Oct 24 12:08:01 AEDT 2018
## This test file will test
##      1. abseqReport
##      2. report and (parameter report = 1)
##      3. (indirectly) the "+" operator

# Load example data -- ex -------------------------------------------------
tmpdir <- tempdir()
exdata <- system.file("extdata", "ex", package = "abseqR")
file.copy(exdata, tmpdir, recursive = TRUE)
ex <- file.path(tmpdir, "ex")

# Test report function ----------------------------------------------------
test_that("report correctly accepts AbSeqRep objects", {
    exp <- abseqReport(ex, report = 0)
    expect_length(exp, 3)
    obs <- report(exp[[1]], outputDir = tempdir(), report = 0)
    # careful not to use exp[[1]] as obs will be a singleton *named* list
    expect_equal(exp[1], obs)
})

test_that("report correctly accepts AbSeqCRep objects", {
    exp <- abseqReport(ex, report = 0)
    expect_length(exp, 3)
    obs <- report(Reduce("+", exp), outputDir = tempdir(), report = 0)
    expect_equal(exp, obs)
})


## TODO: test on report = 2, 3
test_that("report dosen't fall over when reporting one sample", {
    # we're only testing one sample because it will take way to long otherwise
    samples <- abseqReport(ex, report = 0)
    expect_length(samples, 3)
    expect_named(samples, c("PCR1", "PCR2", "PCR3"))
    obs <- report(samples[["PCR2"]], outputDir = tempdir(), report = 1)
    expect_equal(samples[2], obs)
})
