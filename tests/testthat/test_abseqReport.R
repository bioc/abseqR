context("abseqReport function loader and plotter")

## Wed Oct 24 12:08:01 AEDT 2018
## This test file will test
##      1. abseqReport and report = 0 with message output

# Load example data -- ex -------------------------------------------------
tmpdir <- tempdir()
exdata <- system.file("extdata", "ex", package = "abseqR")
file.copy(exdata, tmpdir, recursive = TRUE)
ex <- file.path(tmpdir, "ex")


# Test abseqReport function -----------------------------------------------
test_that("abseqReport generates a message that the directory has been moved", {
    expect_message(abseqReport(ex, report = 0), regex = ".*different.*move.*")
})
