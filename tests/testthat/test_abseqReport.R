context("abseqReport function loader and plotter")

## Wed Oct 24 12:08:01 AEDT 2018
## This test file will test
##      1. abseqReport and report = 1

# Load example data -- ex -------------------------------------------------
tmpdir <- tempdir()
exdata <- system.file("extdata", "ex", package = "abseqR")
file.copy(exdata, tmpdir, recursive = TRUE)
ex <- file.path(tmpdir, "ex")


# Helper functions --------------------------------------------------------
removeSample <- function(path, sampleNames) {

    stopifnot(
        all(
            sampleNames %in%
                lapply(list.dirs(path = file.path(path, RESULT_DIR),
                                 recursive = FALSE),
                       basename)))

    analysisParamPaths <-
        file.path(path, RESULT_DIR, sampleNames, ANALYSIS_PARAMS)

    contents <- lapply(analysisParamPaths, function(app) {
        readLines(app)
    })

    unlink(analysisParamPaths)

    contents
}

addSample <- function(path, sampleNames, contents) {
    stopifnot(
        all(
            sampleNames %in%
                lapply(list.dirs(path = file.path(path, RESULT_DIR),
                                 recursive = FALSE),
                       basename)))


    analysisParamPaths <-
        file.path(path, RESULT_DIR, sampleNames, ANALYSIS_PARAMS)

    lapply(seq_along(analysisParamPaths), function(i) {
        app <- analysisParamPaths[[i]]
        content <- contents[[i]]

        fileConn <- file(app)
        writeLines(content, fileConn)
        close(fileConn)
    })
}


# Test abseqReport function -----------------------------------------------
test_that("abseqReport generates a message that the directory has been moved", {
    expect_message(abseqReport(ex, report = 0), regex = ".*different.*move.*")
})


## TODO: test on report = 2, 3
test_that("abseqReport doesn't fall over when reporting", {
    # we're only testing one sample because it will take way to long otherwise
    remove <- c("PCR1", "PCR3")
    buffer <- removeSample(ex, remove)
    samples <- abseqReport(ex, report = 1)
    expect_length(samples, 1)
    expect_named(samples, "PCR2")
    # not strictly necessary to do this, but just in case there's another
    # test function that wants to re-use the removed samples for analysis
    addSample(ex, remove, buffer)
})
