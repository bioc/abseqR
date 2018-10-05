#' Visualize all analysis conducted by abseqPy
#'
#' @description Plots all samples in the output directory supplied to abseqPy's
#' \code{--outdir} or \code{-o} argument. Users can optionally specify
#' which samples in \code{directory} should be compared. Doing so generates
#' additional plots for clonotype comparison and the usual plots will also
#' conveniently include these samples using additional \code{aes}thetics.
#'
#' Calling this function with a valid \code{directory} will always return a
#' named \code{list} of objects; these individual objects can be
#' combined using the \code{+} operator to form a new comparison, in which the
#' \link{report} function accepts as its first parameter.
#'
#' @include util.R
#' @include AbSeqRep.R
#' @import BiocParallel
#'
#' @param directory string type. directory as specified in
#' \code{-o} or \code{--outdir} in abseqPy. This tells AbSeq where to look for
#' abseqPy's output.
#' @param report (optional) integer type. The possible values are:
#' \itemize{
#'   \item{0 - does nothing (returns named list of \linkS4class{AbSeqRep} objects)}
#'   \item{1 - generates plots for csv files}
#'   \item{2 - generates a report that collates all plots}
#'   \item{3 - generates interactive plots in report (default)}
#' }
#' each higher value also does what the previous values do. For example, \code{report = 2}
#' will return a named list of \linkS4class{AbSeqRep} objects, plot csv files,
#' and generate a (non-interactive)HTML report that collates all the plots together.
#' @param compare (optional) vector of strings. From the samples in found in \code{directory}
#' directory, they can be selected and compared against each other. For example,
#'  to compare "sample1" with "sample2" and "sample3" with "sample4",
#'  \code{compare} should be c("sample1,sample2", "sample3,sample4"). An error
#'  will be thrown if the samples specified in this parameter are not found in
#'  \code{directory}.
#' @param BPPARAM (optional) BiocParallel backend. Configures the parallel implementation.
#' Refer to \href{https://bioconductor.org/packages/release/bioc/html/BiocParallel.html}{BiocParallel}
#' for more information. By default, use all available cores.
#'
#' @return named list. List of \linkS4class{AbSeqRep} objects. The names of
#' the list elements are taken directly from the repertoire object itself.
#' This return value is consistent with the return value of \code{\link{report}}.
#'
#' @seealso \linkS4class{AbSeqRep}
#'
#' @export
#'
#' @seealso \code{\link{report}}. Analogous function, but takes input from
#' an \linkS4class{AbSeqRep} or \linkS4class{AbSeqCRep} object instead.
#'
#' @examples
#' # Use example data from abseqR as abseqPy's output, substitute this
#' # with your own abseqPy output directory
#' abseqPyOutput <- tempdir()
#' file.copy(system.file("extdata", "ex", package = "abseqR"), abseqPyOutput, recursive=TRUE)
#'
#' ### 1. The `report` parameter usage example:
#'
#' # report = 0; don't plot, don't collate a HTML report, don't show anything interactive
#' samples <- abseqReport(file.path(abseqPyOutput, "ex"), report = 0)
#' # samples is now a named list of AbSeqRep objects
#'
#' # report = 1; just plot pngs; don't collate a HTML report; nothing interactive
#' # samples <- abseqReport(file.path(abseqPyOutput, "ex"), report = 1)
#' # samples is now a named list of AbSeqRep objects
#'
#' # report = 2; plot pngs; collate a HTML report; HTML report will NOT be interactive
#' # samples <- abseqReport(file.path(abseqPyOutput, "ex"), report = 2)
#' # samples is now a named list of AbSeqRep objects
#'
#' # report = 3 (default); plot pngs; collate a HTML report; HTML report will be interactive
#' # samples <- abseqReport(file.path(abseqPyOutput, "ex"), report = 3)
#' # samples is now a named list of AbSeqRep objects
#'
#' ### 2. Using the return value of abseqReport:
#'
#' # NOTE, often, this is used to load multiple samples from different directories
#' # using abseqReport (with report = 0), then the samples are added together
#' # before calling the report function. This is most useful when the samples
#' # live in different abseqPy output directory.
#'
#' # Note that the provided example data has PCR1, PCR2, and PCR3
#' # samples contained within the directory
#' stopifnot(names(samples) == c("PCR1", "PCR2", "PCR3"))
#'
#' # as a hypothetical example, say we found something
#' # interesting in PCR1 and PCR3, and we want to isolate them:
#' # we want to explicitly compare PCR1 with PCR3
#' pcr13 <- samples[["PCR1"]] + samples[["PCR3"]]
#'
#' # see abseqR::report for more information.
#' # abseqR::report(pcr13)      # uncomment this line to run
#'
#' ### BPPARAM usage:
#'
#' # 4 core machine, use all cores -  use whatever value that suits you
#' nproc <- 4
#' # samples <- abseqReport(file.path(abseqPyOutput, "ex"),
#' #                        BPPARAM = BiocParallel::MulticoreParam(nproc))
#'
#'
#' # run sequentially - no multiprocessing
#' # samples <- abseqReport(file.path(abseqPyOutput, "ex"),
#' #                        BPPARAM = BiocParallel::SerialParam())
#'
#'# see https://bioconductor.org/packages/release/bioc/html/BiocParallel.html
#'# for more information about how to use BPPARAM and BiocParallel in general.
abseqReport <- function(directory, report, compare, BPPARAM) {
    #  ------ sanitize function arguments ---------
    root <- normalizePath(directory)
    if (!(all(c(RESULT_DIR, AUX_DIR) %in% list.files(root)))) {
        stop(paste("Expected to find", RESULT_DIR, "and", AUX_DIR,
                "in", root,
                "but they are missing. This directory should be the output",
                "directory as specified in abseqPy. Aborting."))
    }

    if (missing(report)) {
        report <- 3
    }

    if (missing(BPPARAM)) {
        BPPARAM <- BiocParallel::bpparam()
    }

    if (missing(compare)) {
        # if user didn't specify which to compare, don't compare, just load
        # all samples in root/RESULT_DIR, EXCLUDING sample comparisons
        compare <- .findRepertoires(file.path(root, RESULT_DIR))
    } else {
        # make sure samples in "compare" actually exists
        availSamples <- .findRepertoires(file.path(root, RESULT_DIR))
        lapply(compare, function(s) {
            user.sample.name <- unlist(lapply(strsplit(s, ","), trimws))
            lapply(user.sample.name, function(si) {
                if (!(si %in% availSamples)) {
                    stop(paste("Sample", si, "in 'compare' argument cannot",
                            "be found in", file.path(root, RESULT_DIR)))
                }
            })
        })
        compare <- unique(c(availSamples, compare))
    }

    # TODO: sanitize 'report' properly using constants/functions
    stopifnot(report %in% c(0, 1, 2, 3))
    # convert number to meaningful variables
    if (report == 0) {
        report <- FALSE
        interactivePlot <- FALSE
        loop <- FALSE
    } else if (report == 1) {
        report <- FALSE
        interactivePlot <- FALSE
        loop <- TRUE
    } else if (report == 2) {
        report <- TRUE
        interactivePlot <- FALSE
        loop <- TRUE
    } else {
        report <- TRUE
        interactivePlot <- TRUE
        loop <- TRUE
    }

    if (loop) {
        bplapply(compare, function(pair) {
            sampleNames <- unlist(lapply(strsplit(pair, ","), trimws))

            # depending on the number of samples requested to plot, outputDir
            # is either a <sample>_vs_<sample> format or just <sample>
            # meanwhile, samples will either be a AbSeqCRep or just AbSeqRep.
            outputDir <- file.path(root, RESULT_DIR,
                          paste(sampleNames, collapse = "_vs_"))
            samples <- .loadSamplesFromString(sampleNames,
                                              root,
                                              warnMove = FALSE)
            # due to the cache/shared user issue (see issue)
            # https://github.com/rstudio/rmarkdown/issues/499
            # we delay the report generation until AFTER the multiprocessing part
            # has completed, that is, we will generate the report one by one
            # TODO: improve this, this is taking way too long to render
            abseqR::report(samples,
                           outputDir,
                           report = 1)   # always plot, but DO NOT GENERATE REPORT!
        }, BPPARAM = BPPARAM)
    }

    # whether or not we generate the report, always return a list of all
    # the samples found in "root"
    if (report) {
        return(.generateDelayedReport(root, compare, interactivePlot))
    } else {
        listOfSamples <- list()
        for (pair in compare) {
            sampleNames <- unlist(lapply(strsplit(pair, ","), trimws))
            if (length(sampleNames) == 1) {
                outputDir <- file.path(root, RESULT_DIR, sampleNames[1])
                sample <- .loadSamplesFromString(sampleNames, root,
                                                  warnMove = TRUE)
                listOfSamples[[sample@name]] <- sample
            }
        }
        return(listOfSamples)
    }
}


#' Collate all HTML reports into a single directory and cretate an entry
#' \code{index.html} file that redirects to all collated HTML files
#'
#' @import flexdashboard png knitr BiocStyle
#' @importFrom rmarkdown render pandoc_available
#'
#' @param reports list/vector type. Collection of strings that are path(s)
#' to <sample>_report.html
#' @param individualSamples list type. list of AbSeqRep objects. Used to
#' extract filtering information and \% read counts.
#' @param outputDirectory string type. Where should the report be placed.
#'
#' @return Nothing
.collateReports <-
    function(reports,
             individualSamples,
             outputDirectory) {
        message("Collating report into index.html")

        # get relative links for HTML report (instead of absolute paths)
        relLinks <- lapply(reports, function(pth) {
            return(file.path(ABSEQ_NESTED_HTML_DIR, basename(pth)))
        })

        # define a mask to filter reports based on comparative or single report
        multiSampleMask <- grepl(".*_vs_.*", names(reports))

        # define rmarkdown params
        chains <- paste(lapply(individualSamples, function(x) {
            return(x@chain)
        }), collapse = ",")

        rawReads <- paste(lapply(individualSamples, function(x) {
            .readSummary(file.path(x@outdir, RESULT_DIR, x@name),
                         ABSEQ_RAW_READ_COUNT_KEY)
        }), collapse = ",")

        annotReads <- paste(lapply(individualSamples, function(x) {
            .readSummary(file.path(x@outdir, RESULT_DIR, x@name),
                         ABSEQ_ANNOT_READ_COUNT_KEY)
        }), collapse = ",")

        filtReads <- paste(lapply(individualSamples, function(x) {
            .readSummary(file.path(x@outdir, RESULT_DIR, x@name),
                         ABSEQ_FILT_READ_COUNT_KEY)
        }), collapse = ",")

        prodReads <- paste(lapply(individualSamples, function(x) {
            .readSummary(file.path(x@outdir, RESULT_DIR, x@name),
                         ABSEQ_PROD_READ_COUNT_KEY)
        }), collapse = ",")

        filterSplitter <- "?"
        filters <- paste(lapply(individualSamples, function(x) {
            paste(
                paste0("Bitscore: ", paste(x@bitscore, collapse = " - ")),
                paste0("V alignment length: ", paste(x@alignlen, collapse = " - ")),
                paste0("Query start: ", paste(x@qstart, collapse = " - ")),
                paste0("Subject start: ", paste(x@sstart, collapse = " - ")),
                sep = ","
            )
        }), collapse = filterSplitter)

        templateFile <-
            system.file("extdata", "index.Rmd", package = "abseqR")
        renderParams <- list(
            singleSamples = paste(names(reports)[!multiSampleMask], collapse = ","),
            multiSamples = paste(names(reports)[multiSampleMask], collapse = ","),
            singleSampleLinks = paste(relLinks[!multiSampleMask], collapse = ","),
            multiSampleLinks = paste(relLinks[multiSampleMask], collapse = ","),
            chains = chains,
            rawReads = rawReads,
            filtReads = filtReads,
            annotReads = annotReads,
            prodReads = prodReads,
            filters = filters,
            filterSplitter = filterSplitter,
            analysisParams = file.path(
                individualSamples[[1]]@outdir,
                RESULT_DIR,
                individualSamples[[1]]@name,
                ANALYSIS_PARAMS
            )
        )
        if (rmarkdown::pandoc_available()) {
            rmarkdown::render(templateFile,
                              output_dir = outputDirectory,
                              params = renderParams)
        } else {
            warning("Pandoc not found in system, will not create index.html")
        }
    }

#' Generates report for all samples in 'compare'
#'
#' @description This function is needed because we are delaying the generation
#' of reports until after all threads/processes have joined. There's currently
#' an issue with rmarkdown::render() in parallel execution, see:
#' https://github.com/rstudio/rmarkdown/issues/499
#'
#' @include util.R
#' @include AbSeqRep.R
#'
#' @param root string, project root directory.
#' @param compare vector of strings, each string is a comparison defined
#' by the user (assumes that this value has been checked).
#' @param interactivePlot logical, whether or not to plot interactive plotly
#' plots.
#'
#' @return a named list of samples, each an AbSeqRep object found in "root"
.generateDelayedReport <- function(root, compare, interactivePlot) {
    # before we start creating reports, create the report directory
    reportDir <- file.path(root, ABSEQ_HTML_DIR)
    if (!dir.exists(reportDir)) {
        dir.create(reportDir)
    } else {
        warning(paste(reportDir, "found, overriding contents."))
    }

    # create nested HTML directory for all html files except index.html
    nestedHTMLdir <- file.path(reportDir, ABSEQ_NESTED_HTML_DIR)
    if (!dir.exists(nestedHTMLdir)) {
        dir.create(nestedHTMLdir)
    } else {
        warning(paste(nestedHTMLdir, "found, overriding contents."))
    }

    individualSamples <- list()
    individualReports <- list()
    # populate individualSample list with samples for user to browse and
    # create report if asked to.
    for (pair in compare) {
        sampleNames <- unlist(lapply(strsplit(pair, ","), trimws))
        outputDir <- file.path(root, RESULT_DIR,
                               paste(sampleNames, collapse = "_vs_"))
        samples <- .loadSamplesFromString(sampleNames, root,
                                          warnMove = TRUE)
        if (length(sampleNames) == 1) {
            individualSamples[[samples@name]] <- samples
        }
        pth <- .generateReport(
            samples,
            root = outputDir,
            outputDir = file.path(root, ABSEQ_HTML_DIR, ABSEQ_NESTED_HTML_DIR),
            interactivePlot = interactivePlot,
            .indexHTML = file.path("..", "index.html")
        )
        if (!is.na(pth)) {
            individualReports[paste(sampleNames, collapse = "_vs_")] <- pth
        }
    }

    if (length(individualReports)) {
        .collateReports(
            individualReports,
            individualSamples,
            outputDirectory = file.path(root, ABSEQ_HTML_DIR)
        )
    }
    return(individualSamples)
}


#' Loads AbSeqCRep or AbSeqRep objects from a list of sampleNames
#'
#' @include util.R
#' @include AbSeqRep.R
#'
#' @param sampleNames vector, singleton or otherwise
#' @param root string type. root directory
#' @param warnMove logical type. Warning message if the directory has been
#' moved?
#'
#' @return AbSeqRep or AbSeqCRep object depending on sampleNames
.loadSamplesFromString <- function(sampleNames, root, warnMove = TRUE) {
    Reduce("+", lapply(sampleNames, function(sname) {
        sample <- .loadAbSeqRepFromParams(file.path(root, RESULT_DIR,
                                                    sname, ANALYSIS_PARAMS))
        if (suppressWarnings(normalizePath(sample@outdir)) !=
            suppressWarnings(normalizePath(root))) {
            if (warnMove) {
                warning(paste("Sample output directory",
                              sample@outdir,
                              "is different from provided path",
                              root,
                              "assuming directory was moved"))
            }
            sample@outdir <- root
        }
        return(sample)
    }))
}
