#' Visualize all analysis conducted by AbSePy
#'
#' @description Plots all samples in the output directory supplied to abseqPy's
#' \code{--outdir} or \code{-o} argument.
#'
#' @include util.R
#' @include AbSeqRep.R
#' @import BiocParallel
#'
#' @param root string type. Root directory as specified in
#' \code{-o} or \code{--outdir} in abseqPy. This tells AbSeq where to look for
#' abseqPy's output.
#' @param report integer type. The possible values are:
#' \itemize{
#'   \item{0 - does nothing (returns named list of \linkS4class{AbSeqRep} objects)}
#'   \item{1 - generates plots for csv files}
#'   \item{2 - generates a report that collates all plots}
#'   \item{3 - generates interactive plots in report (default)}
#' }
#' each higher value also does what the previous values do. For example, \code{report = 2}
#' will return a named list of \linkS4class{AbSeqRep} objects, plot csv files,
#' and generate a (non-interactive)HTML report that collates all the plots together.
#' @param compare vector of strings. From the samples in \code{root}, samples can be
#' selected and compared against each other. For example, to compare "sample1" with "sample2" and
#' "sample3" with "sample4", \code{compare} should be c("sample1,sample2", "sample3,sample4").
#' @param BPPARAM BiocParallel backend. Configures the parallel implementation.
#' Refer to \href{https://bioconductor.org/packages/release/bioc/html/BiocParallel.html}{BiocParallel}
#' for more information.
#'
#' @return named list. List of \linkS4class{AbSeqRep} objects. The names of
#' the list are taken directly from the repertoire object itself. This return
#' value is consistent with the return value of \code{\link{report}}
#'
#' @seealso \linkS4class{AbSeqRep}
#'
#' @export
#'
#' @seealso \code{\link{report}}. Analogous function, but takes input from
#' an \linkS4class{AbSeqRep} or \linkS4class{AbSeqCRep} object instead.
#'
#' @examples
#' \dontrun{
#' # Assuming abseqPy has dumped its output in /path/to/output/directory/
#' # (i.e. the argument to --outdir / -o in abseqPy)
#'
#' ### report parameter usage example:
#'
#' # report = 0; don't plot, don't collate a HTML report, don't show anything interactive
#' samples <- abseqReport("/path/to/output/directory/", report = 0)
#' # samples is now a named list of AbSeqRep objects
#'
#' # report = 1; just plot pngs; don't collate a HTML report; nothing interactive
#' samples <- abseqReport("/path/to/output/directory/", report = 1)
#' # samples is now a named list of AbSeqRep objects
#'
#' # report = 2; plot pngs; collate a HTML report; HTML report will NOT be interactive
#' samples <- abseqReport("/path/to/output/directory/", report = 2)
#' # samples is now a named list of AbSeqRep objects
#'
#' # report = 3 (default); plot pngs; collate a HTML report; HTML report will be interactive
#' samples <- abseqReport("/path/to/output/directory/", report = 3)
#' # samples is now a named list of AbSeqRep objects
#'
#' #### Bonus Section (what to do with the list "samples"?):
#' # assume:
#' # > names(samples)
#' # [1] "Sample1" "Sampel2" "Sample3" "Sample4"
#'
#' # we want to explicitly compare Sample1 with Sample3:
#' new.combination <- samples[["Sample1"]] + samples[["Sample3"]]
#' # see abseqR::report for more information.
#' abseqR::report(new.combination)
#'
#' ### BPPARAM usage:
#'
#' # 4 core machine, use all cores -  use whatever value that suits you
#' nproc <- 4
#' samples <- abseqReport("/path/to/output/directory/",
#'                        BPPARAM = BiocParallel::MulticoreParam(nproc))
#'
#'
#' # run sequentially - no multiprocessing
#' samples <- abseqReport("/path/to/output/directory/",
#'                        BPPARAM = SerialParam())
#'
#'# see https://www.bioconductor.org/packages/devel/bioc/vignettes/BiocParallel/inst/doc/Introduction_To_BiocParallel.pdf
#'# for more information about how to use BPPARAM and BiocParallel in general.
#' }
abseqReport <- function(root, report, compare, BPPARAM) {
    #  ------ sanitize function arguments ---------
    root <- normalizePath(root)
    if (missing(report)) {
        report <- 3
    }
    if (missing(BPPARAM)) {
        BPPARAM <- BiocParallel::bpparam()
    }
    if (missing(compare)) {
        compare <- list.files(file.path(root, RESULT_DIR))
    } else {
        # make sure samples in "compare" actually exists
        availSamples <- list.files(file.path(root, RESULT_DIR))
        lapply(compare, function(s) {
            user.sample.name <- unlist(lapply(strsplit(s, ","), trimws))
            if (length(user.sample.name) > 1) {
                lapply(user.sample.name, function(si) {
                    if (!(si %in% availSamples)) {
                        stop(paste("Sample", si, "in 'compare' argument cannot",
                                   "be found in", file.path(root, RESULT_DIR)))
                    }
                })
            } else {
                if (!(s %in% availSamples)) {
                    stop(paste("Sample", s, "in 'compare' argument cannot",
                               "be found in", file.path(root, RESULT_DIR)))
                }
            }
        })

        compare <- unique(c(availSamples, compare))
    }

    # TODO: sanitize 'report' properly using constants/functions
    stopifnot(report %in% c(0, 1, 2, 3))
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
        #lapply(compare, function(pair) {
        bplapply(compare, function(pair) {
            sampleNames <- unlist(lapply(strsplit(pair, ","), trimws))

            # depending on the number of samples requested to plot, outputDir
            # is either a <sample>_vs_<sample> format or just <sample> meanwhile,
            # samples will either be a AbSeqCRep or just AbSeqRep.
            if (length(sampleNames) > 1) {
                outputDir <- file.path(root, RESULT_DIR, paste(sampleNames, collapse = "_vs_"))
                samples <- Reduce("+", lapply(sampleNames, function(sampleName) {
                    sample_ <- .loadAbSeqRepFromParams(file.path(root, RESULT_DIR, sampleName, ANALYSIS_PARAMS))
                    # sample@outdir should be the same as root
                    if (normalizePath(sample_@outdir) != root) {
                        warning(paste("Sample output directory", sample_@outdir,
                                      "is different from provided path", root,
                                      "assuming directory was moved"))
                        sample_@outdir <- root
                    }
                    return(sample_)
                }))
            } else {
                outputDir <- file.path(root, RESULT_DIR, sampleNames[1])
                samples <- .loadAbSeqRepFromParams(file.path(outputDir, ANALYSIS_PARAMS))
                if (normalizePath(samples@outdir) != root) {
                    warning(paste("Sample output directory", samples@outdir,
                                  "is different from provided path", root,
                                  "assuming directory was moved"))
                    samples@outdir <- root
                }
            }
            # due to the cache/shared user issue (see render() parallel github issue)
            # we delay the report generation until AFTER the multiprocessing part
            # has completed
            abseqR::report(samples,
                           outputDir,
                           report = 1)   # always plot, but DO NOT GENERATE REPORT!

        #})
        }, BPPARAM = BPPARAM)
    }

    # before we start creating reports, create the report directory
    if (report) {
        # create report directory
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
    }

    individualSamples <- list()
    individualReports <- list()
    # populate individualSample list with samples for user to browse and
    # create report if asked to.
    for (pair in compare) {
        sampleNames <- unlist(strsplit(pair, ","))
        if (length(sampleNames) == 1) {
            outputDir <- file.path(root, RESULT_DIR, sampleNames[1])
            samples <- .loadAbSeqRepFromParams(file.path(outputDir, ANALYSIS_PARAMS))
            if (normalizePath(samples@outdir) != root) {
                warning(paste("Sample output directory", samples@outdir,
                              "is different from provided path", root,
                              "assuming directory was moved"))
                samples@outdir <- root
            }
            individualSamples[[samples@name]] <- samples
        } else {
            outputDir <- file.path(root, RESULT_DIR, paste(sampleNames, collapse = "_vs_"))
            samples <- Reduce("+",
                              lapply(sampleNames, function(sampleName) {
                                  tmpsample <- .loadAbSeqRepFromParams(file.path(root, RESULT_DIR, sampleName, ANALYSIS_PARAMS))
                                  if (normalizePath(tmpsample@outdir) != root) {
                                      tmpsample@outdir <- root
                                  }
                                  return(tmpsample)
                              }))
        }
        if (report) {
            pth <- .generateReport(samples, root = outputDir,
                                   outputDir = file.path(root, ABSEQ_HTML_DIR, ABSEQ_NESTED_HTML_DIR),
                                   interactivePlot = interactivePlot,
                                   .indexHTML = file.path("..", "index.html"))
            if (!is.na(pth)) {
                individualReports[paste(sampleNames, collapse = "_vs_")] <- pth
            }
        }
    }

    if (length(individualReports)) {
        .collateReports(individualReports, individualSamples,
                        outputDirectory = file.path(root, ABSEQ_HTML_DIR))
    }
    return(individualSamples)
}


#' Collate all reports into a single directory and cretate an entry
#' \code{index.html} file that redirects to all other HTML files
#'
#' @import rmarkdown
#'
#' @param reports list/vector type. Collection of strings that are path(s)
#' to <sample>_report.html
#' @param individualSamples list type. list of AbSeqRep objects. Used to
#' extract filtering information and % read counts.
#' @param outputDirectory string type. Where should the report be placed.
.collateReports <- function(reports, individualSamples, outputDirectory) {
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
        .readSummary(file.path(x@outdir, RESULT_DIR, x@name), ABSEQ_RAW_READ_COUNT_KEY)
    }), collapse = ",")

    annotReads <- paste(lapply(individualSamples, function(x) {
        .readSummary(file.path(x@outdir, RESULT_DIR, x@name), ABSEQ_ANNOT_READ_COUNT_KEY)
    }), collapse = ",")

    filtReads <- paste(lapply(individualSamples, function(x) {
        .readSummary(file.path(x@outdir, RESULT_DIR, x@name), ABSEQ_FILT_READ_COUNT_KEY)
    }), collapse = ",")

    prodReads <- paste(lapply(individualSamples, function(x) {
        .readSummary(file.path(x@outdir, RESULT_DIR, x@name), ABSEQ_PROD_READ_COUNT_KEY)
    }), collapse = ",")

    filterSplitter <- "?"
    filters <- paste(lapply(individualSamples, function(x) {
        paste(paste0("Bitscore: ", paste(x@bitscore, collapse = " - ")),
              paste0("V alignment length: ", paste(x@alignlen, collapse = " - ")),
              paste0("Query start: ", paste(x@qstart, collapse = " - ")),
              paste0("Subject start: ", paste(x@sstart, collapse = " - ")),
              sep = ",")
    }), collapse = filterSplitter)

    templateFile <- system.file("extdata", "index.Rmd", package = "abseqR")
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
        analysisParams = file.path(individualSamples[[1]]@outdir, RESULT_DIR, individualSamples[[1]]@name, ANALYSIS_PARAMS)
    )
    if (rmarkdown::pandoc_available()) {
        rmarkdown::render(templateFile, output_dir = outputDirectory,
                          params = renderParams)
    } else {
        warning("Pandoc not found in system, will not create index.html")
    }
}
