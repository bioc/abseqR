#' Plots all samples and comparisons as specified in \code{abseq.cfg}
#'
#' @include util.R
#' @include repertoire.R
#' @import BiocParallel
#'
#' @param root string type. Root directory as specified in
#' \code{-o} or \code{--outdir} in AbSeqPy.
#' @param report logical type. Generate HTML report(s) for sample(s)?
#' @param interactivePlot logical type. Use interactive plots (plotly) when generating
#' HTML report? (This argument is ignored if \code{report} is \code{FALSE})
#' @param BPPARAM BiocParallel backend. Configures the parallel implementation.
#' Refer to \href{https://bioconductor.org/packages/release/bioc/html/BiocParallel.html}{BiocParallel}
#' for more information.
#'
#' @return list type. List of \linkS4class{Repertoire} objects
#' @seealso \linkS4class{Repertoire}
#' @export
#'
#' @examples todo
abSeqPlot <- function(root, report = TRUE, interactivePlot = TRUE,
                      BPPARAM = BiocParallel::bpparam()) {
    if (!report && interactivePlot) {
        warning("report is set to FALSE, ignoring interactivePlots argument")
    }

    root <- normalizePath(root)
    metaFile <- file.path(root, ABSEQ_CFG)

    con <- file(metaFile, "r")
    pairings <- tail(readLines(con), n = -1)
    close(con)

    #lapply(pairings, function(pair) {
    BiocParallel::bplapply(pairings, function(pair) {
        sampleNames <- unlist(strsplit(pair, ","))

        # depending on the number of samples requested to plot, outputDir
        # is either a <sample>_vs_<sample> format or just <sample> meanwhile,
        # samples will either be a CompositeRepertoire or just Repertoire.
        if (length(sampleNames) > 1) {
            outputDir <- file.path(root, RESULT_DIR, paste(sampleNames, collapse = "_vs_"))
            samples <- Reduce("+",
                              lapply(sampleNames, function(sampleName) {
                                  sample_ <- .loadRepertoireFromParams(file.path(root, RESULT_DIR, sampleName, ANALYSIS_PARAMS))
                                  # sample@outdir should be the same as root
                                  if (normalizePath(sample_@outdir) != root) {
                                      message(
                                          paste(
                                              "Sample output directory is different from provided",
                                              "path, assuming directory was moved"
                                          )
                                      )
                                      sample_@outdir <- root
                                  }
                                  return(sample_)
                              }))
        } else {
            outputDir <- file.path(root, RESULT_DIR, sampleNames[1])
            samples <- .loadRepertoireFromParams(file.path(outputDir, ANALYSIS_PARAMS))
            if (normalizePath(samples@outdir) != root) {
                message(paste("Sample output directory is different from provided",
                        "path, assuming directory was moved"))
                samples@outdir <- root
            }
        }
        # due to the cache/shared user issue (see render() parallel github issue)
        # we delay the report generation until AFTER the multiprocessing part
        # has completed
        AbSeq::plotRepertoires(samples, outputDir, report = FALSE)
    #})
    }, BPPARAM = BPPARAM)

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
    for (pair in pairings) {
        sampleNames <- unlist(strsplit(pair, ","))
        if (length(sampleNames) == 1) {
            outputDir <- file.path(root, RESULT_DIR, sampleNames[1])
            samples <- .loadRepertoireFromParams(file.path(outputDir, ANALYSIS_PARAMS))
            if (normalizePath(samples@outdir) != root) {
                samples@outdir <- root
            }
            individualSamples <- c(individualSamples, samples)
        } else {
            outputDir <- file.path(root, RESULT_DIR, paste(sampleNames, collapse = "_vs_"))
            samples <- Reduce("+",
                              lapply(sampleNames, function(sampleName) {
                                  tmpsample <- .loadRepertoireFromParams(file.path(root, RESULT_DIR, sampleName, ANALYSIS_PARAMS))
                                  if (normalizePath(tmpsample@outdir) != root) {
                                      tmpsample@outdir <- root
                                  }
                                  return(tmpsample)
                              }))
        }
        if (report) {
            pth <- .generateReport(samples, root = outputDir,
                                   outputDir = file.path(root, ABSEQ_HTML_DIR, ABSEQ_NESTED_HTML_DIR),
                                   interactivePlot = interactivePlot)
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
#' @param individualSamples list type. list of Repertoire objects. Used to
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

    templateFile <- system.file("extdata", "index.Rmd", package = "AbSeq")
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
