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
    # reverse vector to plot single samples first, then multi-samples
    pairings <- rev(tail(readLines(con), n = -1))
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
        AbSeq::plotRepertoires(samples, outputDir, report = FALSE)
    #})
    }, BPPARAM = BPPARAM)

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
            pth <- .generateReport(samples, outputDir, interactivePlot = interactivePlot)
            if (!is.na(pth)) {
                individualReports <- c(individualReports, pth)
            }
        }
    }

    if (length(individualReports)) {
        # todo - collate into index.html
    }

    return(individualSamples)
}
