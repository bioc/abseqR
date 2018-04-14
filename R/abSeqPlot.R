#' Plots all samples (and comparisons) as specified in abseq.cfg
#'
#' @include util.R
#' @include repertoire.R
#' @import BiocParallel
#'
#' @param root string type. Root directory (-o / --outdir in AbSeq)
#' @param BPPARAM BiocParallel backend
#'
#' @return list type. List of Repertoire objects
#' @export
#'
#' @examples todo
abSeqPlot <- function(root, BPPARAM = BiocParallel::bpparam()) {
    root <- normalizePath(root)
    metaFile <- file.path(root, ABSEQ_CFG)

    con <- file(metaFile, "r")
    # reverse vector to plot single samples first, then multi-samples
    pairings <- rev(tail(readLines(con), n = -1))
    close(con)

    lapply(pairings, function(pair) {
    #BiocParallel::bplapply(pairings, function(pair) {
        sampleNames <- unlist(strsplit(pair, ","))

        if (length(sampleNames) > 1) {
            outputDir <- file.path(root,
                                   RESULT_DIR,
                                   paste(sampleNames,
                                         collapse = "_vs_"))
            samples <- Reduce("+",
                              lapply(sampleNames, function(sampleName) {
                                  sample_ <- .loadRepertoireFromParams(file.path(root, RESULT_DIR, sampleName, ANALYSIS_PARAMS))
                                  # sample@outdir should be the same as
                                  if (normalizePath(sample_@outdir) != root) {
                                      message(
                                          paste0(
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
                message(paste0("Sample output directory is different from provided",
                        "path, assuming directory was moved"))
                samples@outdir <- root
            }
        }
        AbSeq::plotRepertoires(samples, outputDir)
    })
    #}, BPPARAM = BPPARAM)

    individualSamples <- list()
    # populate individualSample list with samples for user to browse
    for (pair in pairings) {
        sampleNames <- unlist(strsplit(pair, ","))
        if (length(sampleNames) == 1) {
            outputDir <- file.path(root,
                                   RESULT_DIR,
                                   sampleNames[1])
            samples <-
                .loadRepertoireFromParams(file.path(outputDir, ANALYSIS_PARAMS))
            individualSamples <-
                c(individualSamples, samples)
        }
    }
    return(individualSamples)
}
