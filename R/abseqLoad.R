#' Title
#'
#' @slot outputDirectory character.
#'
#' @return
#' @export AbSeqLoad
#'
#' @examples
AbSeqLoad <- setClass("AbSeqLoad", slots = c(outputDirectory = "character"))


#' Title
#'
#' @import BiocParallel
#'
#' @include repertoire.R
#' @include compositeRepertoire.R
#' @include plotter.R
#' @include util.R
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples
setGeneric(name = "abSeqPlot",
           def = function(object) {
               standardGeneric("abSeqPlot")
           })

setMethod(f = "abSeqPlot",
          signature = "AbSeqLoad",
          definition = function(object) {

              root <- object@outputDirectory

              metaFile <- file.path(root, ABSEQ_CFG)

              con <- file(metaFile, "r")
              # reverse vector to plot single samples first, then multi-samples
              pairings <- rev(tail(readLines(con), n = -1))
              close(con)

              lapply(pairings, function(pair) {
              #BiocParallel::bplapply(pairings, function(pair) {
                  sampleNames <- unlist(strsplit(pair, ","))

                  if (length(sampleNames) > 1) {
                      outputDir <- file.path(root, RESULT_DIR,
                                             paste(sampleNames,
                                                   collapse = "_vs_"))
                      samples <- Reduce("+",
                                        lapply(sampleNames, function(sampleName) {
                                            .loadRepertoireFromParams(
                                                file.path(root, RESULT_DIR, sampleName, ANALYSIS_PARAMS))
                                        }))
                  } else {
                      outputDir <- file.path(root,
                                             RESULT_DIR,
                                             sampleNames[1])
                      samples <- .loadRepertoireFromParams(file.path(outputDir, ANALYSIS_PARAMS))
                  }
                  AbSeq::plotRepertoires(samples, outputDir)
              })

              individualSamples <- list()
              # populate individualSample list with samples for user to browse
              for (pair in pairings) {
                  sampleNames <- unlist(strsplit(pair, ","))
                  if (length(sampleNames) == 1) {
                      outputDir <- file.path(root,
                                             RESULT_DIR,
                                             sampleNames[1])
                      samples <- .loadRepertoireFromParams(file.path(outputDir, ANALYSIS_PARAMS))
                      individualSamples <- c(individualSamples)
                  }
              }
              return(individualSamples)
          })
