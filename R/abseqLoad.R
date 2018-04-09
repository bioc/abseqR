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
              individualSamples <- list()

              root <- object@outputDirectory

              metaFile <- file.path(root, ABSEQ_CFG)

              con <- file(metaFile, "r")
              # reverse vector to plot single samples first, then multi-samples
              pairings <- rev(tail(readLines(con), n = -1))
              close(con)

              for (pair in pairings) {
                  sampleNames <- unlist(strsplit(pair, ","))

                  if (length(sampleNames) > 1) {
                      outputDir <- file.path(root, RESULT_DIR,
                                             paste(sampleNames,
                                                   collapse = "_vs_"))
                      sampleDirectories <- lapply(sampleNames, function(sampleName) {
                          file.path(root, RESULT_DIR, sampleName)
                      })
                      samples <- Reduce("+",
                                        lapply(sampleDirectories, function(directory) {
                                            .loadRepertoireFromParams(
                                                file.path(directory,
                                                          ANALYSIS_PARAMS))
                                        }))
                  } else {
                      outputDir <- file.path(root,
                                             RESULT_DIR,
                                             sampleNames[1])
                      samples <- .loadRepertoireFromParams(file.path(outputDir, ANALYSIS_PARAMS))
                      individualSamples <- c(individualSamples, samples)
                  }
                  AbSeq::plotRepertoires(samples, outputDir)
              }
              return(individualSamples)
          })
