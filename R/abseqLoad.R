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

              outputdir <- object@outputDirectory

              metaFile <- file.path(outputdir, ".rscripts_meta.tmp")

              pairings <- rev(scan(metaFile, character(), quote = "", skip = 1))

              for (i in seq_along(pairings)) {
                  # get pairings, split by "?", you'll have something like
                  # c("dirname1,dirname2,...", "sampleName1, sampleName2, ...")
                  # i.e. pair[1] is a comma separated string of directories
                  # i.e. pair[2] is a comma separated string of sample names
                  pair <- unlist(strsplit(pairings[i], "\\?"))
                  directories <- unlist(strsplit(pair[1], ","))
                  sampleNames <- unlist(strsplit(pair[2], ","))

                  if (length(directories) != length(sampleNames)) {
                      stop(paste("Expected length of directories to be the same as sampleNames",
                           "but got", length(directories), "and", length(sampleNames),
                           "instead"))
                  }

                  # to get the result folder, we just need to look at any one of them,
                  # and take the full path until the penultimate directory
                  # /a/b/c/d/penultimate/<sample_1_dir>
                  # NOTE: using normalizePath to make sure there's no trailing separator
                  decomposed <- unlist(strsplit(normalizePath(directories[1]),
                                         .Platform$file.sep))
                  # /a/b/c/d/penultimate
                  resultFolder <- do.call("file.path",
                                          as.list(head(decomposed, n =
                                                           length(decomposed) - 1)))

                  # different logic in obtaining folder names and sample directory
                  # names depending on sample lengths
                  if (length(sampleNames) > 1) {
                      # /a/b/c/d/penultimate/sampleN_vs_sampleM
                      outputDir <- file.path(resultFolder, paste(sampleNames, collapse = "_vs_"))
                      samples <- Reduce("+",
                                        lapply(directories, function(directory) {
                                            .loadRepertoireFromParams(
                                                file.path(directory,
                                                          ANALYSIS_PARAMS))
                                        }))
                  } else {
                      # /a/b/c/d/penultimate/sampleN
                      outputDir <- directories
                      samples <- .loadRepertoireFromParams(file.path(directories[[1]], ANALYSIS_PARAMS))
                      individualSamples <- c(individualSamples, samples)
                  }
                  AbSeq::plotRepertoires(samples, outputDir)
              }
              return(individualSamples)
          })
