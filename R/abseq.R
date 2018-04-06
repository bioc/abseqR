
#' AbSeq S4 class
#'
#' @slot outputDirectory character type. Output directory of AbSeq analysis
#'
#' @return none
#' @export AbSeq
#'
#' @examples
AbSeq <- setClass("AbSeq",
                  slots = c(
                      outputDirectory = "character",
                      primer5File = "character",
                      primer3File = "character",
                      upstreamStart = "character",
                      upstreamEnd = "character"
                  ))


#' Todo
#'
#' @param object AbSeq object
#'
#' @include abundanceAnalysis.R
#' @include annotationAnalysis.R
#' @include distributions.R
#' @include diversityAnalysis.R
#' @include primerAnalysis.R
#' @include productivityAnalysis.R
#' @include upstreamAnalysis.R
#' @include util.R
#'
#' @import ggplot2
#'
#' @return None
#' @export
#'
#' @examples
#' todo
setGeneric(name = "abSeqPlot",
           def = function(object) {
               standardGeneric("abSeqPlot")
           })

setMethod(f = "abSeqPlot",
          signature(object = "AbSeq"),
          definition = function(object) {
              outputdir <- object@outputDirectory
              primer5File <- object@primer5File
              primer3File <- object@primer3File
              upstreamStart <- object@upstreamStart
              upstreamEnd <- object@upstreamEnd

              metaFile <- file.path(outputdir, ".rscripts_meta.tmp")

              pairings <- rev(scan(metaFile, character(), quote = "", skip = 1))

              # see what analysis has been done (list of strings, each a name of an
              # abseq analysis)
              analysis <- .inferAnalyzed((unlist(strsplit(pairings[1], "\\?")[1]))[1])

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

                  mashedNames <- paste(sampleNames, collapse = "_")
                  combinedNames <- paste(sampleNames, collapse = ", ")

                  # different logic in obtaining folder names and sample directory
                  # names depending on sample lengths
                  if (length(sampleNames) > 1) {
                      # /a/b/c/d/penultimate/sampleN_vs_sampleM
                      outputDir <- file.path(resultFolder,
                                             paste(sampleNames, collapse = "_vs_"))
                      dir.create(outputDir)
                  } else {
                      # /a/b/c/d/penultimate/sampleN
                      outputDir <- directories
                  }

                  # a little explaination: - the only variables you need to care about:

                  # directories - a vector of directories for this pairing
                  # (could be a singleton vector)
                  # (root directory for each sample within a pairing)
                  # eg: c("/a/b/c/PCR1_BHZ123-AGGGA-ACGTA_L001",
                  # "/a/b/cPCR2_BHZ123-ACGT_L001", "/a/b/cPCR3_BHZ123_AGCT-GGACT_L001")

                  # sampleNames - canonical names used by AbSeq
                  #               (a vector of them, sample length as directories)
                  # eg (for the above directories vector):
                  #       c("PCR1_L001", "PCR2_L001", "PCR3_L001")

                  # mashedNames - collapsing sampleNames with underscores for output png/pdfs

                  # Plotting begins now - a little more walkthrough.
                  # for each analysis, we:
                  #   1) create <analysis>Out as the output directory for the specific analysis
                  #   2) specialized directories to <analysis>Directories by pasting
                  #       "/<analysis>/" to each root sample directory
                  #   3) using the specialized directores, we use listFilesInOrder
                  #       (list.files returns alphabetically, but our samples
                  #       may be jumbled up when the
                  #       provided ordering isn't alphabetical order
                  #       (i.e. sampleNames isn't 1-1 with dataframes)) with
                  #       regex to easily (and safely) find the right samples for the right
                  #       function (within the analysis)
                  #   4) plot!

                  ##################################################
                  #                                                #
                  #               ANNOT PLOTS                      #
                  #                                                #
                  ##################################################
                  if ('annot' %in% analysis) {
                      annotOut <- file.path(outputDir, "annot")
                      if (!file.exists(annotOut)) {
                          dir.create(annotOut)
                      }
                      annotDirectories <- sapply(directories, file.path, "annot",
                                                 USE.NAMES = F)
                      .annotAnalysis(annotDirectories, annotOut, sampleNames, mashedNames)
                  }

                  ##################################################
                  #                                                #
                  #               ABUNDANCE PLOTS                  #
                  #                                                #
                  ##################################################
                  if ('abundance' %in% analysis) {
                      abunOut <- file.path(outputDir, "abundance")
                      if (!file.exists(abunOut)) {
                          dir.create(abunOut)
                      }
                      abundanceDirectories <- sapply(directories, file.path, "abundance",
                                                     USE.NAMES = F)
                      .abundanceAnalysis(abundanceDirectories, abunOut, sampleNames,
                                         combinedNames, mashedNames)
                  }

                  ##################################################
                  #                                                #
                  #              PRODUCTIVITY PLOTS                #
                  #                                                #
                  ##################################################
                  if ('productivity' %in% analysis) {
                      prodOut <- file.path(outputDir, "productivity")
                      if (!file.exists(prodOut)) {
                          dir.create(prodOut)
                      }
                      productivityDirectories <- sapply(directories, file.path,
                                                        "productivity", USE.NAMES = F)
                      .productivityAnalysis(productivityDirectories, prodOut,
                                            sampleNames, combinedNames, mashedNames)
                  }

                  ##################################################
                  #                                                #
                  #               DIVERSITY PLOTS                  #
                  #                                                #
                  ##################################################
                  if ('diversity' %in% analysis) {
                      diversityOut <- file.path(outputDir, "diversity")
                      if (!file.exists(diversityOut)) {
                          dir.create(diversityOut)
                      }
                      diversityDirectories <- sapply(directories, file.path, "diversity",
                                                     USE.NAMES = F)
                      .diversityAnalysis(diversityDirectories, diversityOut,
                                         sampleNames, mashedNames)
                  }

                  ##################################################
                  #                                                #
                  #               PRIMER.S  PLOTS                  #
                  #                                                #
                  ##################################################
                  if ('primer_specificity' %in% analysis) {
                      primerOut <- file.path(outputDir, "primer_specificity")
                      if (!file.exists(primerOut)) {
                          dir.create(primerOut)
                      }
                      primerDirectories <- sapply(directories, file.path,
                                                  "primer_specificity", USE.NAMES = F)
                      args <- commandArgs(trailingOnly = TRUE)
                      if (primer5File != "None" || primer3File != "None") {
                          .primerAnalysis(primer5File, primer3File, primerDirectories,
                                          primerOut, sampleNames,
                                          combinedNames, mashedNames)
                      }
                  }

                  ##################################################
                  #                                                #
                  #               UPSTREAM 5UTR PLOTS              #
                  #                                                #
                  ##################################################
                  if ('utr5' %in% analysis) {
                      utr5Out <- file.path(outputDir, "utr5")
                      if (!file.exists(utr5Out)) {
                          dir.create(utr5Out)
                      }
                      utr5Directories <- sapply(directories, file.path, "utr5",
                                                USE.NAMES = F)

                      if (upstreamStart != "None" && upstreamEnd != "None") {
                          expectedLength <- as.integer(upstreamEnd) - as.integer(upstreamStart) + 1
                          # NA => upstreamEnd is Inf
                          if (is.na(expectedLength)) {
                              upstreamRange <- paste0(upstreamStart, "_", 'inf')
                          } else {
                              upstreamRange <- paste0(expectedLength, "_", expectedLength)
                          }
                          .upstreamDist(utr5Directories, utr5Out, expectedLength,
                                        paste0(upstreamStart, "_", upstreamEnd),
                                        sampleNames, combinedNames, mashedNames, FALSE)
                          .upstreamAnalysis(utr5Directories, utr5Out, expectedLength,
                                            upstreamRange, sampleNames, combinedNames,
                                            mashedNames, FALSE)
                      }
                  }

                  ##################################################
                  #                                                #
                  #               UPSTREAM SEC.S PLOTS             #
                  #                                                #
                  ##################################################
                  if ("secretion" %in% analysis) {
                      secOut <- file.path(outputDir, "secretion")
                      if (!file.exists(secOut)) {
                          dir.create(secOut)
                      }
                      secDirectories <- sapply(directories, file.path, "secretion",
                                               USE.NAMES = F)
                      if (upstreamStart != "None" && upstreamEnd != "None") {
                          expectedLength <- as.integer(upstreamEnd) - as.integer(upstreamStart) + 1
                          # NA => upstream is Inf
                          if (is.na(expectedLength)) {
                              upstreamRange <- paste0(upstreamStart, "_", 'inf')
                          } else {
                              upstreamRange <- paste0(expectedLength, "_", expectedLength)
                              upstreamRangeTrimmed <- paste0("1_", expectedLength - 1)
                          }
                          .upstreamDist(secDirectories, secOut, expectedLength,
                                        paste0(upstreamStart, "_", upstreamEnd),
                                        sampleNames, combinedNames, mashedNames, TRUE)
                          .upstreamAnalysis(secDirectories, secOut, expectedLength,
                                            upstreamRange, sampleNames, combinedNames,
                                            mashedNames, TRUE)
                          if (!is.na(expectedLength)) {
                              .upstreamAnalysis(secDirectories, secOut, expectedLength,
                                                upstreamRangeTrimmed, sampleNames,
                                                combinedNames, mashedNames, TRUE)
                          }
                      }
                  }
              }
          })
