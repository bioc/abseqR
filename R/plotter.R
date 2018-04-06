#'
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
#'
#' @param sampleNames
#' @param directories
#' @param analysis
#' @param outputDir
#' @param primer5File
#' @param primer3File
#' @param upstreamStart
#' @param upstreamEnd
#'
#' @return None
.plotSamples <-
    function(sampleNames,
             directories,
             analysis,
             outputDir,
             primer5File,
             primer3File,
             upstreamStart,
             upstreamEnd) {
        mashedNames <- paste(sampleNames, collapse = "_")
        combinedNames <- paste(sampleNames, collapse = ", ")
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
            abundanceDirectories <-
                sapply(directories, file.path, "abundance",
                       USE.NAMES = F)
            .abundanceAnalysis(abundanceDirectories,
                               abunOut,
                               sampleNames,
                               combinedNames,
                               mashedNames)
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
            .productivityAnalysis(productivityDirectories,
                                  prodOut,
                                  sampleNames,
                                  combinedNames,
                                  mashedNames)
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
            diversityDirectories <-
                sapply(directories, file.path, "diversity",
                       USE.NAMES = F)
            .diversityAnalysis(diversityDirectories,
                               diversityOut,
                               sampleNames,
                               mashedNames)
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
            primerDirectories <- sapply(directories,
                                        file.path,
                                        "primer_specificity",
                                        USE.NAMES = F)
            if (primer5File != "None" || primer3File != "None") {
                .primerAnalysis(
                    primer5File,
                    primer3File,
                    primerDirectories,
                    primerOut,
                    sampleNames,
                    combinedNames,
                    mashedNames
                )
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
                expectedLength <-
                    as.integer(upstreamEnd) - as.integer(upstreamStart) + 1
                # NA => upstreamEnd is Inf
                if (is.na(expectedLength)) {
                    upstreamRange <- paste0(upstreamStart, "_", 'inf')
                } else {
                    upstreamRange <- paste0(expectedLength, "_", expectedLength)
                }
                .upstreamDist(
                    utr5Directories,
                    utr5Out,
                    expectedLength,
                    paste0(upstreamStart, "_", upstreamEnd),
                    sampleNames,
                    combinedNames,
                    mashedNames,
                    FALSE
                )
                .upstreamAnalysis(
                    utr5Directories,
                    utr5Out,
                    expectedLength,
                    upstreamRange,
                    sampleNames,
                    combinedNames,
                    mashedNames,
                    FALSE
                )
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
            secDirectories <-
                sapply(directories, file.path, "secretion",
                       USE.NAMES = F)
            if (upstreamStart != "None" && upstreamEnd != "None") {
                expectedLength <-
                    as.integer(upstreamEnd) - as.integer(upstreamStart) + 1
                # NA => upstream is Inf
                if (is.na(expectedLength)) {
                    upstreamRange <- paste0(upstreamStart, "_", 'inf')
                } else {
                    upstreamRange <- paste0(expectedLength, "_", expectedLength)
                    upstreamRangeTrimmed <-
                        paste0("1_", expectedLength - 1)
                }
                .upstreamDist(
                    secDirectories,
                    secOut,
                    expectedLength,
                    paste0(upstreamStart, "_", upstreamEnd),
                    sampleNames,
                    combinedNames,
                    mashedNames,
                    TRUE
                )
                .upstreamAnalysis(
                    secDirectories,
                    secOut,
                    expectedLength,
                    upstreamRange,
                    sampleNames,
                    combinedNames,
                    mashedNames,
                    TRUE
                )
                if (!is.na(expectedLength)) {
                    .upstreamAnalysis(
                        secDirectories,
                        secOut,
                        expectedLength,
                        upstreamRangeTrimmed,
                        sampleNames,
                        combinedNames,
                        mashedNames,
                        TRUE
                    )
                }
            }
        }
    }
