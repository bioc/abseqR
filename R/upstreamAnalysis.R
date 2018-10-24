#' Plot upstream distribution
#'
#' @include util.R
#' @include distributions.R
#'
#' @import ggplot2
#'
#' @param upstreamDirectories list type. List of sample directories
#' @param upstreamOut string type. Output directory
#' @param expectedLength int type. Expected length of upstream sequences.
#' (i.e. upstream_end - upstream_start + 1).
#' @param upstreamLengthRange string type. start_end format
#' @param sampleNames vector type. 1-1 with upstream directories
#' @param combinedNames string type. Title friendly "combined" sample names
#' @param mashedNames string type. File friendly "mashed-up" sample names
#' @param .save logical type. Save Rdata?
#'
#' @return None
.plotUpstreamLength <-
    function(upstreamDirectories, upstreamOut, expectedLength,
             upstreamLengthRange, sampleNames,
             combinedNames, mashedNames, .save = TRUE) {
        message("Plotting upstream distributions for ", combinedNames)

        # full lengthed upstream sequences and upstream seqs
        # that are shorter than the expected length
        if (!is.infinite(expectedLength)) {
            lengths <- c("", "_short")
        } else {
            # Expected length is Inf
            lengths <- c("")
        }

        lapply(lengths, function(lengthType) {
            .plotUpstreamLengthDist(upstreamDirectories, upstreamOut,
                                    upstreamLengthRange, lengthType,
                                    sampleNames, combinedNames,
                                    mashedNames, .save = .save)

            .plotIGVUpstreamLenDist(upstreamDirectories, upstreamOut,
                                    upstreamLengthRange, lengthType,
                                    sampleNames, combinedNames,
                                    mashedNames, .save = .save)

            .plotIGVUpstreamLenDistDetailed(upstreamDirectories, upstreamOut,
                                            upstreamLengthRange, lengthType,
                                            sampleNames, combinedNames,
                                            mashedNames, .save = .save)
        })
}


#' Plots the validity of upstream sequences
#'
#' @description Plots the distribution of valid, faulty, and missing start
#' codon in IGV germlines (repeated for gene and family levels).
#'
#' @include util.R
#' @include distributions.R
#' @import ggplot2
#'
#' @param upstreamDirectories list type. List of sample directories
#' @param upstreamOut string type. Output directory
#' @param expectedLength int type. Expected length of upstream sequences.
#' (i.e. upstream_end - upstream_start + 1). If this is infinite, no plots
#' will be generated.
#' @param upstreamLengthRange string type. start_end format
#' @param sampleNames vector type. 1-1 with upstream directories
#' @param combinedNames string type. Title friendly "combined" sample names
#' @param mashedNames string type. File friendly "mashed-up" sample names
#' @param .save logical type. Save Rdata?
#'
#' @return None
.analyzeUpstreamValidity <- function(upstreamDirectories, upstreamOut,
                                     expectedLength, upstreamLengthRange,
                                     sampleNames, combinedNames, mashedNames,
                                     .save = TRUE) {

    message("Starting upstream analysis on ", combinedNames)

    if (!is.infinite(expectedLength)) {
        lapply(c("gene", "family"), function(lvl) {
            lapply(c("valid", "faulty", "no_atg"), function(status) {
                requestedFiles <-
                    .listFilesInOrder(path = upstreamDirectories,
                                      pattern = paste0( ".*_(secsig|5utr)_",
                                                        upstreamLengthRange,
                                                        "_", status, "_", lvl,
                                                        '\\.csv(\\.gz)?$'))
                if (length(requestedFiles)) {
                    fname <-  if (grepl("secsig", requestedFiles[[1]])) {
                        "secsig"
                    } else {
                        stopifnot(grepl("5utr", requestedFiles[[1]]))
                        "5utr"
                    }
                    if (status == "valid") {
                        if (fname == "secsig") {
                            title <- "Valid Secretion Signals"
                        } else {
                            title <- "Valid 5'-UTRs"
                        }
                    } else if (status == 'faulty') {
                        title <- "Faulty Translations"
                    } else {
                        title <- "Upstream sequences without start codon"
                    }

                    subtitle <- paste("Total is ",
                                      paste(lapply(requestedFiles, function(x) {
                                          as.integer(.getTotal(x))
                                      }), collapse = ", "))

                    plotVert <- .checkVert(requestedFiles[[1]])
                    g <- .plotDist(lapply(requestedFiles, read.csv, skip = 1),
                                   sampleNames, title, vert = plotVert,
                                   subs = subtitle)
                    if (plotVert) {
                        width <- V_WIDTH
                        height <- V_HEIGHT
                    } else {
                        width <- H_WIDTH
                        height <- H_WIDTH
                    }
                    saveName <-
                        file.path(upstreamOut,
                                  paste0(mashedNames, "_", fname, "_",
                                         upstreamLengthRange, "_", status, "_",
                                         lvl, ".png"))

                    ggsave(saveName, plot = g, width = width, height = height)
                    .saveAs(.save, saveName, plot = g)
                } else {
                    warning("Could not find ",
                            status,
                            " distribution for ",
                            lvl,
                            " for ",
                            combinedNames)
                }
            })
        })
    }
}


#' 5' UTR analysis
#'
#' @description Generates all the required plots for 5' UTR analysis. This
#' includes upstream length distributions and upstream sequence validity.
#'
#' @param utr5Directories list type. 5UTR directories where files are located
#' @param utr5Out string type. Where to dump output
#' @param sampleNames vector type. 1-1 with utr5Directories
#' @param combinedNames string type. Title friendly string
#' @param mashedNames string type. File name friendly string
#' @param upstreamRanges list type. Upstream ranges for each sample.
#' If length(utr5Directories) > 1, the plots will only be generated for
#' upstream ranges that are present in ALL samples. (i.e the intersection)
#' @param .save logical type, save Rdata?
#'
#' @return none
.UTR5Analysis <- function(utr5Directories,
                          utr5Out,
                          sampleNames,
                          combinedNames,
                          mashedNames,
                          upstreamRanges,
                          .save = TRUE) {
    message("Starting 5'UTR analysis on samples ", combinedNames)
    upstreamRange <- unique(upstreamRanges)
    if (length(upstreamRange) == 1 && is.numeric(upstreamRange[[1]])) {
        upRange <- upstreamRange[[1]]
        START <- 1
        END <- 2

        if (length(upRange) != 2) {
            stop("Expected range to only have start and stop values, but got ",
                 upRange, " instead.")
        }

        expectedLength <-
            as.numeric(upRange[END]) - as.numeric(upRange[START]) + 1

        if (is.infinite(expectedLength)) {
            upRangeString <- paste0(upRange[START], "_", "inf")
        } else {
            upRangeString <- paste0(expectedLength, "_", expectedLength)
        }

        .plotUpstreamLength(utr5Directories, utr5Out, expectedLength,
                            paste0(upRange[START], "_", upRange[END]),
                            sampleNames, combinedNames,
                            mashedNames, .save = .save)

        .analyzeUpstreamValidity(utr5Directories, utr5Out, expectedLength,
                                 upRangeString, sampleNames, combinedNames,
                                 mashedNames, .save = .save)

    } else if (length(upstreamRange) > 1) {
        warning("Found multiple different upstream ranges for samples ",
                combinedNames,
                ": ",
                upstreamRange,
                ". Will not plot comparisons.")
    } else {
        warning("UpstreamRange is not numeric for sample ", combinedNames)
    }
}



#' Secretion signal analysis
#'
#' @description Generates all the required plots for Secretion signal analysis.
#' This includes upstream length distributions and upstream sequence validity.
#'
#' @param secDirectories list type. Secretion signal directories where files
#' are located
#' @param secOut string type. Where to dump output
#' @param sampleNames vector type. 1-1 with secDirectories
#' @param combinedNames string type. Title friendly string
#' @param mashedNames string type. File name friendly string
#' @param upstreamRanges list type. Upstream ranges for each sample.
#' If length(secDirectories) > 1, the plots will only be generated for
#' upstream ranges that are present in ALL samples. (i.e. the intersection)
#' @param .save logical type, save Rdata?
#'
#' @return none
.secretionSignalAnalysis <-
    function(secDirectories,
             secOut,
             sampleNames,
             combinedNames,
             mashedNames,
             upstreamRanges,
             .save = TRUE) {
        message("Starting secretion signal analysis on samples ", combinedNames)
        upstreamRange <- unique(upstreamRanges)
        if (length(upstreamRange) == 1 && is.numeric(upstreamRange[[1]])) {

            upRange <- upstreamRange[[1]]
            START <- 1
            END <- 2

            if (length(upRange) != 2) {
                stop("Expected range to only have start and stop values,",
                     " but got ", upRange, " instead.")
            }

            expectedLength <-
                as.numeric(upRange[END]) - as.numeric(upRange[START]) + 1

            if (is.infinite(expectedLength)) {
                upRangeString <- paste0(upRange[START], "_", "inf")
            } else {
                upRangeString <- paste0(expectedLength, "_", expectedLength)
                upRangeTrimmedString <- paste0("1_", expectedLength - 1)
            }

            .plotUpstreamLength(secDirectories, secOut, expectedLength,
                                paste0(upRange[START], "_", upRange[END]),
                                sampleNames, combinedNames,
                                mashedNames, .save = .save)
            .analyzeUpstreamValidity(secDirectories, secOut, expectedLength,
                                     upRangeString, sampleNames, combinedNames,
                                     mashedNames, .save = .save)
            if (!is.infinite(expectedLength)) {
                .analyzeUpstreamValidity(secDirectories, secOut,
                                         expectedLength, upRangeTrimmedString,
                                         sampleNames, combinedNames,
                                         mashedNames, .save = .save)
            }

        } else if (length(upstreamRange) > 1) {
            warning("Found multiple different upstream ranges for samples ",
                    combinedNames,
                    ": ",
                    upstreamRange,
                    ". Will not plot comparisons.")
        } else {
            warning("UpstreamRange is not numeric for sample ", combinedNames)
        }
    }


#' Plot upstream sequence length distribution for upstream sequences (5'UTR or
#' secretion signal) for a given \code{upstreamLengthRange}
#'
#' @description Given an upstream length range, plot the
#' distribution of upstream sequence lengths.
#'
#' @param upstreamDirectories list type. List of sample directories
#' @param upstreamOut string type. Output directory
#' @param upstreamLengthRange The range of upstream sequences to be included in
#' this plot. This is usually determined by abseqPy and the format should be as
#' follows: "min_max", e.g.: 1_15 means range(1, 15) inclusive.string type.
#' @param lengthType string type. "" (the empty string) denotes everything
#' and "_short" denotes a short sequence. abseqPy dictates this because it's
#' used for locating the files.
#' @param sampleNames vector type. 1-1 with upstream directories
#' @param combinedNames string type. Title friendly "combined" sample names
#' @param mashedNames string type. File friendly "mashed-up" sample names
#' @param .save logical type. Save Rdata?
#'
#' @return None
.plotUpstreamLengthDist <- function(upstreamDirectories, upstreamOut,
                                    upstreamLengthRange, lengthType,
                                    sampleNames, combinedNames,
                                    mashedNames, .save) {
    seqLengthFiles <- .listFilesInOrder(path = upstreamDirectories,
                                        pattern = paste0(".*_(secsig|5utr)_",
                                                         upstreamLengthRange,
                                                         "_dist", lengthType,
                                                         "\\.csv(\\.gz)?$"))
    if (length(seqLengthFiles)) {
        # fname is either secsig or 5utr for now
        fname <- if (grepl("secsig", seqLengthFiles[[1]])) {
            "secsig"
        } else {
            stopifnot(grepl("5utr", seqLengthFiles[[1]]))
            "5utr"
        }

        g <- .plotSpectratype(lapply(seqLengthFiles, read.csv),
                              sampleNames, title = "Sequence lengths",
                              xlabel = "Sequence Length(bp)",
                              ylabel = "Distribution")
        saveName <- file.path(upstreamOut,
                              paste0(mashedNames, "_", fname,
                                     "_", upstreamLengthRange,
                                     "_dist", lengthType, ".png"))
        ggsave(saveName , plot = g, width = V_WIDTH, height = V_HEIGHT)
        .saveAs(.save, saveName, plot = g)
    } else {
        warning("Can't find upstream dist file for range ",
                 upstreamLengthRange,
                " for samples ",
                 combinedNames)
    }
}


#' Plot IGV family distribution for a given \code{upstreamLengthRange}
#'
#' @description Given an upstream length range,
#' plot the distributions of IGV family without showing their actual lengths. If
#' their actual lengths matter, refer
#' to \code{\link{.plotIGVUpstreamLenDistDetailed}}.
#'
#' @param upstreamDirectories list type. List of sample directories
#' @param upstreamOut string type. Output directory
#' @param upstreamLengthRange The range of upstream sequences to be included in
#' this plot. This is usually determined by abseqPy and the format should be as
#' follows: "min_max", e.g.: 1_15 means range(1, 15) inclusive.string type.
#' @param lengthType string type. "" (the empty string) denotes everything
#' and "_short" denotes a short sequence. abseqPy dictates this because it's
#' used for locating the files.
#' @param sampleNames vector type. 1-1 with upstream directories
#' @param combinedNames string type. Title friendly "combined" sample names
#' @param mashedNames string type. File friendly "mashed-up" sample names
#' @param .save logical type. Save Rdata?
#'
#' @return None
.plotIGVUpstreamLenDist <- function(upstreamDirectories, upstreamOut,
                                    upstreamLengthRange, lengthType,
                                    sampleNames, combinedNames,
                                    mashedNames, .save = TRUE) {

    seqClassLengthFiles <-
        .listFilesInOrder(path = upstreamDirectories,
                          pattern = paste0( ".*_(secsig|5utr)_",
                                            upstreamLengthRange,
                                            "_dist", lengthType, "_class",
                                            "\\.csv(\\.gz)?$"))
    if (length(seqClassLengthFiles)) {
        fname <- if (grepl("secsig", seqClassLengthFiles[[1]])) {
            "secsig"
        } else {
            stopifnot(grepl("5utr", seqClassLengthFiles[[1]]))
            "5utr"
        }

        subtitle <- paste("Total is ", paste(
            lapply(seqClassLengthFiles,
                   function(x) {
                       as.integer(.getTotal(x))
                   }),
            collapse = ", "
        ))

        plotVert <- .checkVert(seqClassLengthFiles[[1]])

        g <- .plotDist(lapply(seqClassLengthFiles, read.csv, skip = 1),
                       sampleNames,
                       paste("IGV Abundance in Sample", combinedNames),
                       plotVert, subs = subtitle)
        if (plotVert) {
            plotWidth <- V_WIDTH
            plotHeight <- V_HEIGHT
        } else {
            plotWidth <- H_WIDTH
            plotHeight <- H_HEIGHT
        }
        saveName <-
            file.path(upstreamOut, paste0(mashedNames, "_", fname, "_",
                                          upstreamLengthRange,
                                          "_dist", lengthType,
                                          "_class", ".png"))
        ggsave(saveName, plot = g,
               width = plotWidth, height = plotHeight)
        .saveAs(.save, saveName, plot = g)
    } else {
        warning("Can't find upstream dist file for samples ", combinedNames)
    }
}



#' Plots the detailed length distribution for IGV families
#'
#' @description A boxplot for each IGV families showing the IQR of upstream
#' lengths. In contrast to \code{\link{.plotIGVUpstreamLenDist}},
#' which only shows the distribution of IGV families
#' over \code{upstreamLengthRange}.
#'
#' @param upstreamDirectories list type. List of sample directories
#' @param upstreamOut string type. Output directory
#' @param upstreamLengthRange The range of upstream sequences to be included in
#' this plot. This is usually determined by abseqPy and the format should be as
#' follows: "min_max", e.g.: 1_15 means range(1, 15) inclusive.string type.
#' @param lengthType string type. "" (the empty string) denotes everything
#' and "_short" denotes a short sequence. abseqPy dictates this because it's
#' used for locating the files.
#' @param sampleNames vector type. 1-1 with upstream directories
#' @param combinedNames string type. Title friendly "combined" sample names
#' @param mashedNames string type. File friendly "mashed-up" sample names
#' @param .save logical type. Save Rdata?
#'
#' @return None
.plotIGVUpstreamLenDistDetailed <- function(upstreamDirectories,
                                            upstreamOut, upstreamLengthRange,
                                            lengthType, sampleNames,
                                            combinedNames, mashedNames,
                                            .save = TRUE) {

    # Plot upstream lengths for each IGV family ----------------------
    seqClassLengthBoxFiles <-
        .listFilesInOrder(path = upstreamDirectories,
                          pattern = paste0(".*_(5utr|secsig)_",
                                           upstreamLengthRange, "_dist",
                                           lengthType, "_class",
                                           "_box\\.csv(\\.gz)?$"))
    if (length(seqClassLengthBoxFiles)) {
        fname <- if (grepl("secsig", seqClassLengthBoxFiles[[1]])) {
            "secsig"
        } else {
            stopifnot(grepl("5utr", seqClassLengthBoxFiles[[1]]))
            "5utr"
        }

        g <- .boxPlot(lapply(seqClassLengthBoxFiles, read.csv),
                      sampleNames, paste("Sequence Lengths in",
                                         combinedNames))

        saveName <- file.path(upstreamOut,
                              paste0(mashedNames, "_", fname, "_",
                                     upstreamLengthRange, "_dist", lengthType,
                                     "_class", "_box.png"))

        ggsave(saveName, plot = g, width = V_WIDTH, height = V_HEIGHT)
        .saveAs(.save, saveName, plot = g)
    } else {
        warning("Can't find upstream dist files for samples ", combinedNames)
    }
}
