#' Conducts upstream analysis
#'
#' @include util.R
#' @include distributions.R
#' @import ggplot2
#'
#' @param upstreamDirectories list type. List of sample directories
#' @param upstreamOut string type. Output directory
#' @param expectedLength int type. Expected length of upstream sequences.
#' Can be infinite
#' @param upstreamLengthRange string type. start_end format
#' @param sampleNames vector type. 1-1 with upstream directories
#' @param combinedNames string type. Title friendly "combined" sample names
#' @param mashedNames string type. File friendly "mashed-up" sample names
#' @param secsig boolean type. True if upstream is for secretion signal,
#' False if it's for 5UTR
#'
#' @return None
.upstreamAnalysis <- function(upstreamDirectories,
                              upstreamOut,
                              expectedLength,
                              upstreamLengthRange,
                              sampleNames,
                              combinedNames,
                              mashedNames,
                              secsig) {
    fname <-  if (secsig)
        "secsig"
    else
        "5utr"

    message(paste("Starting upstream analysis on", combinedNames))

    if (!is.infinite(expectedLength)) {
        level <- c("gene", "family")
        for (lvl in level) {
            statuses <- c("valid", "faulty", "no_atg")
            for (status in statuses) {
                requestedFiles <-
                    .listFilesInOrder(
                        path = upstreamDirectories,
                        pattern = paste0(
                            ".*_",
                            fname,
                            "_",
                            upstreamLengthRange,
                            "_",
                            status,
                            "_",
                            lvl,
                            '\\.csv(\\.gz)?$'
                        )
                    )
                if (length(requestedFiles) > 0) {
                    if (status == "valid") {
                        if (secsig) {
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
                                      paste(
                                          lapply(requestedFiles, function(x) {
                                              as.integer(.getTotal(x))
                                          }),
                                          collapse = ", "
                                      ))

                    plotVert <- .checkVert(requestedFiles[[1]])

                    g <-
                        .plotDist(
                            lapply(requestedFiles, read.csv, skip = 1),
                            sampleNames,
                            title,
                            vert = plotVert,
                            subs = subtitle
                        )
                    if (plotVert) {
                        width <- V_WIDTH

                        height <- V_HEIGHT

                    } else {
                        width <- H_WIDTH

                        height <- H_WIDTH

                    }
                    ggsave(
                        file.path(
                            upstreamOut,
                            paste0(
                                mashedNames,
                                "_",
                                fname,
                                "_",
                                upstreamLengthRange,
                                "_",
                                status,
                                "_",
                                lvl,
                                ".png"
                            )
                        ),
                        plot = g,
                        width = width,
                        height = height
                    )
                }
            }
        }
    }

}

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
#' Can be infinite
#' @param upstreamLengthRange string type. start_end format
#' @param sampleNames vector type. 1-1 with upstream directories
#' @param combinedNames string type. Title friendly "combined" sample names
#' @param mashedNames string type. File friendly "mashed-up" sample names
#' @param secsig boolean type. True if upstream is for secretion signal,
#' False if it's for 5UTR
#'
#' @return None
.upstreamDist <-
    function(upstreamDirectories,
             upstreamOut,
             expectedLength,
             upstreamLengthRange,
             sampleNames,
             combinedNames,
             mashedNames,
             secsig) {
        message(paste("Plotting upstream distributions for", combinedNames))

        fname <-  if (secsig)
            "secsig"
        else
            "5utr"
        # full lengthed upstream sequences and upstream seqs
        # that are shorter than the expected length
        if (!is.infinite(expectedLength)) {
            lengths <- c("", "_short")
        } else {
            # Expected length is Inf
            lengths <- c("")
        }

        for (len in lengths) {
            seqLengthFiles <- .listFilesInOrder(
                path = upstreamDirectories,
                pattern = paste0(
                    ".*_",
                    fname,
                    "_",
                    upstreamLengthRange,
                    "_dist",
                    len,
                    "\\.csv(\\.gz)?$"
                )
            )
            if (length(seqLengthFiles) > 0) {
                g <- .plotSpectratype(
                    lapply(seqLengthFiles, read.csv),
                    sampleNames,
                    title = "Sequence lengths",
                    xlabel = "Sequence Length(bp)",
                    ylabel = "Distribution"
                )
                ggsave(
                    file.path(
                        upstreamOut,
                        paste0(
                            mashedNames,
                            "_",
                            fname,
                            "_",
                            upstreamLengthRange,
                            "_dist",
                            len,
                            ".png"
                        )
                    ),
                    plot = g,
                    width = V_WIDTH,
                    height = V_HEIGHT
                )
            }

            # class level
            seqClassLengthFiles <-
                .listFilesInOrder(
                    path = upstreamDirectories,
                    pattern = paste0(
                        ".*_",
                        fname,
                        "_",
                        upstreamLengthRange,
                        "_dist",
                        len,
                        "_class",
                        "\\.csv(\\.gz)?$"
                    )
                )
            if (length(seqClassLengthFiles) > 0) {
                subtitle <- paste("Total is ", paste(
                    lapply(seqClassLengthFiles,
                           function(x) {
                               as.integer(.getTotal(x))
                           }),
                    collapse = ", "
                ))

                plotVert <- .checkVert(seqClassLengthFiles[[1]])

                g <-
                    .plotDist(
                        lapply(seqClassLengthFiles, read.csv, skip = 1),
                        sampleNames,
                        paste("IGV Abundance in Sample", combinedNames),
                        plotVert,
                        subs = subtitle
                    )

                if (plotVert) {
                    ggsave(
                        file.path(
                            upstreamOut,
                            paste0(
                                mashedNames ,
                                "_",
                                fname,
                                "_",
                                upstreamLengthRange,
                                "_dist",
                                len,
                                "_class",
                                ".png"
                            )
                        ),
                        plot = g,
                        width = V_WIDTH,
                        height = V_HEIGHT
                    )
                } else {
                    ggsave(
                        file.path(
                            upstreamOut,
                            paste0(
                                mashedNames,
                                "_",
                                fname,
                                "_",
                                upstreamLengthRange,
                                "_dist",
                                len,
                                "_class",
                                ".png"
                            )
                        ),
                        plot = g,
                        width = H_WIDTH,
                        height = H_HEIGHT
                    )
                }
            }

            # box plot for class level
            seqClassLengthBoxFiles <-
                .listFilesInOrder(
                    path = upstreamDirectories,
                    pattern = paste0(
                        ".*_",
                        fname,
                        "_",
                        upstreamLengthRange,
                        "_dist",
                        len,
                        "_class",
                        "_box\\.csv(\\.gz)?$"
                    )
                )
            if (length(seqClassLengthBoxFiles) > 0) {
                g <- .boxPlot(
                    lapply(seqClassLengthBoxFiles, read.csv),
                    sampleNames,
                    paste("Sequence Lengths in", combinedNames)
                )
                ggsave(
                    file.path(
                        upstreamOut,
                        paste0(
                            mashedNames,
                            "_",
                            fname,
                            "_",
                            upstreamLengthRange,
                            "_dist",
                            len,
                            "_class",
                            "_box.png"
                        )
                    ),
                    plot = g,
                    width = V_WIDTH,
                    height = V_HEIGHT
                )
            }
        }
    }



#' Title
#'
#' @param utr5Directories list type. 5UTR directories where files are located
#' @param utr5Out string type. Where to dump output
#' @param sampleNames vector type. 1-1 with utr5Directories
#' @param combinedNames string type. Title friendly string
#' @param mashedNames string type. File name friendly string
#' @param upstreamRanges list type. Upstream ranges for each sample.
#' If length(utr5Directories) > 1, the plots will only be generated for
#' upstream ranges that are present in ALL samples. (i.e the intersection)
#'
#' @return none
.UTR5Analysis <- function(utr5Directories,
                          utr5Out,
                          sampleNames,
                          combinedNames,
                          mashedNames,
                          upstreamRanges) {
    message(paste("Starting 5'UTR analysis on samples", combinedNames))
    upstreamRange <- unique(upstreamRanges)
    if (length(upstreamRange) == 1 &&
        is.numeric(upstreamRange[[1]])) {
        upRange <- upstreamRange[[1]]
        START <- 1
        END <- 2
        if (length(upRange) != 2) {
            stop(
                paste(
                    "Expected range to only have start and stop values, but got",
                    upRange,
                    "instead"
                )
            )
        }
        expectedLength <-
            as.numeric(upRange[END]) - as.numeric(upRange[START]) + 1
        if (is.infinite(expectedLength)) {
            upRangeString <- paste0(upRange[START], "_", "inf")
        } else {
            upRangeString <- paste0(expectedLength, "_", expectedLength)
        }
        .upstreamDist(
            utr5Directories,
            utr5Out,
            expectedLength,
            paste0(upRange[START], "_", upRange[END]),
            sampleNames,
            combinedNames,
            mashedNames,
            FALSE
        )

        .upstreamAnalysis(
            utr5Directories,
            utr5Out,
            expectedLength,
            upRangeString,
            sampleNames,
            combinedNames,
            mashedNames,
            FALSE
        )
    } else if (length(upstreamRange) > 1) {
        warning(
            paste(
                "Found multiple different upstream ranges for samples",
                combinedNames,
                ":",
                upstreamRange,
                "Will not plot comparisons."
            )
        )
    } else {
        warning(paste("UpstreamRange is not numeric for sample", combinedNames))
    }
}



#' Title
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
#'
#' @return none
.secretionSignalAnalysis <-
    function(secDirectories,
             secOut,
             sampleNames,
             combinedNames,
             mashedNames,
             upstreamRanges) {
        message(paste(
            "Starting secretion signal analysis on samples",
            combinedNames
        ))
        upstreamRange <- unique(upstreamRanges)
        if (length(upstreamRange) == 1 &&
            is.numeric(upstreamRange[[1]])) {
            upRange <- upstreamRange[[1]]
            START <- 1
            END <- 2
            if (length(upRange) != 2) {
                stop(
                    paste(
                        "Expected range to only have start and stop values, but got",
                        upRange,
                        "instead"
                    )
                )
            }
            expectedLength <-
                as.numeric(upRange[END]) - as.numeric(upRange[START]) + 1
            if (is.infinite(expectedLength)) {
                upRangeString <- paste0(upRange[START], "_", "inf")
            } else {
                upRangeString <- paste0(expectedLength, "_", expectedLength)
                upRangeTrimmedString <- paste0("1_", expectedLength - 1)
            }
            .upstreamDist(
                secDirectories,
                secOut,
                expectedLength,
                paste0(upRange[START], "_", upRange[END]),
                sampleNames,
                combinedNames,
                mashedNames,
                TRUE
            )
            .upstreamAnalysis(
                secDirectories,
                secOut,
                expectedLength,
                upRangeString,
                sampleNames,
                combinedNames,
                mashedNames,
                TRUE
            )
            if (!is.infinite(expectedLength)) {
                .upstreamAnalysis(
                    secDirectories,
                    secOut,
                    expectedLength,
                    upRangeTrimmedString,
                    sampleNames,
                    combinedNames,
                    mashedNames,
                    TRUE
                )
            }
        } else if (length(upstreamRange) > 1) {
            warning(
                paste(
                    "Found multiple different upstream ranges for samples",
                    combinedNames,
                    ":",
                    upstreamRange,
                    "Will not plot comparisons."
                )
            )
        } else {
            warning(paste("UpstreamRange is not numeric for sample", combinedNames))
        }
    }
