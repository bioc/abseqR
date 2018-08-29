#' Collect primer names from FASTA
#'
#' @param primerFile string type. Path to primer file
#'
#' @return vector of primer names as seen in primerFile
.allPrimerNames <- function(primerFile) {
    primers <- c()
    fp = file(primerFile, open = "r")
    while (TRUE) {
        line = readLines(fp, n = 1)
        if (!length(line)) {
            break
        }
        if (startsWith(line, ">")) {
            primers <- c(primers, gsub("^>\\s*", "", line))
        }
    }
    close(fp)
    return(primers)
}

#' Convert file names to human friendly text
#'
#' @include util.R
#'
#' @param str string type
#'
#' @return string
.canonicalizeTitle <- function(str) {
    if (str == "outframe") {
        return("Out-of-frame")
    } else if (str == "indel_pos") {
        return("Abundance of Indel Positions in")
    } else if (str == "indelled") {
        return("Abundance of Indelled")
    } else {
        return(.capitalize(str))
    }
}

#' Collect the intersection of all primer names within a collection
#' of primer files
#'
#' @param primerFiles list / vector type. Collection of primer files
#'
#' @return vector type. Vector of primerNames that are present in ALL
#' primerFiles. NULL if there's no intersection at all
.commonPrimerNames <- function(primerFiles) {
    unionPrimers <- lapply(primerFiles, function(f) {
        if (file.exists(f)) {
            return(.allPrimerNames(f))
        } else {
            return(c())
        }
    })
    return(Reduce(intersect, unionPrimers))
}


#' Conducts primer specificity analysis
#'
#' @include util.R
#' @include distributions.R
#' @import ggplot2
#'
#' @param primerDirectories string type. Path to primer analysis directory
#' @param primer5Files vector / list type. 5' end primer files
#' @param primer3Files vector / list type. 3' end primer files
#' @param primerOut string type. output directory
#' @param sampleNames vector type. 1-1 with primerDirectories
#' @param combinedNames string type. Title friendly "combined" sample names
#' @param mashedNames string type. File friendly "mashed-up" sample names
#' @param .save logical type. Save Rdata?
#'
#' @return None
.primerAnalysis <- function(primerDirectories,
                            primer5Files,
                            primer3Files,
                            primerOut,
                            sampleNames,
                            combinedNames,
                            mashedNames,
                            .save = TRUE) {
    message(paste("Starting primer analysis for", combinedNames))

    primer5 <- .commonPrimerNames(primer5Files)
    primer3 <- .commonPrimerNames(primer3Files)

    allPrimers <- Filter(function(x) !is.null(x),
                         list("5" = primer5, "3" = primer3))

    category <- c("all", "productive", "outframe")
    primerIntegrityTypes <- c("stopcodon", "integrity", "indelled", "indel_pos")

    # all primer plots are repeated for productive and out-of-frame sequences
    lapply(category, function(ctg) {
        # each primer analysis is split to 5' and 3' ends (if exists)
        lapply(seq_along(allPrimers), function(i) {
            # get vector of primer names
            primerNames <- allPrimers[[i]]
            pend <- names(allPrimers)[[i]]
            message(paste("Primers found across all samples",
                          combinedNames,
                          "for",
                          paste0(pend, "'"),
                          "end are: ",
                          paste(primerNames, collapse = ", ")))

            lapply(X = primerNames, FUN = .plotPrimerIGVStatus, pend, ctg,
                   primerDirectories, sampleNames,
                   primerOut, combinedNames, mashedNames, .save)

            lapply(X = primerIntegrityTypes, FUN = .plotPrimerIntegrity, pend,
                   ctg, primerDirectories, sampleNames,
                   primerOut, combinedNames, mashedNames, .save)
        })
    })
}

#' Plots, for a given \code{category} and \code{pend}, the \code{primer}
#' IGV indelled distribution in a bar plot
#'
#' @description Plots the abundace of indelled primers relative to IGV germlines
#'
#' @include util.R
#' @include distributions.R
#' @import ggplot2
#'
#' @param primer string, primer name
#' @param pend string, either 3 or 5 (primer end)
#' @param category string, either "all", "productive", or "outframe"
#' @param primerDirectories string type. Path to primer analysis directory
#' @param sampleNames vector type. 1-1 with primerDirectories
#' @param primerOut string type. output directory
#' @param combinedNames string type. Title friendly "combined" sample names
#' @param mashedNames string type. File friendly "mashed-up" sample names
#' @param .save logical type. Save Rdata?
#'
#' @return None
.plotPrimerIGVStatus <- function(primer, pend, category,
                                 primerDirectories, sampleNames,
                                 primerOut, combinedNames, mashedNames,
                                 .save = TRUE) {
    files <-
        .listFilesInOrder(path = primerDirectories,
                          pattern = paste0( ".*_", category, "_",
                                      pend, "end_", primer,
                                      "_igv_dist\\.csv(\\.gz)?$"))
    # there is a slight chance that the user provided primer
    # has no hit, hence, not dist files
    # (AbSEq (obviously) doesn't plot something empty)
    if (length(files) == 0) {
        # nothing left to do
        return()
    }
    subtitle <- paste("Total is",
                      paste(lapply(files, function(x) {
                          as.integer(.getTotal(x))
                      }), collapse = ', '))
    vertical <- .checkVert(files[[1]])
    if (vertical) {
        plotWidth <- V_WIDTH
        plotHeight <- V_HEIGHT
    } else {
        plotWidth <- H_WIDTH
        plotHeight <- H_HEIGHT
    }
    primPlot <- .plotDist(lapply(files, read.csv, skip = 1),
                          sampleNames,
                          paste(paste0(
                              "IGV Abundance of indelled ",
                              primer, " (",
                              .canonicalizeTitle(category), ") in "),
                              combinedNames ),
                          vertical,
                          subs = subtitle)
    saveName <-
        file.path(primerOut, paste0(mashedNames, "_", category, "_",
                                    pend, "end_", primer, "_igv_dist.png"))
    ggsave(saveName, plot = primPlot, width = plotWidth, height = plotHeight)
    .saveAs(.save, saveName, plot = primPlot)
}

#' Plots the distribution of primer integrity for a given \code{category} and
#' 5' or 3' \code{pend}
#'
#' @include util.R
#' @include distributions.R
#' @import ggplot2
#'
#' @param primerIntegrity string. One of "stopcodon", "integrity",
#' "indelled", "indel_pos"
#' @param pend string, either 3 or 5 (primer end)
#' @param category string, either "all", "productive", or "outframe"
#' @param primerDirectories string type. Path to primer analysis directory
#' @param sampleNames vector type. 1-1 with primerDirectories
#' @param primerOut string type. output directory
#' @param combinedNames string type. Title friendly "combined" sample names
#' @param mashedNames string type. File friendly "mashed-up" sample names
#' @param .save logical type. Save Rdata?
#'
#' @return None
.plotPrimerIntegrity <- function(primerIntegrity, pend, category,
                                 primerDirectories, sampleNames, primerOut,
                                 combinedNames, mashedNames, .save = TRUE) {
    files <-
        .listFilesInOrder(path = primerDirectories,
                          pattern = paste0( ".*_", category, "_",
                                            pend, "end_",
                                            primerIntegrity,
                                            "_dist\\.csv(\\.gz)?$"))
    # there is a slight chance that the user provided
    # primer has no hit, hence, not dist files
    # (AbSEq (obviously) doesn't plot something empty)
    if (length(files) == 0) {
        # nothing left to do
        return()
    }

    subtitle <- paste("Total is",
                      paste(lapply(files, function(x) {
                          as.integer(.getTotal(x))
                      }), collapse = ', '))

    vertical <- .checkVert(files[[1]])
    if (vertical) {
        plotWidth <- V_WIDTH
        plotHeight <- V_HEIGHT
    } else {
        plotWidth <- H_WIDTH
        plotHeight <- H_HEIGHT
    }

    primPlot <- .plotDist(lapply(files, read.csv, skip = 1),
                          sampleNames,
                          paste(paste0(
                              .canonicalizeTitle(primerIntegrity),
                              " ", pend,
                              "'-end Primer Sequence (",
                              .canonicalizeTitle(category),
                              ") in" ),
                              combinedNames ),
                          vertical,
                          subs = subtitle)
    saveName <- file.path(primerOut,
                          paste0(mashedNames, "_", category, "_", pend,
                                 "end_", primerIntegrity, "_dist.png"))
    ggsave(saveName, plot = primPlot, width = plotWidth, height = plotHeight)
    .saveAs(.save, saveName, plot = primPlot)
}
