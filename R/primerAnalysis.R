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

#' Convert short-form names to human friendly speech
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
#'
#' @return None
.primerAnalysis <- function(primerDirectories, primer5Files, primer3Files,
                           primerOut, sampleNames, combinedNames, mashedNames) {
    message(paste("Starting primer analysis for", combinedNames))

    primer5 <- .commonPrimerNames(primer5Files)
    primer3 <- .commonPrimerNames(primer3Files)

    allPrimers <- list(primer5, primer3)
    category <- c("all", "productive", "outframe")
    analysisType <- c("stopcodon", "integrity", "indelled", "indel_pos")


    for (ctg in category) {

        for (i in seq_along(allPrimers)) {
            # what end this is, 3? 5?
            primerNames <- allPrimers[[i]]
            if (!is.null(primerNames)) {
                pend <- if (i == 1) '5' else '3'

                message(paste("Primers found across all samples", combinedNames
                              ,"for", paste0(pend,"'"), "end are: ",
                              paste(primerNames, collapse = ", ")))

                # individual primer IGV abundance
                for (j in seq_along(primerNames)) {
                    files <-
                        .listFilesInOrder(path = primerDirectories,
                                          pattern = paste0(".*_", ctg,
                                                           "_", pend, "end_",
                                                           primerNames[j],
                                                           "_igv_dist\\.csv(\\.gz)?$"))

                    # there is a slight chance that the user provided primer
                    # has no hit, hence, not dist files
                    # (AbSEq (obviously) doesn't plot something empty)
                    if (length(files)) {
                        subtitle <- paste("Total is",
                                          paste(lapply(files, function(x) {
                                              as.integer(.getTotal(x))
                                          }), collapse = ', '))
                        vertical <- .checkVert(files[[1]])
                        primPlot <- .plotDist(
                            lapply(files, read.csv, skip = 1),
                            sampleNames,
                            paste(paste0("IGV Abundance of indelled ",
                                         primerNames[j]," (",
                                         .canonicalizeTitle(ctg), ") in "),
                                  combinedNames),
                            vertical,
                            subs = subtitle
                        )
                        if (vertical) {
                            ggsave(file.path(primerOut,
                                             paste0(mashedNames,"_", ctg,
                                                    "_", pend, "end_",
                                                    primerNames[j], "_igv_dist.png")),
                                   plot = primPlot, width = V_WIDTH,
                                   height = V_HEIGHT)
                        } else {
                            ggsave(file.path(primerOut,
                                             paste0(mashedNames,"_", ctg, "_",
                                                    pend, "end_",
                                                    primerNames[j], "_igv_dist.png")),
                                   plot = primPlot,
                                   width = H_WIDTH, height = H_HEIGHT)
                        }
                    }  # endof if length(files)
                }   # endof individual primer IGV abundance


                # stopcodon, integrity, indelled and indel_pos analysis
                # on all categories
                for (j in seq_along(analysisType)) {
                    files <-
                        .listFilesInOrder(path = primerDirectories,
                                          pattern = paste0(".*_", ctg, "_",
                                                           pend, "end_",
                                                           analysisType[j],
                                                           "_dist\\.csv(\\.gz)?$"))

                    # there is a slight chance that the user provided
                    # primer has no hit, hence, not dist files
                    # (AbSEq (obviously) doesn't plot something empty)
                    if (length(files)) {
                        subtitle <- paste("Total is",
                                          paste(lapply(files, function(x) {
                                              as.integer(.getTotal(x))
                                          }), collapse = ', '))

                        vertical <- .checkVert(files[[1]])
                        primPlot <- .plotDist(
                            lapply(files, read.csv, skip = 1),
                            sampleNames,
                            paste(paste0(.canonicalizeTitle(analysisType[j]),
                                         " ", pend,
                                         "'-end Primer Sequence (",
                                         .canonicalizeTitle(ctg), ") in"),
                                  combinedNames),
                            vertical,
                            subs = subtitle
                        )
                        if (vertical) {
                            ggsave(file.path(primerOut,
                                             paste0(mashedNames,"_", ctg,
                                                    "_", pend, "end_",
                                                    analysisType[j], "_dist.png")),
                                   plot = primPlot, width = V_WIDTH,
                                   height = V_HEIGHT)
                        } else {
                            ggsave(file.path(primerOut,
                                             paste0(mashedNames,"_", ctg, "_",
                                                    pend, "end_", analysisType[j],
                                                    "_dist.png")),
                                   plot = primPlot,
                                   width = H_WIDTH, height = H_HEIGHT)
                        }
                    }  # endof if length(files)
                } # endof "type" analysis
            } # endof if !is.null(primerNames)
        } # end of for (i in allPrimers)
    } # end of category
}
