# single plot uses blue
BLUEHEX <- "#56B4E9"

# plot sizes
V_WIDTH <- 8
V_HEIGHT <- 5
H_WIDTH <- 5
H_HEIGHT <- 8
VENN_WIDTH <- 8
VENN_HEIGHT <- 8
V_WIDTH_L <- 12
V_HEIGHT_L <- 7.5

RESULT_DIR <- "report"
AUX_DIR <- "auxiliary"
ANALYSIS_PARAMS <- "analysis.params"


.checkVert <- function(filename) {
    f <- file(filename, "r")
    res <- grepl("vert", readLines(f, n = 1), fixed = TRUE)
    close(f)
    return(res)
}

.getTotal <- function(filename) {
    f <- file(filename, "r")
    res <- unlist(strsplit(readLines(f, n = 1), "="))[2]
    close(f)
    return(res)
}

.listFilesInOrder <- function(path, pattern, expectedRet = c(1)) {
    # Returns files in order of path.
    # i.e overrides list.files' behaviour of sorted files
    # This is crucial because we are assuming sampleNames is
    # in 1-1 correspondance with dataframes.
    # This can be achieved by manually iterating over all sample path
    # in the provided vector of path and appending to a vector.
    # Sometimes, (like abundance) list.files will return more that one
    # matching file (gene, family, variant), so expectedRet is there to ensure
    # that we aren't doing something silly. expectRet is a vector because
    # sometimes there might be more than one possible configuration.
    # EG: abundance plot may or maynot have D gene when analyzing heavy/light
    # chains
    # Returns: ordered vector of files (according to provided path's ordering)
    orderedFiles <- c()
    for (p in path) {
        retval <- list.files(path = p, pattern = pattern,
                             full.names = TRUE, recursive = TRUE)
        if (length(retval) == 0) {
            return(c())
        }
        if (!(length(retval) %in% expectedRet)) {
            stop(paste("Expected either", paste(expectedRet, collapse = ", ")),
                 "files to be found, but only", length(retval), "were found.")
        }
        orderedFiles <- c(orderedFiles,  retval)
    }
    return(orderedFiles)
}

.inferAnalyzed <- function(sampleDirectory) {
    everything <- list.files(sampleDirectory)
    fullPath <- lapply(everything, function(x) {
        file.path(sampleDirectory, x)
    })
    return(everything[unlist(lapply(fullPath, dir.exists))])
}

.capitalize <- function(str) {
    firstLetter <- substr(str, 1, 1)
    rest <- substr(str, 2, nchar(str))
    return(paste0(toupper(firstLetter), rest))
}

#' Summary of dataframe
#'
#' @description Gives count, mean, standard deviation,
#' standard error of the mean, and confidence interval (default 95%).
#'
#' adapted from http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper functions
#'
#' @import plyr
#'
#' @param data a data frame.
#' @param measurevar the name of a column that contains the variable to be summariezed
#' @param groupvars a vector containing names of columns that contain grouping variables
#' @param na.rm a boolean that indicates whether to ignore NA's
#' @param conf.interval the percent range of the confidence interval (default is 95%)
#' @param .drop NA
#'
#' @return dataframe
.summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm = FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                       c(N    = length2(xx[[col]], na.rm=na.rm),
                         mean = mean   (xx[[col]], na.rm=na.rm),
                         sd   = sd     (xx[[col]], na.rm=na.rm)
                       )
                   },
                   measurevar
    )

    # Rename the "mean" column
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

.getLineTypes <- function(regions) {
    if (length(regions) > 6) {
        stop("No line types for regions with length > 6 ")
    }
    regions <- unlist(lapply(regions, toupper))
    # order of importance: min -> max
    lvls <- c("FR4", "FR1", "FR2", "FR3", "CDR1", "CDR2", "CDR3", "V")
    # order of importance: max -> min
    lines <- c("solid", "twodash", "dotted", "dotdash", "longdash", "dashed")

    factorRegions <- factor(regions, levels = lvls)
    return(lines[order(factorRegions, decreasing = T)])
}
