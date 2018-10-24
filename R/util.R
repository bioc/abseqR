# single plot uses light blue colour
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

# AbSeq's output directory structure(s)
RESULT_DIR <- "auxiliary"
AUX_DIR <- "hdf"
# AbSeq's analysis directory names
ABSEQ_DIR_ANNOT <- "annot"
ABSEQ_DIR_PROD <- "productivity"
ABSEQ_DIR_ABUN <- "abundance"
ABSEQ_DIR_DIV <- "diversity"
ABSEQ_DIR_PAIR <- "clonotype_analysis"
ABSEQ_DIR_PRIM <- "primer_specificity"
ABSEQ_DIR_5UTR <- "utr5"
ABSEQ_DIR_SEC <- "secretion"
ABSEQ_HTML_DIR <- "report"
ABSEQ_NESTED_HTML_DIR <- "html_files"

# parameter file from AbSeq's run
ANALYSIS_PARAMS <- "analysis.params"

# AbSeq's summary file about the repertoire - raw/annot/prod counts
ABSEQ_SUMMARY <- "summary.txt"

# These are the "key"s from ABSEQ_SUMMARY file
# They should be in the format of: (for example)
# RawReads:<number>
# AnnotatedReads:<number>
ABSEQ_RAW_READ_COUNT_KEY <- "RawReads"
ABSEQ_ANNOT_READ_COUNT_KEY <- "AnnotatedReads"
ABSEQ_FILT_READ_COUNT_KEY <- "FilteredReads"
ABSEQ_PROD_READ_COUNT_KEY <- "ProductiveReads"


#' Checks if abseqPy has a metadata line that suggests
#' the orientation
#'
#' @param filename csv filename
#'
#' @return True if CSV metadata says "plot vertically"
.checkVert <- function(filename) {
    f <- file(filename, "r")
    res <- grepl("vert", readLines(f, n = 1), fixed = TRUE)
    close(f)
    return(res)
}

#' Get total number of samples (n)
#'
#' @description Often enough, the CSV values supplied do not contain
#' raw counts but percentages (so this value will let us know exactly
#' the sample size).
#'
#' @param filename csv filename
#'
#' @return string, sample size.
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
        retval <- list.files(
            path = p,
            pattern = pattern,
            full.names = TRUE,
            recursive = TRUE
        )
        # short circuit once at least one of the paths don't have the required
        # file. (Can't compare all the samples if one of them is missing!)
        if (length(retval) == 0) {
            return(c())
        }
        if (!(length(retval) %in% expectedRet)) {
            stop("Expected either ", paste(expectedRet, collapse = ", "),
                 " files to be found, but only ",
                 length(retval),
                 " were found.")
        }
        orderedFiles <- c(orderedFiles,  retval)
    }
    return(orderedFiles)
}

#' Returns all samples found under \code{sampleDirectory}
#'
#' @param sampleDirectory string, path to sample directory.
#'
#' @return un-normalized path to all samples under \code{sampleDirectory}
.inferAnalyzed <- function(sampleDirectory) {
    everything <- list.files(sampleDirectory)
    fullPath <- lapply(everything, function(x) {
        file.path(sampleDirectory, x)
    })
    return(everything[unlist(lapply(fullPath, dir.exists))])
}

#' Helper function to capitalize the first letter of \code{str}
#'
#' @param str string type
#'
#' @return string, \code{str} capitalized
.capitalize <- function(str) {
    firstLetter <- substr(str, 1, 1)
    rest <- substr(str, 2, nchar(str))
    return(paste0(toupper(firstLetter), rest))
}

#' Summary of dataframe
#'
#' @description Gives count, mean, standard deviation,
#' standard error of the mean, and confidence interval (default 95\%).
#'
#' adapted from http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper functions
#'
#' @param data a data frame.
#' @param measurevar the name of a column that contains the variable to be summariezed
#' @param groupvars a vector containing names of columns that contain grouping variables
#' @param na.rm a boolean that indicates whether to ignore NA's
#' @param conf.interval the percent range of the confidence interval (default is 95\%)
#' @param .drop logical.
#'
#' @import plyr stats
#'
#' @return dataframe
.summarySE <-
    function(data = NULL,
             measurevar,
             groupvars = NULL,
             na.rm = FALSE,
             conf.interval = .95,
             .drop = TRUE) {
        # New version of length which can handle NA's: if na.rm==T, don't count them
        length2 <- function(x, na.rm = FALSE) {
            if (na.rm) {
                sum(!is.na(x))
            } else {
                length(x)
            }
        }

        # This does the summary. For each group's data frame, return a vector with
        # N, mean, and sd
        datac <- ddply(
            data,
            groupvars,
            .drop = .drop,
            .fun = function(xx, col) {
                c(
                    N    = length2(xx[[col]], na.rm = na.rm),
                    mean = mean   (xx[[col]], na.rm = na.rm),
                    sd   = sd     (xx[[col]], na.rm = na.rm)
                )
            },
            measurevar
        )

        # Rename the "mean" column
        datac <- rename(datac, c("mean" = measurevar))

        datac$se <-
            datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

        # Confidence interval multiplier for standard error
        # Calculate t-statistic for confidence interval:
        # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
        ciMult <- qt(conf.interval / 2 + .5, datac$N - 1)
        datac$ci <- datac$se * ciMult

        return(datac)
    }

#' Helper function to return line types by importance based on provided
#' CD/Fs regions
#'
#' @description In the aesthetics of diversity plots (rarefaction, recapture,
#' and duplication), the line types should emphasise the most important
#' antibody region, they're ranked in ascending order of:
#' "FR4", "FR1", "FR2", "FR3", "CDR1", "CDR2", "CDR3", "V".
#'
#' @param regions a list/vector of strings (regions)
#'
#' @return vector of strings, each corresponding to the appropriate line type
#' for \code{regions}.
.getLineTypes <- function(regions) {
    if (length(regions) > 6) {
        stop("No line types for regions with length > 6 ")
    }
    regions <- unlist(lapply(regions, toupper))
    # order of importance: min -> max
    lvls <-
        c("FR4", "FR1", "FR2", "FR3", "CDR1", "CDR2", "CDR3", "V")
    # order of importance: max -> min
    lines <-
        c("solid",
          "twodash",
          "dotted",
          "dotdash",
          "longdash",
          "dashed")

    factorRegions <- factor(regions, levels = lvls)
    return(lines[order(factorRegions, decreasing = TRUE)])
}


#' Saves ggplot object as a Rdata file.
#'
#' @description It's a convinient function that does the check and saves
#' at the same time, for brevity within other areas of the code (to eliminate
#' repeated if checks).
#'
#' @import tools
#'
#' @param .save logical type. Whether or not we should save.
#' @param filename string.
#' @param plot ggplot object.
#'
#' @return nothing
.saveAs <- function(.save, filename, plot) {
    if (.save) {
        fname <- sub(tools::file_ext(filename), "Rdata", filename)
        save(file = fname, list = c("plot"))
    }
}


#' Return value specifed by key from AbSeq's summary file
#'
#' @param sampleRoot sample's root directory. For example,
#' \code{/path/to/<outputdir>/reports/<sample_name>}.
#' @param key character type. Possible values are (though they might change)
#' \itemize{
#'    \item{RawReads}
#'    \item{AnnotatedReads}
#'    \item{FilteredReads}
#'    \item{ProductiveReads}
#' }
#' @return value associated with key from summary file. "NA" (in string) if the field
#' is not available refer to util.R for the key values
.readSummary <- function(sampleRoot, key) {
    fname <- file.path(sampleRoot, ABSEQ_SUMMARY)
    con <- file(fname, "r")
    lines <- readLines(con)
    close(con)
    for (line in lines) {
        if (grepl(key, line, fixed = TRUE)) {
            return(strsplit(line, ":")[[1]][2])
        }
    }
    return("NA")
}


#' Substitutes the first occurance of `key` with `value` in `filename`
#'
#' @param filename character type
#' @param key character type
#' @param value character type
#' @param fixed logical type
#'
#' @return None
.substituteStringInFile <-
    function(filename, key, value, fixed = FALSE) {
        con <- file(filename, "r")
        lines <- readLines(con)
        close(con)
        lines <- sub(key, value, lines, fixed = fixed)
        con <- file(filename, "w")
        cat(lines, file = con, sep = "\n")
        close(con)
    }


#' Creates and returns an empty plot
#'
#' @import ggplot2
#'
#' @return empty ggplot2 object
.emptyPlot <- function() {
    # placeholder plot - prints nothing at all
    # https://www.r-bloggers.com/ggplot2-cheatsheet-for-visualizing-distributions/
    # https://github.com/mikessh/vdjtools/blob/master/src/main/resources/rscripts/intersect_pair_scatter.r
    g <- ggplot() + geom_point(aes(1, 1), colour = "white") +
        theme(
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            plot.margin = unit(c(3, -5.5, 4, 3), "mm")
        )
    return(g)
}


#' Given a dataframe with the columns "from", "to", and value.var, return
#' a symmetric matrix (with diagonal values = diag). I.e. a call to
#' isSymmetric(return_value_of_this_function) will always be TRUE.
#'
#'
#' @import reshape2
#'
#' @param dataframe dataframe with 3 required columns, namely:
#' +---------------------------------------+
#' | from | to | value.var | ...           |
#' +---------------------------------------+
#' |      |    |           |               |
#' +---------------------------------------+
#' where value.var is the string provided in the function parameter
#' @param value.var the column to use as the matrix value
#' @param diag what should the diagonal values be if the dataframe doesn't provide them
#' @param unidirectional logical type. If the dataframe provided has the reverse
#' pairs (i.e. a from-to pair AND a to-from pair with the save values in the
#' value.var column), then this should be FALSE. Otherwise, this function will
#' flip the from-to columns to generate a symmetric dataframe (and hence, a
#' symmetric matrix).
#'
#' @return a symmetric matrix with rownames(mat) == colnames(mat)
#' The diagonal values are filled with diag if the dataframe itself doesn't have
#' diagonal data
.loadMatrixFromDF <-
    function(dataframe,
             value.var,
             diag,
             unidirectional = TRUE) {
        if (unidirectional) {
            # swap the columns "from" and "to", while the others remain the same
            df.r <-
                dataframe[, c(2, 1, tail(seq_along(names(dataframe)), -2))]
            # rename the columns (after swapping, it's to - from, need it to be from - to)
            names(df.r) <- names(dataframe)
            # rowbind the dataframes into one
            df.f <- rbind(dataframe, df.r)
        } else {
            # bidirectional dataframe doesn't require mirror-ing
            df.f <- dataframe
        }
        mat <-
            reshape2::acast(df.f, from ~ to, value.var = value.var, fill = diag)
        # make sure the matrix is symmetric
        mat <- mat[, rownames(mat)]
        stopifnot(isSymmetric(mat))
        mat
    }


#' Given a directory = <abseqPy_outputdir>/RESULT_DIR/, returns the directories (repositories) in
#' 'directory'. That is, will not return any sample_vs_sample directories.
#' This is done by asserting that a 'repository'
#' must have an (analysis.params) file, and a summary.txt file.
#'
#' A sample_vs_sample directory will not have these files.
#'
#' @param directory string. Path up until <abseqPy_outputdir>/RESULT_DIR/
#'
#' @return vector of strings that are samples in 'directory', note, this is NOT
#' a full path, but just the sample/repertoire name itself
.findRepertoires <- function(directory) {
    repos <- list.files(directory, full.names = TRUE)
    # given a directory (d), return True if d is a repository
    .isRepo <- function(d) {
        all(c("analysis.params", "summary.txt") %in% list.files(d))
    }
    vapply(Filter(.isRepo, repos), basename, USE.NAMES = FALSE, FUN.VALUE = "")
}
