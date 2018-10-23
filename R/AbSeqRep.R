#' S4 class - AbSeqRepertoire analysis object
#'
#' @import methods
#'
#' @description The AbSeqRep object contains all metadata associated with
#' the AbSeq (python backend) run conducted on it. This S4 class represents
#' a single sample(repertoire) and it can be "combined" with other samples
#' by using the \code{+} operator to create an \linkS4class{AbSeqCRep} object.
#' This value, in turn, can be used as the first argument to
#' \link{report} which generates a comparison between all samples included
#' in the \code{+} operation.
#'
#' Users do not manually construct this class, but rather indirectly
#' obtain this class object as a return value
#' from the \link{abseqReport} and \link{report} functions.
#'
#' @slot f1 character. Path to FASTA/FASTQ file 1.
#' @slot f2 character. Path to FASTA/FASTQ file 2.
#' @slot chain character. Type of chain, possible values:
#' \itemize{
#'   \item{hv}
#'   \item{lv}
#'   \item{kv}
#'   \item{klv}
#' }
#' each representing \bold{H}eavy, \bold{L}ambda and \bold{K}appa respectively.
#' @slot task character. Type of analysis conducted, possible values:
#' \itemize{
#'   \item{all}
#'   \item{annotate}
#'   \item{abundance}
#'   \item{diversity}
#'   \item{productivity}
#'   \item{fastqc}
#'   \item{primer}
#'   \item{5utr}
#'   \item{rsasimple}
#'   \item{seqlen}
#'   \item{secretion}
#'   \item{seqlenclass}
#' }
#' @slot name character. Name of analysis.
#' @slot bitscore numeric. Part of filtering criteria: V gene bitscore filter value.
#' @slot qstart numeric. Part of filtering criteria: V gene query start filter value.
#' @slot sstart numeric. Part of filtering criteria: V gene subject start filter value.
#' @slot alignlen numeric. Part of filtering criteria: V gene alignment length filter value.
#' @slot clonelimit numeric. Number of clones to export into csv file. This is
#' only relevant in \code{-t all} or \code{-t diversity} where clonotypes
#' are exported into \code{<outdir>/<name>/diversity/clonotypes}
#' @slot detailedComposition logical. Plots composition logo by IGHV families if
#' set to true, otherwise, plots logos by FR/CDRs.
#' @slot log character. Path to log file.
#' @slot merger character. Merger used to merge paired-end reads.
#' @slot fmt character. File format of \code{file1} and (if present) \code{file2}.
#' Possible values are FASTA or FASTQ.
#' @slot sites character. Path to restriction sites \code{txt} file.
#' This option is only used if \code{-t rsasimple}
#' @slot primer5end ANY. Path to 5' end primer FASTA file.
#' @slot primer3end ANY. Path to 3' end primer FASTA file.
#' @slot trim5 numeric. Number of nucleotides to trimd at the 5' end;
#' @slot trim3 numeric. Number of nucleotides to trimd at the 3' end;
#' @slot outdir character. Path to output directory
#' @slot primer5endoffset numeric. Number of nucleotides to offset before aligning
#' 5' end primers in \code{primer5end} FASTA file.
#' @slot threads numeric. Number of threads to run.
#' @slot upstream character. Index (range) of upstream nucleotides to analyze.
#' This option is only used if \code{-t 5utr} or \code{-t secretion}.
#' @slot seqtype character. Sequence type, possible values are either \code{dna}
#' or \code{protein}.
#' @slot database character. Path to IgBLAST database.
#' @slot actualqstart numeric. Query sequence's starting index (indexing starts from 1).
#' This value overrides the inferred query start position by AbSeq.
#' @slot fr4cut logical. The end of FR4 is marked as the end of the sequence if
#' set to TRUE, otherwise the end of the sequence is either the end of the read
#' itself, or trimmed to \code{--trim3 <num>}.
#' @slot domainSystem character. Domain system to use in IgBLAST, possible
#' values are either \code{imgt} or \code{kabat}.
#'
#' @seealso \code{\link{abseqReport}} returns a \code{list} of \code{AbSeqRep}
#' objects.
#' @return AbSeqRep
#' @export
#'
#' @examples
#' # this class is not directly constructed by users, but as a return
#' # value from the abseqReport method.
#'
#' # Use example data from abseqR as abseqPy's output, substitute this
#' # with your own abseqPy output directory
#' abseqPyOutput <- tempdir()
#' file.copy(system.file("extdata", "ex", package = "abseqR"), abseqPyOutput, recursive=TRUE)
#' samples <- abseqReport(file.path(abseqPyOutput, "ex"), report = 0)
#'
#'
#' # gives the name of the first repertoire object returned by abseqReport
#' # samples[[1]]@name
AbSeqRep <- setClass(
    "AbSeqRep",
    slots = c(
        f1 = "character",
        f2 = "character",
        chain = "character",
        task = "character",
        name = "character",
        bitscore = "numeric",
        qstart = "numeric",
        sstart = "numeric",
        alignlen = "numeric",
        clonelimit = "numeric",
        detailedComposition = "logical",
        log = "character",
        merger = "character",
        fmt = "character",
        sites = "character",
        primer5end = "ANY",
        primer3end = "ANY",
        trim5 = "numeric",
        trim3 = "numeric",
        outdir = "character",
        primer5endoffset = "numeric",
        threads = "numeric",
        upstream = "ANY",
        seqtype = "character",
        database = "character",
        actualqstart = "numeric",
        fr4cut = "logical",
        domainSystem = "character"
    )
)


.loadAbSeqRepFromParams <- function(analysisParams) {
    con <- file(analysisParams, "r")
    lines <- readLines(con)
    close(con)
    params <- list(Class = "AbSeqRep")

    # parameters that do not provide any information - ignore them so they
    # wont be added as slots
    skip <- c("yaml")

    for (line in lines) {
        # read in parameters
        if (startsWith(line, "Parameter")) {
            tokens <- unlist(strsplit(line, "\t"))
            parameter <- trimws(strsplit(tokens[1], ":")[[1]][2])
            value <- trimws(strsplit(tokens[2], ":")[[1]][2])
            if (!(parameter %in% skip)) {
                if (!is.na(as.logical(value))) {
                    params[[parameter]] <- as.logical(value)
                } else if (!is.na(suppressWarnings(as.numeric(value)))) {
                    params[[parameter]] <- as.numeric(value)
                } else if (grepl("[", value, fixed = TRUE)) {
                    # value = [start, end]
                    newValue <- gsub("\\[|\\]", "", value)
                    # to -> c(start, end)
                    params[[parameter]] <-
                        as.numeric(unlist(strsplit(newValue, ",")))
                } else {
                    params[[parameter]] <- value
                }
            }
        }
    }
    return(do.call("new", params))
}


#' Combines 2 \linkS4class{AbSeqRep} objects together for comparison
#'
#' @include AbSeqCRep.R
#'
#' @param e1 AbSeqRep object.
#' @param e2 AbSeqRep object.
#'
#' @return \linkS4class{AbSeqCRep} object. Calling \code{abseqR}'s
#' functions on this object will always result in a comparison.
#'
#' @export
#'
#' @seealso \code{\link{abseqReport}} returns a \code{list} of \code{AbSeqRep}s
#'
#' @examples
#' # Use example data from abseqR as abseqPy's output, substitute this
#' # with your own abseqPy output directory
#' abseqPyOutput <- tempdir()
#' file.copy(system.file("extdata", "ex", package = "abseqR"), abseqPyOutput, recursive=TRUE)
#' samples <- abseqReport(file.path(abseqPyOutput, "ex"), report = 0)
#'
#' # The provided example data has PCR1, PCR2, and PCR3 samples contained within
#' # pcr1 and pcr2 are instances of AbSeqRep
#' pcr1 <- samples[["PCR1"]]
#' pcr2 <- samples[["PCR2"]]
#'
#' # pcr12 is an instance of AbSeqCRep
#' pcr12 <- pcr1 + pcr2
#'
#' # you can now call the report function on this object
#' # report(pcr12)           # uncomment this line to execute report
setMethod("+", signature(e1 = "AbSeqRep", e2 = "AbSeqRep"), function(e1, e2) {
    new("AbSeqCRep", repertoires = list(e1, e2))
})

#' Combines a \linkS4class{AbSeqCRep} object with
#' a \linkS4class{AbSeqRep} object together for comparison
#'
#' @include AbSeqCRep.R
#'
#' @param e1 AbSeqCRep.
#' @param e2 AbSeqRep.
#'
#' @return \linkS4class{AbSeqCRep} object. Calling \code{abseqR}'s
#' functions on this object will always result in a comparison.
#'
#' @export
#'
#' @seealso \code{\link{abseqReport}} returns a \code{list} of \code{AbSeqRep}s
#'
#' @examples
#' # Use example data from abseqR as abseqPy's output, substitute this
#' # with your own abseqPy output directory
#' abseqPyOutput <- tempdir()
#' file.copy(system.file("extdata", "ex", package = "abseqR"), abseqPyOutput, recursive=TRUE)
#' samples <- abseqReport(file.path(abseqPyOutput, "ex"), report = 0)
#'
#' # The provided example data has PCR1, PCR2, and PCR3 samples contained within
#' # pcr12 is an instance of AbSeqCRep
#' pcr12 <- samples[["PCR1"]] + samples[["PCR2"]]
#' # pcr3 is instance of AbSeqRep
#' pcr3 <- samples[["PCR3"]]
#'
#' # pcr123 is an instance of AbSeqCRep
#' pcr123 <- pcr12 + pcr3
#'
#' # you can now call the report function on this object
#' # report(pcr123)           # uncomment this line to execute report
setMethod("+", signature(e1 = "AbSeqCRep", e2 = "AbSeqRep"), function(e1, e2) {
    new("AbSeqCRep", repertoires = unique(c(e1@repertoires, e2)))
})

#' Combines a \linkS4class{AbSeqRep} object with
#' a \linkS4class{AbSeqCRep} object together for comparison
#'
#' @include AbSeqCRep.R
#'
#' @param e1 AbSeqRep.
#' @param e2 AbSeqCRep.
#'
#' @return \linkS4class{AbSeqCRep} object. Calling \code{abseqR}'s
#' functions on this object will always result in a comparison.
#'
#' @export
#'
#' @seealso \code{\link{abseqReport}} returns a \code{list} of \code{AbSeqRep}s
#'
#' @examples
#' # Use example data from abseqR as abseqPy's output, substitute this
#' # with your own abseqPy output directory
#' abseqPyOutput <- tempdir()
#' file.copy(system.file("extdata", "ex", package = "abseqR"), abseqPyOutput, recursive=TRUE)
#' samples <- abseqReport(file.path(abseqPyOutput, "ex"), report = 0)
#'
#' # The provided example data has PCR1, PCR2, and PCR3 samples contained within
#' # pcr1 is an instance of AbSeqRep
#' pcr1 <- samples[["PCR1"]]
#' # pcr23 is instance of AbSeqCRep
#' pcr23 <- samples[["PCR2"]] + samples[["PCR3"]]
#'
#' # pcr123 is an instance of AbSeqCRep
#' pcr123 <- pcr1 + pcr23
#'
#' # you can now call the report function on this object
#' # report(pcr123)           # uncomment this line to execute report
setMethod("+", signature(e1 = "AbSeqRep", e2 = "AbSeqCRep"), function(e1, e2) {
    new("AbSeqCRep", repertoires = unique(c(e1, e2@repertoires)))
})

#' Plots \linkS4class{AbSeqRep} or
#' \linkS4class{AbSeqCRep} object to the specfied directory
#'
#' @description Plots all samples in the \code{object} argument
#' and saves the analysis in \code{outputDir}.
#' Users can optionally specify which samples
#' in \code{object} should be compared. Doing so generates
#' additional plots for clonotype comparison and
#' the usual plots will also conveniently include these samples
#' using additional \code{aes}thetics.
#'
#' This method is analogous to \code{\link{abseqReport}}.
#' The only difference is that this method accepts \linkS4class{AbSeqRep} or
#' \linkS4class{AbSeqCRep} objects as its first parameter, and the
#' \code{outputDir} specifies where to store the result.
#'
#'
#' @include util.R
#' @include plotter.R
#' @include AbSeqCRep.R
#'
#' @param object AbSeqRep or AbSeqCRep object to plot.
#' @param outputDir string type. Directory where analysis will be saved to.
#' @param report (optional) integer type. The possible values are:
#' \itemize{
#'   \item{0 - does nothing (returns named list of \linkS4class{AbSeqRep} objects)}
#'   \item{1 - generates plots for csv files}
#'   \item{2 - generates a report that collates all plots}
#'   \item{3 - generates interactive plots in report} (default)
#' }
#' each value also does what the previous values do. For example, \code{report = 2}
#' will return a named list of \linkS4class{AbSeqRep} objects, plot csv files,
#' and generate a (non-interactive)HTML report that collates all the plots together.
#' @return named list. List of \linkS4class{AbSeqRep} objects. The names of
#' the list elements are taken directly from the repertoire object itself.
#' This return value is consistent with the return value of \code{\link{abseqReport}}.
#'
#' @export
#'
#' @seealso \code{\link{abseqReport}}. Analogus function, but takes input from
#' a string that signifies the output directory of abseqPy as the first
#' arugment instead.
#'
#' @seealso \linkS4class{AbSeqRep}
#' @seealso \linkS4class{AbSeqCRep}
#'
#' @rdname report
#'
#' @examples
#' # Use example data from abseqR as abseqPy's output, substitute this
#' # with your own abseqPy output directory
#' abseqPyOutput <- tempdir()
#' file.copy(system.file("extdata", "ex", package = "abseqR"), abseqPyOutput, recursive=TRUE)
#' samples <- abseqReport(file.path(abseqPyOutput, "ex"), report = 0)
#'
#'
#' # The provided example data has PCR1, PCR2, and PCR3 samples contained within
#' # We can use the + operator to combine samples, thus requesting the
#' # report function to compare them:
#' pcr12 <- samples[["PCR1"]] + samples[["PCR2"]]
#'
#' # generate plots and report for this new comparison
#' # report(pcr12, "PCR1_vs_PCR2")
#'
#' # generate plots only
#' # report(pcr12, "PCR1_vs_PCR2", report = 1)
#'
#' # generate plots, and a non-interactive report
#' # report(pcr12, "PCR1_vs_PCR2", report = 2)
#'
#' # generate plots, and an interactive report
#' # report(pcr12, "PCR1_vs_PCR2", report = 3)   # this is the default
setGeneric(
    name = "report",
    def = function(object, outputDir, report = 3) {
        standardGeneric("report")
    }
)



#' @rdname report
setMethod(
    f = "report",
    signature = "AbSeqRep",
    definition = function(object, outputDir, report = 3) {
        skip <- FALSE
        stopifnot(report %in% c(0, 1, 2, 3))

        # give meaningful names to numerical values
        if (report == 0) {
            report <- FALSE
            interactivePlot <- FALSE
            skip <- TRUE
        } else if (report == 1) {
            report <- FALSE
            interactivePlot <- FALSE
        } else if (report == 2) {
            report <- TRUE
            interactivePlot <- FALSE
        } else {
            report <- TRUE
            interactivePlot <- TRUE
        }

        if (!skip) {
            analysisDirectories = c(file.path(object@outdir,
                                              RESULT_DIR, object@name))
            # get analysis that were conducted by looking at the directory
            # structure
            analyses <-
                unlist(.inferAnalyzed(analysisDirectories[1]))
            sampleNames <- c(object@name)
            primer5Files <- list(object@primer5end)
            primer3Files <- list(object@primer3end)
            upstreamRanges <- list(object@upstream)

            .plotSamples(
                sampleNames,
                analysisDirectories,
                analyses,
                outputDir,
                primer5Files,
                primer3Files,
                upstreamRanges,
                skipDgene = (object@chain != "hv")
            )

            if (report) {
                .generateReport(
                    object,
                    root = outputDir,
                    outputDir = outputDir,
                    interactivePlot = interactivePlot
                )
            }
        }

        # to be consistent with the return value of abseqR::abseqReport()
        lst <- list()
        lst[[object@name]] <- object
        return(lst)
    }
)

#' @rdname report
setMethod(
    f = "report",
    signature = "AbSeqCRep",
    definition = function(object, outputDir, report = 3) {
        stopifnot(report %in% c(0, 1, 2, 3))
        skip <- FALSE
        if (report == 0) {
            report <- FALSE
            interactivePlot <- FALSE
            skip <- TRUE
        } else if (report == 1) {
            report <- FALSE
            interactivePlot <- FALSE
        } else if (report == 2) {
            report <- TRUE
            interactivePlot <- FALSE
        } else {
            report <- TRUE
            interactivePlot <- TRUE
        }

        sampleNames <-
            unlist(lapply(object@repertoires, function(x) {
                x@name
            }))

        if (!skip) {
            analysisDirectories <- unlist(lapply(object@repertoires, function(x) {
                # /a/b/c/RESULT_DIR/sample_name has all the analysis folders
                file.path(x@outdir, RESULT_DIR, x@name)
            }))

            # get all analyses conducted by each repertoire by looking
            # at the directory structure created by abseqPy
            analysesConducted <- lapply(analysisDirectories, .inferAnalyzed)

            # find the intersection of all analysis, because it makes no
            # sense to compare samples if they don't have the same analysis
            similarAnalyses <- unlist(Reduce(intersect, analysesConducted))


            primer5Files <-
                unlist(lapply(object@repertoires, function(x) {
                    x@primer5end
                }))
            primer3Files <-
                unlist(lapply(object@repertoires, function(x) {
                    x@primer3end
                }))

            upstreamRanges <- lapply(object@repertoires,
                                     function(x) {
                                         x@upstream
                                     })

            allChains <-
                lapply(object@repertoires, function(x) {
                    x@chain
                })
            # only include D gene plots if all chains only contain "hv"
            skipD <- (("kv" %in% allChains) || ("lv" %in% allChains))

            .plotSamples(
                sampleNames,
                analysisDirectories,
                similarAnalyses,
                outputDir,
                primer5Files,
                primer3Files,
                upstreamRanges,
                skipDgene = skipD
            )

            if (report) {
                .generateReport(
                    object,
                    root = outputDir,
                    outputDir = outputDir,
                    interactivePlot = interactivePlot
                )
            }
        }

        lst <- object@repertoires
        names(lst) <- sampleNames
        return(lst)
    }
)


#' Generates HTML report from \code{\linkS4class{AbSeqRep}} and
#' \code{\linkS4class{AbSeqCRep}} ojects
#'
#' @importFrom rmarkdown pandoc_available render
#' @include util.R
#'
#' @param object AbSeqCRep type.
#' @param root string type. Root directory of the sample(s)
#' @param outputDir string type. The path where the HTML will be generated
#' @param interactivePlot logical type. Interactive or not
#' @param .indexHTML character type. The back button will redirect to this link.
#' This is typically used to redirect users back to index.html page
#'
#' @return path (including HTML name) where the report (HTML file) was saved to
setGeneric(
    name = ".generateReport",
    def = function(object,
                   root,
                   outputDir,
                   interactivePlot = TRUE,
                   .indexHTML = "#") {
        standardGeneric(".generateReport")
    }
)

setMethod(
    f = ".generateReport",
    signature = "AbSeqCRep",
    definition = function(object,
                          root,
                          outputDir,
                          interactivePlot = TRUE,
                          .indexHTML = "#") {
        if (rmarkdown::pandoc_available()) {
            analysisDirectories = unlist(lapply(object@repertoires, function(x) {
                # /a/b/c/RESULT_DIR/sample_name has all the analysis folders
                file.path(x@outdir, RESULT_DIR, x@name)
            }))
            sampleNames <-
                unlist(lapply(object@repertoires, function(x) {
                    x@name
                }))

            message("Generating HTML report for ",
                    paste(sampleNames, collapse = "_vs_"))

            # get all analyses conducted by each repertoire by looking
            # at the directory structure created by abseqPy
            analysesConducted <- lapply(analysisDirectories, .inferAnalyzed)

            # find the intersection of all analysis, because it makes no
            # sense to compare samples if they don't have the same analysis
            similarAnalyses <- unlist(Reduce(intersect, analysesConducted))

            # define variables for template.Rmd's param list
            allChains <-
                lapply(object@repertoires, function(x) {
                    x@chain
                })
            # only include D gene plots if all chains only contain "hv"
            includeD <- !(("kv" %in% allChains) || ("lv" %in% allChains))


            # template.Rmd param list
            bitFilters <-
                lapply(object@repertoires, function(x) {
                    paste(x@bitscore, collapse = " - ")
                })
            alFilters <-
                lapply(object@repertoires, function(x) {
                    paste(x@alignlen, collapse = " - ")
                })
            ssFilters <-
                lapply(object@repertoires, function(x) {
                    paste(x@sstart, collapse = " - ")
                })
            qsFilters <-
                lapply(object@repertoires, function(x) {
                    paste(x@qstart, collapse = " - ")
                })
            rawReadCounts <- lapply(X = analysisDirectories,
                                    FUN = .readSummary,
                                    ABSEQ_RAW_READ_COUNT_KEY)
            annotReadCounts <- lapply(X = analysisDirectories,
                                      FUN = .readSummary,
                                      ABSEQ_ANNOT_READ_COUNT_KEY)
            filteredReadCounts <- lapply(X = analysisDirectories,
                                         FUN = .readSummary,
                                         ABSEQ_FILT_READ_COUNT_KEY)
            productiveReadCounts <- lapply(X = analysisDirectories,
                                           FUN = .readSummary,
                                           ABSEQ_PROD_READ_COUNT_KEY)
            renderParams <- list(
                rootDir = normalizePath(root),
                single = FALSE,
                interactive = interactivePlot,
                inclD = includeD,
                hasAnnot = (ABSEQ_DIR_ANNOT %in% similarAnalyses),
                hasAbun = (ABSEQ_DIR_ABUN %in% similarAnalyses),
                hasProd = (ABSEQ_DIR_PROD %in% similarAnalyses),
                hasDiv = (ABSEQ_DIR_DIV %in% similarAnalyses),
                name = paste(sampleNames, collapse = "_vs_"),
                bitfilters = paste(bitFilters, collapse = ","),
                alignfilters = paste(alFilters, collapse = ","),
                sstartfilters = paste(ssFilters, collapse = ","),
                qstartfilters = paste(qsFilters, collapse = ","),
                rawReads = paste(rawReadCounts, collapse = ","),
                annotReads = paste(annotReadCounts, collapse = ","),
                filteredReads = paste(filteredReadCounts, collapse = ","),
                productiveReads = paste(productiveReadCounts, collapse = ",")
            )

            # https://github.com/rstudio/rmarkdown/issues/861#issuecomment-326637323
            # paraphrased: just don't use output_file or output_dir
            # instead, copy to output_dir and remove it if you want to
            # as a safety measure, rename the template.Rmd file to
            # sample name so that the cache will not clash with other
            # samples it it wasn't cleared in time
            tmpTemplate <-
                file.path(outputDir, paste0(paste(sampleNames, collapse = "_vs_"), ".Rmd"))
            file.copy(
                system.file("extdata", "template.Rmd", package = "abseqR"),
                tmpTemplate,
                overwrite = TRUE
            )
            # substitute href to '#' to index.html: only when there's a landing
            # page that we do this. Calling abseqR::report() will leave href
            # as-is
            .substituteStringInFile(tmpTemplate,
                                    "href: \"#\"",
                                    paste0("href: \"", .indexHTML, "\""),
                                    fixed = TRUE)
            rmarkdown::render(tmpTemplate, params = renderParams)
            file.remove(tmpTemplate)
            return(sub(".Rmd", ".html", tmpTemplate))
        } else {
            warning("Pandoc cannot be detected on system, skipping HTML report")
            return(NA)
        }

    }
)


setMethod(
    f = ".generateReport",
    signature = "AbSeqRep",
    definition = function(object,
                          root,
                          outputDir,
                          interactivePlot = TRUE,
                          .indexHTML = "#") {
        if (rmarkdown::pandoc_available()) {
            message("Generating HTML report for ", object@name)
            analysisDirectory = file.path(object@outdir, RESULT_DIR, object@name)
            analyses <-
                unlist(.inferAnalyzed(analysisDirectory))

            # define parameters for template.Rmd's param list
            bitFilter <-
                paste(object@bitscore, collapse = " - ")
            qsFilter <- paste(object@qstart, collapse = " - ")
            ssFilter <- paste(object@sstart, collapse = " - ")
            alFilter <-
                paste(object@alignlen, collapse = " - ")
            rawReadCount <- .readSummary(analysisDirectory,
                                         ABSEQ_RAW_READ_COUNT_KEY)
            annotReadCount <- .readSummary(analysisDirectory,
                                           ABSEQ_ANNOT_READ_COUNT_KEY)
            filteredReadCount <-
                .readSummary(analysisDirectory,
                             ABSEQ_FILT_READ_COUNT_KEY)
            productiveReadCount <-
                .readSummary(analysisDirectory,
                             ABSEQ_PROD_READ_COUNT_KEY)

            renderParams <- list(
                rootDir = normalizePath(root),
                single = TRUE,
                interactive = interactivePlot,
                inclD = (object@chain == "hv"),
                hasAnnot = (ABSEQ_DIR_ANNOT %in% analyses),
                hasAbun = (ABSEQ_DIR_ABUN %in% analyses),
                hasProd = (ABSEQ_DIR_PROD %in% analyses),
                hasDiv = (ABSEQ_DIR_DIV %in% analyses),
                name = object@name,
                bitfilters = bitFilter,
                qstartfilters = qsFilter,
                sstartfilters = ssFilter,
                alignfilters = alFilter,
                rawReads = rawReadCount,
                annotReads = annotReadCount,
                filteredReads = filteredReadCount,
                productiveReads = productiveReadCount
            )

            # https://github.com/rstudio/rmarkdown/issues/861#issuecomment-326637323
            # paraphrased: just don't use output_file or output_dir
            # instead, copy to output_dir and remove it if you want to
            # as a safety measure, rename the template.Rmd file to
            # sample name so that the cache will not clash with other
            # samples it it wasn't cleared in time
            tmpTemplate <-
                file.path(outputDir, paste0(object@name, ".Rmd"))
            file.copy(
                system.file("extdata", "template.Rmd", package = "abseqR"),
                tmpTemplate,
                overwrite = TRUE
            )
            .substituteStringInFile(tmpTemplate,
                                    "href: \"#\"",
                                    paste0("href: \"", .indexHTML, "\""),
                                    fixed = TRUE)
            rmarkdown::render(tmpTemplate, params = renderParams)
            file.remove(tmpTemplate)
            return(sub(".Rmd", ".html", tmpTemplate))
        } else {
            warning("Pandoc cannot be detected on system, skipping HTML report")
            return(NA)
        }
    }
)
