#' AbSeq analysis object.
#'
#' @description The AbSeqRep object contains all metadata associated with the AbSeq (python backend)
#' run conducted on it. For further information, refer to AbSeq's python help.
#'
#' @slot f1 character. Path to FASTA/FASTQ file 1.
#' @slot f2 character. Path to FASTA/FASTQ file 2.
#' @slot chain character. Type of chain, possible values:
#' \itemize{
#'   \item{hv}
#'   \item{lv}
#'   \item{kv}
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
#' @slot primer numeric. Dummy value - not implemented yet.
#' @seealso \code{\link{abseqReport}} returns a \code{list} of \code{AbSeqRep}
#' objects.
#' @return none
#' @export
#'
#' @examples
#' \dontrun{
#' # this class is (usually) not directly constructed by users, but as a return
#' # value from the abseqReport method.
#' samples <- abseqReport("/path/to/output/directory/")
#' samples[[1]]@name     # gives the name of the first repertoire object returned by abseqReport
#' }
AbSeqRep <- setClass("AbSeqRep", slots = c(
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
    domainSystem = "character",
    primer = "numeric"
))


.loadAbSeqRepFromParams <- function(analysisParams) {
    con <- file(analysisParams, "r")
    lines <- readLines(con)
    close(con)
    params <- list(Class = "AbSeqRep")
    skip <- c("report_interim", "yaml")

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
                } else if (grepl("[", value, fixed = T)) {
                    # value = [start, end]
                    newValue <- gsub("\\[|\\]", "", value)
                    # to -> c(start, end)
                    params[[parameter]] <- as.numeric(unlist(strsplit(newValue, ",")))
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
#' \dontrun{
#' # 'load' AbSeqRep objects using abseqReport (ignoring plots, report, etc..)
#' # also assumes there's a result/ directory in current working directory,
#' # where result/ is the same argument passed to abseqPy's -o or --outdir parameter
#' samples <- abseqReport("results")
#'
#' # assuming there are samples named "Sample1" and "Sample3"
#' S1S3 <- samples[["Sample1"]] + samples[["Sample3"]]
#'
#' # generate plots and report for this new comparison
#' report(S1S3, "s1_vs_s3")
#' }
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
#' \dontrun{
#' # 'load' AbSeqRep objects using abseqReport (ignoring plots, report, etc..)
#' # also assumes there's a result/ directory in current working directory,
#' # where result/ is the same argument passed to abseqPy's -o or --outdir parameter
#' samples <- abseqReport("results")
#'
#' # assuming there are samples named "Sample1", "Sample3", and "Sample4"
#' S1S3 <- samples[["Sample1"]] + samples[["Sample3"]]
#' S4 <- samples[["Sample4"]]
#'
#' S1S3S4 <- S1S3 + S4
#'
#' # generate plots and report for this new comparison
#' report(S1S3S4, "s1_vs_s3_vs_s4")
#' }
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
#' \dontrun{
#' # 'load' AbSeqRep objects using abseqReport (ignoring plots, report, etc..)
#' # also assumes there's a result/ directory in current working directory,
#' # where result/ is the same argument passed to abseqPy's -o or --outdir parameter
#' samples <- abseqReport("results")
#'
#' # assuming there are samples named "Sample1", "Sample3", and "Sample4"
#' S1S3 <- samples[["Sample1"]] + samples[["Sample3"]]
#' S4 <- samples[["Sample4"]]
#'
#' S4S1S3 <- S4 + S1S3
#'
#' # generate plots and report for this new comparison
#' report(S4S1S3, "s4_vs_s1_vs_s3")
#' }
setMethod("+", signature(e1 = "AbSeqRep", e2 = "AbSeqCRep"), function(e1, e2) {
    new("AbSeqCRep", repertoires = unique(c(e1, e2@repertoires)))
})

#' Plots \linkS4class{AbSeqRep} or
#' \linkS4class{AbSeqCRep} object to the specfied directory
#'
#' @import rmarkdown
#'
#' @include util.R
#' @include plotter.R
#' @include AbSeqCRep.R
#'
#' @param object AbSeqRep or AbSeqCRep object to plot
#' @param outputDir string type. where to save the files to
#' @param report logical type. Should HTML report be generated
#' @param interactivePlot logical type. Should the HTML report's plots be
#' interactive (using plotly)? This argument will be ignored if \code{report}
#' is \code{FALSE}
#'
#' @return nothing
#' @export
#'
#' @seealso \code{\link{abseqReport}} returns a \code{list} of \code{AbSeqRep}s
#' @seealso \linkS4class{AbSeqRep}
#' @seealso \linkS4class{AbSeqCRep}
#'
#' @examples todo
setGeneric(name = "report",
           def = function(object, outputDir, report = TRUE, interactivePlot = TRUE) {
               standardGeneric("report")
           })



setMethod(f = "report",
          signature = "AbSeqRep",
          definition = function(object, outputDir, report = TRUE, interactivePlot = TRUE) {
              if (!report && interactivePlot) {
                  warning("report is FALSE, ignoring interactivePlot argument")
              }

              analysisDirectories = c(file.path(object@outdir,
                                                RESULT_DIR, object@name))
              # get analysis that were conducted by looking at the directory
              # structure - is this the best way?
              analyses <- unlist(.inferAnalyzed(analysisDirectories[1]))
              sampleNames <- c(object@name)
              primer5Files <- list(object@primer5end)
              primer3Files <- list(object@primer3end)
              upstreamRanges <- list(object@upstream)

              .plotSamples(sampleNames, analysisDirectories, analyses, outputDir,
                           primer5Files, primer3Files, upstreamRanges,
                           skipDgene = (object@chain != "hv"))

              if (report) {
                  .generateReport(object, root = outputDir,
                                  outputDir = outputDir,
                                  interactivePlot = interactivePlot)
              }
          })

setMethod(f = "report",
          signature = "AbSeqCRep",
          definition = function(object, outputDir, report = TRUE, interactivePlot = TRUE) {
              if (!report && interactivePlot) {
                  warning("report is FALSE, ignoring interactivePlot argument")
              }
              analysisDirectories = unlist(lapply(object@repertoires, function(x) {
                  # /a/b/c/RESULT_DIR/sample_name has all the analysis folders
                  file.path(x@outdir, RESULT_DIR, x@name)
              }))

              # get all analyses conducted by each repertoire by looking
              # at the directory structure created by AbSeqPy
              analysesConducted <- lapply(analysisDirectories, .inferAnalyzed)

              # find the intersection of all analysis, because it makes no
              # sense to compare samples if they don't have the same analysis
              similarAnalyses <- unlist(Reduce(intersect, analysesConducted))

              sampleNames <- unlist(lapply(object@repertoires, function(x) {
                  x@name
              }))

              primer5Files <- unlist(lapply(object@repertoires, function(x) {
                  x@primer5end
              }))
              primer3Files <- unlist(lapply(object@repertoires, function(x) {
                  x@primer3end
              }))

              upstreamRanges <- lapply(object@repertoires,
                                       function(x) {
                                           x@upstream
                                       })

              allChains <- lapply(object@repertoires, function(x) { x@chain })
              # only include D gene plots if all chains only contain "hv"
              skipD <- (("kv" %in% allChains) || ("lv" %in% allChains))

              .plotSamples(sampleNames, analysisDirectories, similarAnalyses,
                           outputDir, primer5Files, primer3Files, upstreamRanges,
                           skipDgene = skipD)

              if (report) {
                  .generateReport(object, root = outputDir,
                                  outputDir = outputDir,
                                  interactivePlot = interactivePlot)
              }
          })


#' Title
#'
#' @import rmarkdown
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
setGeneric(name = ".generateReport",
           def = function(object, root, outputDir, interactivePlot = TRUE, .indexHTML = "#") {
               standardGeneric(".generateReport")
           })

setMethod(f = ".generateReport",
          signature = "AbSeqCRep",
          definition = function(object, root, outputDir, interactivePlot = TRUE, .indexHTML = "#") {
              if (rmarkdown::pandoc_available()) {
                  analysisDirectories = unlist(lapply(object@repertoires, function(x) {
                      # /a/b/c/RESULT_DIR/sample_name has all the analysis folders
                      file.path(x@outdir, RESULT_DIR, x@name)
                  }))
                  sampleNames <- unlist(lapply(object@repertoires, function(x) {
                      x@name
                  }))
                  message(paste("Generating HTML report for", paste(sampleNames, collapse = "_vs_")))

                  # get all analyses conducted by each repertoire by looking
                  # at the directory structure created by AbSeqPy
                  analysesConducted <- lapply(analysisDirectories, .inferAnalyzed)

                  # find the intersection of all analysis, because it makes no
                  # sense to compare samples if they don't have the same analysis
                  similarAnalyses <- unlist(Reduce(intersect, analysesConducted))

                  # define variables for template.Rmd's param list
                  allChains <- lapply(object@repertoires, function(x) {
                      x@chain
                  })
                  # only include D gene plots if all chains only contain "hv"
                  includeD <- !(("kv" %in% allChains) || ("lv" %in% allChains))


                  # template.Rmd param list
                  bitFilters <- lapply(object@repertoires, function(x) {
                      paste(x@bitscore, collapse = " - ")
                  })
                  alFilters <- lapply(object@repertoires, function(x) {
                      paste(x@alignlen, collapse = " - ")
                  })
                  ssFilters <- lapply(object@repertoires, function(x) {
                      paste(x@sstart, collapse = " - ")
                  })
                  qsFilters <- lapply(object@repertoires, function(x) {
                      paste(x@qstart, collapse = " - ")
                  })
                  rawReadCounts <- lapply(analysisDirectories, function(pth) {
                      .readSummary(pth, ABSEQ_RAW_READ_COUNT_KEY)
                  })
                  annotReadCounts <- lapply(analysisDirectories, function(pth) {
                      .readSummary(pth, ABSEQ_ANNOT_READ_COUNT_KEY)
                  })
                  filteredReadCounts <- lapply(analysisDirectories, function(pth) {
                      .readSummary(pth, ABSEQ_FILT_READ_COUNT_KEY)
                  })
                  productiveReadCounts <- lapply(analysisDirectories, function(pth) {
                      .readSummary(pth, ABSEQ_PROD_READ_COUNT_KEY)
                  })

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
                  tmpTemplate <- file.path(outputDir, paste0(paste(sampleNames, collapse = "_vs_"), ".Rmd"))
                  file.copy(system.file("extdata", "template.Rmd", package = "abseqR"), tmpTemplate, overwrite = T)
                  .substituteStringInFile(tmpTemplate, "href: \"#\"", paste0("href: \"", .indexHTML, "\""), fixed = T)
                  rmarkdown::render(tmpTemplate, params = renderParams)
                  file.remove(tmpTemplate)
                  return(sub(".Rmd", ".html", tmpTemplate))
              } else {
                  warning("Pandoc cannot be detected on system, skipping HTML report")
                  return(NA)
              }

          })


setMethod(f = ".generateReport",
          signature = "AbSeqRep",
          definition = function(object, root, outputDir, interactivePlot = TRUE, .indexHTML = "#") {
              if (rmarkdown::pandoc_available()) {
                  message(paste("Generating HTML report for", object@name))
                  analysisDirectory = file.path(object@outdir, RESULT_DIR, object@name)
                  analyses <- unlist(.inferAnalyzed(analysisDirectory))

                  # define parameters for template.Rmd's param list
                  bitFilter <- paste(object@bitscore, collapse = " - ")
                  qsFilter <- paste(object@qstart, collapse = " - ")
                  ssFilter <- paste(object@sstart, collapse = " - ")
                  alFilter <- paste(object@alignlen, collapse = " - ")
                  rawReadCount <- .readSummary(analysisDirectory,
                                               ABSEQ_RAW_READ_COUNT_KEY)
                  annotReadCount <- .readSummary(analysisDirectory,
                                                 ABSEQ_ANNOT_READ_COUNT_KEY)
                  filteredReadCount <- .readSummary(analysisDirectory,
                                                    ABSEQ_FILT_READ_COUNT_KEY)
                  productiveReadCount <- .readSummary(analysisDirectory,
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
                  tmpTemplate <- file.path(outputDir, paste0(object@name, ".Rmd"))
                  file.copy(system.file("extdata", "template.Rmd", package = "abseqR"), tmpTemplate, overwrite = T)
                  .substituteStringInFile(tmpTemplate, "href: \"#\"", paste0("href: \"", .indexHTML, "\""), fixed = T)
                  rmarkdown::render(tmpTemplate, params = renderParams)
                  file.remove(tmpTemplate)
                  return(sub(".Rmd", ".html", tmpTemplate))
              } else {
                  warning("Pandoc cannot be detected on system, skipping HTML report")
                  return(NA)
              }
         })
