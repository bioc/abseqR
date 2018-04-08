#' Title
#'
#' @slot f1 character.
#' @slot f2 character.
#' @slot chain character.
#' @slot task character.
#' @slot name character.
#' @slot bitscore numeric.
#' @slot qstart numeric.
#' @slot sstart numeric.
#' @slot alignlen numeric.
#' @slot clonelimit numeric.
#' @slot log character.
#' @slot merger character.
#' @slot fmt character.
#' @slot sites character.
#' @slot primer5end ANY.
#' @slot primer3end ANY.
#' @slot trim5 numeric.
#' @slot trim3 numeric.
#' @slot outdir character.
#' @slot primer5endoffset numeric.
#' @slot threads numeric.
#' @slot upstream character.
#' @slot seqtype character.
#' @slot database character.
#' @slot actualqstart numeric.
#' @slot fr4cut logical.
#' @slot domainSystem character.
#' @slot primer numeric.
#'
#' @return
#' @export
#'
#' @examples
Repertoire <- setClass("Repertoire", slots = c(
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


.loadRepertoireFromParams <- function(analysisParams) {
    con <- file(analysisParams, "r")
    lines <- tail(readLines(con), n = -2)
    close(con)
    params <- list(Class = "Repertoire")
    skip <- c("report_interim", "rscripts")
    for (line in lines) {
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
    return(do.call("new", params))
}


#' Title
#'
#' @include compositeRepertoire.R
#'
#' @param e1 Repertoire.
#' @param e2 Repertoire.
#'
#' @return
#' @export
#'
#' @examples
setMethod("+", signature(e1 = "Repertoire", e2 = "Repertoire"), function(e1, e2) {
    new("CompositeRepertoire", repertoires = list(e1, e2))
})

#' Title
#'
#' @include compositeRepertoire.R
#'
#' @param e1 CompositeRepertoire.
#' @param e2 Repertoire.
#'
#' @return
#' @export
#'
#' @examples
setMethod("+", signature(e1 = "CompositeRepertoire", e2 = "Repertoire"), function(e1, e2) {
    new("CompositeRepertoire", repertoires = unique(c(e1@repertoires, e2)))
})

#' Title
#'
#' @include compositeRepertoire.R
#'
#' @param e1 Repertoire.
#' @param e2 CompositeRepertoire.
#'
#' @return
#' @export
#'
#' @examples
setMethod("+", signature(e1 = "Repertoire", e2 = "CompositeRepertoire"), function(e1, e2) {
    new("CompositeRepertoire", repertoires = unique(c(e1, e2@repertoires)))
})
