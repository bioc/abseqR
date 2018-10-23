#' Accessor for the \code{name} slot
#'
#' @param object AbSeqRep object
#'
#' @return character type, the sample name of this object.
.asRepertoireName <- function(object) {
    object@name
}


#' Accessor for the \code{outdir} slot
#'
#' @param object AbSeqRep object
#'
#' @return character type, the output directory of this object
.asRepertoireDir <- function(object) {
    object@outdir
}


#' Accessor for \linkS4class{AbSeqCRep}'s list of \linkS4class{AbSeqRep} objects
#'
#' @param object AbSeqCRep object
#'
#' @return list type, list of \linkS4class{AbSeqRep} objects that together,
#' compose this \linkS4class{AbSeqCRep} object.
.asRepertoireList <- function(object) {
    object@repertoires
}


#' Accessor for \code{chain} slot
#'
#' @param object AbSeqRep object
#'
#' @return character type, the chain type of this sample
.asRepertoireChain <- function(object) {
    object@chain
}


#' Accessor for \code{bitscore} slot
#'
#' @param object AbSeqRep object
#' @param collapse character type, collapse the range using this string.
#'
#' @return character type. If collapse is a string, then the ranges are represented
#' as `start - end` in a string, if collapse is NULL,
#' returns a character vector of length
#' two, denoting the start and end value respectively.
.asRepertoireBitscore <- function(object, collapse = " - ") {
    if (is.null(collapse)) {
        return(object@bitscore)
    }
    return(paste(object@bitscore, collapse = collapse))
}


#' Accessor for \code{qstart} slot
#'
#' @param object AbSeqRep object
#' @param collapse character type, collapse the range using this string.
#'
#' @return character type. If collapse is a string, then the ranges
#' are represented as `start - end` in a string, if collapse is NULL,
#' returns a character vector of length
#' two, denoting the start and end value respectively.
.asRepertoireQueryStart <- function(object, collapse = " - ") {
    if (is.null(collapse)) {
        return(object@qstart)
    }
    return(paste(object@qstart, collapse = collapse))
}


#' Accessor for \code{alignlen} slot
#'
#' @param object AbSeqRep object
#' @param collapse character type, collapse the range using this string.
#'
#' @return character type. If collapse is a string, then the ranges
#' are represented as `start - end` in a string, if collapse is NULL,
#' returns a character vector of length
#' two, denoting the start and end value respectively.
.asRepertoireAlignLen <- function(object, collapse = " - ") {
    if (is.null(collapse)) {
        return(object@alignlen)
    }
    return(paste(object@alignlen, collapse = collapse))
}


#' Accessor for \code{sstart} slot
#'
#' @param object AbSeqRep object
#' @param collapse character type, collapse the range using this string.
#'
#' @return character type. If collapse is a string, then the ranges
#' are represented as `start - end` in a string, if collapse is NULL,
#' returns a character vector of length
#' two, denoting the start and end value respectively.
.asRepertoireSubjectStart <- function(object, collapse = " - ") {
    if (is.null(collapse)) {
        return(object@sstart)
    }
    return(paste(object@sstart, collapse = collapse))
}


#' Accessor for the \code{primer5end} slot
#'
#' @param object AbSeqRep object
#'
#' @return character type, the FASTA file name for primer 5' end sequences
.asRepertoirePrimer5 <- function(object) {
    object@primer5end
}

#' Accessor for the \code{primer3end} slot
#'
#' @param object AbSeqRep object
#'
#' @return character type, the FASTA file name for primer 3' end sequences
.asRepertoirePrimer3 <- function(object) {
    object@primer3end
}

#' Accessor for the \code{upstream} slot
#'
#' @param object AbSeqRep object
#'
#' @return character type
.asRepertoireUpstream <- function(object) {
    object@upstream
}
