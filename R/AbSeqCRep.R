#' S4 class - AbSeqCompositeRepertoire analysis object
#'
#' @description  AbSeqCRep is a collection of \linkS4class{AbSeqRep} S4 objects.
#' This S4 class contains multiple samples(repertoires) and it can be
#' "combined" with other samples by using the \code{+} operator to
#' create an extended \linkS4class{AbSeqCRep} object.
#' This value, in turn, can be used as the first argument to
#' \link{report} which generates a comparison between all samples included
#' in the \code{+} operation.
#'
#' Users do not manually construct this class, but rather indirectly obtain
#' this class object as a return value from the \code{+} operation between two
#' \linkS4class{AbSeqRep} objects, which are in turn, obtained indirectly from
#' \link{abseqReport} and \link{report} functions. It is also possible to
#' obtain this class object by \code{+} (adding) \linkS4class{AbSeqCRep} objects.
#'
#' @slot repertoires list of \linkS4class{AbSeqRep} objects.
#'
#' @return AbSeqCRep
#'
#' @export
#'
#' @seealso \linkS4class{AbSeqRep}
#'
#' @examples
#' # Use example data from abseqR as abseqPy's output, substitute this
#' # with your own abseqPy output directory
#' abseqPyOutput <- tempdir()
#' file.copy(system.file("extdata", "ex", package = "abseqR"), abseqPyOutput, recursive=TRUE)
#' samples <- abseqReport(file.path(abseqPyOutput, "ex"), report = 0)
#'
#' # The provided example data has PCR1, PCR2, and PCR3 samples contained within
#' # pcr12 and pcr13 are instances of AbSeqCRep
#' pcr12 <- samples[["PCR1"]] + samples[["PCR2"]]
#' pcr13 <- samples[["PCR1"]] + samples[["PCR3"]]
#'
#' # all_S is also an instance of AbSeqCRep
#' all_S <- pcr12 + pcr13
AbSeqCRep <-
    setClass("AbSeqCRep", slots = c(repertoires = "list"))


#' Combines 2 \linkS4class{AbSeqCRep} objects together for comparison
#'
#' @include accessors-AbSeq.R
#'
#' @param e1 AbSeqCRep.
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
#' # pcr12 and pcr13 are instances of AbSeqCRep
#' pcr12 <- samples[["PCR1"]] + samples[["PCR2"]]
#' pcr13 <- samples[["PCR1"]] + samples[["PCR3"]]
#'
#' # all_S is also an instance of AbSeqCRep
#' all_S <- pcr12 + pcr13
#'
#' # you can now call the report function on this object
#' # report(all_S)           # uncomment this line to execute report
setMethod("+",
          signature(e1 = "AbSeqCRep",
                    e2 = "AbSeqCRep"),
          function(e1, e2) {
              new("AbSeqCRep",
                  repertoires = unique(c(.asRepertoireList(e1),
                                         .asRepertoireList(e2))))

          })
