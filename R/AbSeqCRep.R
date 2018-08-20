#' AbSeqCompositeRepertoire
#'
#' @description  AbSeqCRep is a collection of \linkS4class{AbSeqRep} S4 objects
#'
#' @slot repertoires list.
#'
#' @return AbSeqCRep
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
                  repertoires = unique(c(e1@repertoires, e2@repertoires)))

          })
