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
#' \dontrun{
#' # 'load' AbSeqRep objects using abseqReport (ignoring plots, report, etc..)
#' # also assumes there's a result/ directory in current working directory,
#' # where result/ is the same argument passed to abseqPy's -o or --outdir parameter
#' samples <- abseqReport("results", report = 0)
#'
#' # assuming there are samples named "Sample1", "Sample2", "Sample3", and "Sample4"
#'
#' # S1S3 is an instance of AbSeqCRep
#' S1S3 <- samples[["Sample1"]] + samples[["Sample3"]]
#'
#' # S2S4 is an instance of AbSeqCRep
#' S2S4 <- samples[["Sample2"]] + samples[["Sample4"]]
#'
#' # all.S is an instance of AbSeqCRep
#' all.S <- S1S3 + S2S4
#' }
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
#' \dontrun{
#' # 'load' AbSeqRep objects using abseqReport (ignoring plots, report, etc..)
#' # also assumes there's a result/ directory in current working directory,
#' # where result/ is the same argument passed to abseqPy's -o or --outdir parameter
#' samples <- abseqReport("results", report = 0)
#'
#' # assuming there are samples named "Sample1", "Sample2", "Sample3", and "Sample4"
#' S1S3 <- samples[["Sample1"]] + samples[["Sample3"]]
#' S2S4 <- samples[["Sample2"]] + samples[["Sample4"]]
#'
#' all.S <- S1S3 + S2S4
#'
#' # generate plots and report for this new comparison
#' report(all.S, "s1_vs_s2_vs_s3_vs_s4")
#' }
setMethod("+",
          signature(e1 = "AbSeqCRep",
                    e2 = "AbSeqCRep"),
          function(e1, e2) {
              new("AbSeqCRep",
                  repertoires = unique(c(e1@repertoires, e2@repertoires)))

          })
