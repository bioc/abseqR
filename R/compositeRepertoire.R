#' Title
#'
#' @slot repertoires list.
#'
#' @return
#' @export
#'
#' @examples
CompositeRepertoire <-
    setClass("CompositeRepertoire", slots = c(repertoires = "list"))


#' Combines 2 \linkS4class{CompositeRepertoire} objects together for comparison
#'
#' @param e1 CompositeRepertoire.
#' @param e2 CompositeRepertoire.
#'
#' @return \linkS4class{CompositeRepertoire} object. Calling \code{abseqR}'s
#' functions on this object will always result in a comparison.
#'
#' @export
#'
#' @seealso \code{\link{abseqReport}} returns a \code{list} of \code{Repertoire}s
#'
#' @examples
#' \dontrun{
#' # 'load' Repertoire objects using abseqReport (ignoring plots, report, etc..)
#' # also assumes there's a result/ directory in current working directory,
#' # where result/ is the same argument passed to abseqPy's -o or --outdir parameter
#' samples <- abseqReport("results")
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
          signature(e1 = "CompositeRepertoire",
                    e2 = "CompositeRepertoire"),
          function(e1, e2) {
              new("CompositeRepertoire",
                  repertoires = unique(c(e1@repertoires, e2@repertoires)))

          })
