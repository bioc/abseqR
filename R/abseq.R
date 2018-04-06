
#' AbSeq S4 class
#'
#' @slot outputDirectory character type. Output directory of AbSeq analysis
#'
#' @return none
#' @export AbSeq
#'
#' @examples
AbSeq <- setClass("AbSeq",
                  slots = c(
                      outputDirectory = "character"
                  ))
