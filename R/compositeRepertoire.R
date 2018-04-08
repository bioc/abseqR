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


#' Title
#'
#' @param e1 CompositeRepertoire.
#' @param e2 CompositeRepertoire.
#'
#' @return
#' @export
#'
#' @examples
setMethod("+",
          signature(e1 = "CompositeRepertoire",
                    e2 = "CompositeRepertoire"),
          function(e1, e2) {
              new("CompositeRepertoire",
                  repertoires = unique(c(e1@repertoires, e2@repertoires)))

          })
