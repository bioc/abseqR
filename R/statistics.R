#' Conducts pearson and spearman correlation analysis on dataframe
#'
#' @param df dataframe with at least the following 2 columns:
#' +-----------------+
#' | prop.x | prop.y |
#' +-----------------+
#' |.....   | ....   |
#' +-----------------+
#' where prop.x and prop.y are normalized counts (i.e. frequencies) of the clones
#' They may contain 0 in a column to denote it being missing from sample x or y.
#' @return named list of pearson, pearson.p, spearman, spearman.p
.correlationTest <- function(df) {
    df <- df[, c("prop.x", "prop.y")]
    pearson <- cor.test(df$prop.x, df$prop.y, method = "pearson")
    spearman <- cor.test(df$prop.x, df$prop.y, method = "spearman")
    # return values:
    # [1] "statistic"   "parameter"   "p.value"     "estimate"    "null.value"
    # [6] "alternative" "method"      "data.name"
    list(
        pearson = as.numeric(pearson$estimate),
        pearson.p = as.numeric(pearson$p.value),
        spearman = as.numeric(spearman$estimate),
        spearman.p = as.numeric(spearman$p.value)
    )
}


#' Computes the distance between pariwise samples
#'
#' @param df dataframe with at least the following 2 columns:
#' +-----------------+
#' | prop.x | prop.y |
#' +-----------------+
#' |.....   | ....   |
#' +-----------------+
#' where prop.x and prop.y are normalized counts (i.e. frequencies) of the clones
#' They may contain 0 in a column to denote it being missing from sample x or y.
#'
#' @return named list of bray.curtis, jaccard, and morisita.horn
.distanceMeasure <- function(df) {
    # at the very minimum
    # jaccard
    # sorensen
    # morisita-horn

    df.prop <- t(df[, c("prop.x", "prop.y")])

    # sorensen:
    bray.curtis <- as.numeric(vegdist(df.prop, method = "bray"))
    # jaccard: see ?vegandist to understand this:
    jaccard <- (2 * (bray.curtis)) / (1 + bray.curtis)
    # morisita-horn:
    morisita.horn <- as.numeric(vegdist(df.prop, method = "horn"))

    list(
        bray.curtis = bray.curtis,
        jaccard = jaccard,
        morisita.horn = morisita.horn
    )

}
