#' Annotation analysis
#'
#' @import ggplot2
#' @include util.R
#' @include distributions.R
#'
#' @param annotDirectories list type. List of sample directories
#' @param annotOut string type. Output directory
#' @param sampleNames vector type. 1-1 with annotDirectories
#' @param mashedNames string type. File output "mashed" sample names
#' @param .save logical type. Saves ggplot object
#'
#' @return none
.annotAnalysis <- function(annotDirectories, annotOut,
                          sampleNames, mashedNames, .save = TRUE) {

    # with outliers
    searchFiles <-
        .listFilesInOrder(path = annotDirectories,
                          pattern = ".*_all_clones_len_dist\\.csv(\\.gz)?$")
    message("Starting annotation plot")
    # plotSpectratype names the sample(s) for us in the title
    if (length(searchFiles) > 0) {
        g <- .plotSpectratype(
            lapply(searchFiles, read.csv),
            sampleNames,
            title = "Sequence lengths",
            xlabel = "Sequence Length (bp)",
            ylabel = "Proportion",
            showLabel = F
        )
        fname <- file.path(annotOut,
                           paste0(mashedNames, "_all_clones_len_dist.png"))
        ggsave(fname, plot = g, width = V_WIDTH, height = V_HEIGHT)
        .saveAs(.save, fname, g)
    } else {
        warning(paste("Cannot find clone length distribution file from samples",
                      paste(sampleNames, collapse = ", ")))
    }

    # without outliers
    searchFiles <-
        .listFilesInOrder(path = annotDirectories,
                          pattern = ".*_all_clones_len_dist_no_outliers\\.csv(\\.gz)?$")

    # plotSpectratype names the sample(s) for us in the title
    if (length(searchFiles) > 0) {
        g <- .plotSpectratype(
            lapply(searchFiles, read.csv),
            sampleNames,
            title = "Sequence lengths",
            xlabel = "Sequence Length (bp)",
            ylabel = "Proportion",
            showLabel = F
        )
        fname <- file.path(annotOut,
                           paste0(mashedNames,
                                  "_all_clones_len_dist_no_outliers.png"))
        ggsave(fname, plot = g, width = V_WIDTH, height = V_HEIGHT)
        .saveAs(.save, fname, g)
    } else {
        warning(paste(
            "Cannot find clone length (no outliers)",
            "distribution file from samples",
            paste(sampleNames, collapse = ", ")
        ))
    }
}
