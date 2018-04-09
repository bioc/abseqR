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
#'
#' @return None
.annotAnalysis <- function(annotDirectories, annotOut,
                          sampleNames, mashedNames) {

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
            ylabel = "Proportion"
        )
        ggsave(file.path(annotOut, paste0(mashedNames,
                      "_all_clones_len_dist.png")), plot = g,
               width = V_WIDTH, height = V_HEIGHT)
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
            ylabel = "Proportion"
        )
        ggsave(file.path(annotOut, paste0(mashedNames,
                      "_all_clones_len_dist_no_outliers.png")),
               plot = g, width = V_WIDTH, height = V_HEIGHT)
    } else {
        warning(paste(
            "Cannot find clone length (no outliers)",
            "distribution file from samples",
            paste(sampleNames, collapse = ", ")
        ))
    }
}
