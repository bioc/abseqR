#' Title Abundance distribution
#'
#' @import ggplot2
#' @include util.R
#'
#' @param files list type. list of files in abundance directory
#' @param sampleNames vector type. 1-1 correspondance to files
#' @param outputDir string type
#'
#' @return None
#'
#' @examples
.abundancePlot <- function(files, sampleNames, outputDir) {
    for (expression in c("family", "gene")) {
        for (gene in c("v", "d", "j")) {
            # correction, J has no "gene" but rather variant
            if (gene == 'j' && expression == 'gene') {
                expression <- 'variant'
            }
            reg <- paste0(".*_ig", gene, "_dist_",
                          expression, "_level\\.csv(\\.gz)?$")
            selectedFiles <- files[grepl(reg, files)]
            if (length(selectedFiles)) {
                vert <- .checkVert(selectedFiles[1])
                if (vert) {
                    width <- V_WIDTH
                    height <- V_HEIGHT
                } else {
                    width <- H_WIDTH
                    height <- H_HEIGHT
                }
                mashedName <- paste(sampleNames, collapse = ", ")
                dataframes <- lapply(selectedFiles, read.csv,
                                     stringsAsFactors = FALSE, skip = 1)
                subtitle <- paste("Total is",
                                  paste(lapply(selectedFiles, function(x) {
                                        as.integer(.getTotal(x))
                                  }), collapse = ", "))
                p <- .plotDist(dataframes, sampleNames,
                               paste0("IG", toupper(gene),
                                      " abundance in ", mashedName),
                               vert, subs = subtitle)
                ggsave(paste0(outputDir,
                              paste(sampleNames, collapse = "_"), "_ig",
                              gene, "_dist_", expression, "_level.png"),
                       plot = p, width = width, height = height)
            }
        }
    }
}

#' Title
#'
#' @import circlize
#' @import RColorBrewer
#'
#' @param sampleName
#' @param path
#'
#' @return
#'
#' @examples
.plotCirclize <- function(sampleName, path) {
    filename <- paste0(path, sampleName, "_vjassoc.csv")

    if (file.exists(filename)) {
        df <- read.csv(filename)

        # output file
        png(gsub(".csv", ".png", filename))

        # circos theme setup
        #if (length(unique(df[[1]]))-1 < 8 && length(unique(df[[2]]))-1 < 8)  {
        #    circos.par(gap.after = c(rep(5, length(unique(df[[1]]))-1), 15,
        #                             rep(5, length(unique(df[[2]]))-1), 15))
        #}

        row = rep(brewer.pal(12, "Paired"), nrow(df))[1:length(unique(df[[1]]))]
        col = rep(rev(brewer.pal(12, "Paired")),
                  nrow(df))[1:length(unique(df[[2]]))]

        chordDiagram(df, annotationTrack = "grid",
                     preAllocateTracks = list(track.height = 0.2),
                     grid.col = c(row,col))
        title(sampleName, cex = 0.8)
        circos.trackPlotRegion(track.index = 1, bg.border = NA,
                               panel.fun = function(x, y) {
                                   sector.name =
                                       get.cell.meta.data("sector.index")
                                   xlim = get.cell.meta.data("xlim")
                                   ylim = get.cell.meta.data("ylim")
                                   circos.text(mean(xlim), ylim[1],
                                               sector.name,
                                               facing = "clockwise",
                                               adj = c(0, 1.5))
                               }
        )
        circos.clear()
        dev.off()
    }
}


#' Title
#'
#' @import ggplot2
#' @include util.R
#'
#' @param abundanceDirectories
#' @param abunOut
#' @param sampleNames
#' @param combinedNames
#' @param mashedNames
#'
#' @return
#'
#' @examples
.abundanceAnalysis <- function(abundanceDirectories, abunOut,
                               sampleNames, combinedNames, mashedNames) {
    # where to find the files
    # 3 each from V and D, then 2 from J (no gene) or 5 (exclude the 3 from D)
    searchFiles <-
        .listFilesInOrder(path = abundanceDirectories,
                          pattern = ".*ig[vdj]_dist_[family|gene|variant].*\\.csv(\\.gz)?$",
                          expectedRet = c(8, 5))
    if (length(searchFiles) > 0) {
        .abundancePlot(
            searchFiles,
            # what are the sample names (in-order)
            sampleNames,
            # output directory
            abunOut
        )
    }

    # plot igv mismatches distribution
    abunIgvMismatchFiles <-
        .listFilesInOrder(path = abundanceDirectories,
                          pattern = ".*_igv_mismatches_dist\\.csv(\\.gz)?$")

    if (length(abunIgvMismatchFiles) > 0) {
        subtitle <- paste("Total is",
                          paste(lapply(abunIgvMismatchFiles, function(x) {
                              as.integer(.getTotal(x))
                          }), collapse = ", "))

        abunIgvMismatches <- .plotDist(
            lapply(abunIgvMismatchFiles, read.csv, skip = 1),
            sampleNames,
            paste("Number of mismatches in V gene in", combinedNames),
            .checkVert(abunIgvMismatchFiles[[1]]),
            subs = subtitle
        )
        ggsave(paste0(abunOut, mashedNames, "_igv_mismatches_dist.png"),
               plot = abunIgvMismatches, width = V_WIDTH, height = V_HEIGHT)
    }

    # plot igv gaps distribution
    abunIgvGapsFiles <-
        .listFilesInOrder(path = abundanceDirectories,
                          pattern = ".*_igv_gaps_dist\\.csv(\\.gz)?$")
    if (length(abunIgvGapsFiles) > 0) {
        subtitle <- paste("Total is",
                          paste(lapply(abunIgvGapsFiles, function(x) {
                              as.integer(.getTotal(x))
                          }), collapse = ", "))
        abunIgvGaps <- .plotDist(
            lapply(abunIgvGapsFiles, read.csv, skip = 1),
            sampleNames,
            paste("Number of gaps in V gene in ", combinedNames),
            .checkVert(abunIgvGapsFiles[[1]]),
            subs = subtitle
        )
        ggsave(paste0(abunOut, mashedNames, "_igv_gaps_dist.png"),
               plot = abunIgvGaps, width = V_WIDTH, height = V_HEIGHT)
    }

    if (length(sampleNames) == 1) {
        # we can plot circlize if there's only one sample
        .plotCirclize(sampleNames[1], abunOut)
    }
}
