#' Abundance distribution
#'
#' @import ggplot2
#' @include util.R
#' @include distributions.R
#'
#' @param files list type. list of files in abundance directory
#' @param sampleNames vector type. 1-1 correspondance to files
#' @param outputDir string type.
#' @param skipDgene logical type. Skip D germline abundance plot if TRUE.
#' @param .save logical type. Save Rdata ggplot item
#'
#' @return None
.abundancePlot <- function(files, sampleNames, outputDir,
                           skipDgene = FALSE, .save = TRUE) {

    vdj <- if (skipDgene) c("v", "j") else c("v", "d", "j")

    for (expression in c("family", "gene")) {
        for (gene in vdj) {
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
                filename <- file.path(outputDir, paste0(paste(sampleNames, collapse = "_"), "_ig", gene, "_dist_", expression, "_level.png"))
                ggsave(filename, plot = p, width = width, height = height)
                .saveAs(.save, filename, plot = p)
            }
        }
    }
}

#' V-J association plot
#'
#' @import circlize
#' @import RColorBrewer
#'
#' @param sampleName string type
#' @param path string type. Path to _vjassoc.csv
#' @param outputdir string type
#'
#' @return None
.plotCirclize <- function(sampleName, path, outputdir) {
    filename <- file.path(path, paste0(sampleName, "_vjassoc.csv"))

    message(paste("Plotting V-J association for", sampleName))

    if (file.exists(filename)) {
        df <- read.csv(filename)

        # output file
        outputFileName <- file.path(outputdir, gsub(".csv",
                                                    ".png",
                                                    basename(filename),
                                                    fixed = T))

        png(outputFileName, width = V_WIDTH, height = V_HEIGHT,
            units = "in", res = 1200, pointsize = 10)

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
        title(sampleName, cex = 8)
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
    } else {
        warning(paste("Could not find file", filename,
                      "skipping V-J association plot"))
    }
}


#' Conducts abundance analysis
#'
#' @import ggplot2
#' @include util.R
#' @include distributions.R
#'
#' @param abundanceDirectories list type. List of sample directories
#' @param abunOut string type. Output directory
#' @param sampleNames vector type. 1-1 correspondence with abundanceDirectories
#' @param combinedNames string type. Title "combined" sample names
#' @param mashedNames string type. File "mashed" names - avoid special chars
#' @param skipDgene logical type. Skip D gene plots?
#' @param .save logical type. Save ggplot as Rdata
#'
#' @return None
.abundanceAnalysis <- function(abundanceDirectories, abunOut,
                               sampleNames, combinedNames, mashedNames,
                               skipDgene = FALSE, .save = TRUE) {
    # where to find the files
    # 3 each from V and D, then 2 from J (no gene) or 5 (exclude the 3 from D)
    searchFiles <-
        .listFilesInOrder(path = abundanceDirectories,
                          pattern = ".*ig[vdj]_dist_[family|gene|variant].*\\.csv(\\.gz)?$",
                          expectedRet = c(8, 5))
    if (length(searchFiles) > 0) {
        message(paste("Plotting V(D)J abundance distribution for",
                      combinedNames))
        .abundancePlot(
            searchFiles,
            # what are the sample names (in-order)
            sampleNames,
            # output directory
            abunOut,
            # whether or not to ignore D gene plots
            skipDgene = skipDgene
        )
    } else {
        warning(paste("Could not find V(D)J abundance CSV files in",
                      paste(abundanceDirectories, collapse = ",")))
    }

    # plot igv mismatches distribution
    abunIgvMismatchFiles <-
        .listFilesInOrder(path = abundanceDirectories,
                          pattern = ".*_igv_mismatches_dist\\.csv(\\.gz)?$")

    if (length(abunIgvMismatchFiles) > 0) {
        message(paste("Plotting IGV mismatches distribution for", combinedNames))
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
        saveName <- file.path(abunOut, paste0(mashedNames, "_igv_mismatches_dist.png"))
        ggsave(saveName, plot = abunIgvMismatches, width = V_WIDTH, height = V_HEIGHT)
        .saveAs(.save, saveName, plot = abunIgvMismatches)
    } else {
        warning("Could not find IGV mismatches distribution CSV files in",
                paste(abundanceDirectories, collapse = ","))
    }

    # plot igv gaps distribution
    abunIgvGapsFiles <-
        .listFilesInOrder(path = abundanceDirectories,
                          pattern = ".*_igv_gaps_dist\\.csv(\\.gz)?$")
    if (length(abunIgvGapsFiles) > 0) {
        message(paste("Plotting IGV indels distribution for", combinedNames))
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
        saveName <- file.path(abunOut, paste0(mashedNames, "_igv_gaps_dist.png"))
        ggsave(saveName, plot = abunIgvGaps, width = V_WIDTH, height = V_HEIGHT)
        .saveAs(.save, saveName, plot = abunIgvGaps)
    } else {
        warning("Could not find IGV indels distribution CSV files in",
                paste(abundanceDirectories, collapse = ","))
    }

    if (length(sampleNames) == 1) {
        # we can plot circlize if there's only one sample
        .plotCirclize(sampleNames[1], abundanceDirectories[[1]], abunOut)
    }
}
