#' Returns a list of nC2 pairwings
#'
#' @import BiocParallel
#'
#' @param dataframes
#'
#' @return
.pairwiseDataFrames <- function(dataframes) {
    # TODO: make this parallel
    # nC2 operation
    sampleNames <- names(dataframes)
    dfs <-  lapply(head(seq_along(dataframes), n = -1), function(i) {
        pairs <- lapply((i + 1):length(dataframes), function(j) {
            df1 <- dataframes[[i]]
            df2 <- dataframes[[j]]

            df.union <- merge(df1, df2, by = "Clonotype", all.y = T, all.x = T)
        })
        names(pairs) <- unlist(lapply((i + 1):length(dataframes), function(j) {
            paste0(sampleNames[i], "_vs_", sampleNames[j])
        }))
        return(pairs)
    })
    # to get intersecting clones, it's as easy as finding rows with
    # no NAs:
    # intersect.clones <- df.union[complete.cases(df.union), "Clonotype"]
    return(unlist(dfs, recursive = F))
}

#' Title Creates a scatter plot
#'
#' @import ggplot2
#' @include util.R
#'
#' @param df1 dataframe for sample 1
#' @param df2 dataframe for sample 2
#' @param name1 string type, Sample 1 name
#' @param name2 string type. Sample 2 name
#' @param cloneClass string type.
#' What region was used to classify clonotypes - appears in title. For example,
#' CDR3 or V region
#'
#' @return ggplot2 object
.scatterPlot <- function(df1, df2, name1, name2, cloneClass) {

    message(paste("Generating scatter plot for", name1, "and", name2))

    df.union <- merge(df1, df2, by = "Clonotype", all.y = TRUE, all.x = TRUE)

    # replace NaN with 0
    df.union[is.na(df.union)] <- 0

    # plot!
    gg <- ggplot(df.union, aes(x = Count.x, y = Count.y)) +
        geom_point() +
        theme_bw() +
        labs(subtitle = paste(name2, " vs ", name1,
                              " plot based on clonotype counts"),
             y = name2,
             x = name1,
             title = paste("Scatter plot of ", cloneClass, " clonotype counts"))
    return(gg)
}

#' Creates a complex scatter plot
#'
#' @import ggplot2
#' @import gridExtra grid
#' @include util.R
#'
#' @param df.union a 'lossless' dataframe created by intersecting sample1 and
#' sample2's dataframes. It should contain NAs where clones that appear in one
#' sample doesn't appear in the other. For example:
#'
#' +-------------------------------------------------+
#' | Clonotype | prop.x | prop.y | Count.x | Count.y |
#' +-------------------------------------------------+
#' | ABCDEF       NA       0.01      NA        210   |
#' | ......                                          |
#' +-------------------------------------------------+
#'
#' @param df1 dataframe for sample 1
#' @param df2 dataframe for sample 2
#' @param name1 string type, Sample 1 name
#' @param name2 string type. Sample 2 name
#' @param cloneClass string type.
#' What region was used to classify clonotypes - appears in title. For example,
#' CDR3 or V region
#'
#' this plotting techique was
#' shamelessly plagarised from
#' https://github.com/mikessh/vdjtools/blob/master/src/main/resources/rscripts/intersect_pair_scatter.r
#' (VDJTools) with minor modifications
#'
#' @return ggplot2 object
.scatterPlotComplex <- function(df.union, df1, df2, name1, name2, cloneClass) {
    message(paste("Generating scatter plot for", name1, "and", name2))

    # find what the smallest percentage is so we know what to use as the
    # minimum (since log10(0) wont cut it)
    smallestPercentage <- min(df1$prop, df2$prop)

    # log10 scale the proportions (-inf, 0)
    df.union$prop.x <- log10(df.union$prop.x)
    df.union$prop.y <- log10(df.union$prop.y)
    df1$prop <- log10(df1$prop)
    df2$prop <- log10(df2$prop)


    intersectingClones <- df.union[complete.cases(df.union), "Clonotype"]

    # convert NAs to x and y intercept values
    df.union[is.na(df.union)] <- log10(smallestPercentage * 5e-1)

    xmin <- min(df.union$prop.x)
    ymin <- min(df.union$prop.y)

    sample1.margin <- .cloneDistMarginal(df1, intersectingClones, xmin, flip = F)
    sample2.margin <- .cloneDistMarginal(df2, intersectingClones, ymin, flip = T)

    g <- ggplot(df.union, aes(x = prop.x, y = prop.y)) +
        theme_bw() +
        geom_point(aes(size = (prop.x + prop.y / 2)),
                   color = "black",
                   alpha = 0.4,
                   pch = 21,
                   fill = BLUEHEX) +
        scale_size_continuous(guide = "none", range = c(1, 10)) +
        scale_x_continuous(limits = c(xmin, 0), expand = c(0, 0.1)) +
        scale_y_continuous(limits = c(ymin, 0), expand = c(0, 0.1)) +
        # geom_smooth(method = "lm", se = T, fullrange = T) +
        # stat_smooth(data = df.union, aes(prop.x, prop.y,
        #                                  weight = 10^((prop.x + prop.y) / 2)),
        #             color = "blue", method = "lm", fullrange = T, se = T) +
        labs(y = name2, x = name1)
    grid.arrange(sample1.margin, .emptyPlot(), g, sample2.margin,
                 ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4),
                 top = paste("Scatter plot of", cloneClass, "clonotype frequencies\n", name2, "vs", name1)
    )
}



#' Marginal density graph of clonotypes (blue for shared, grey for total)
#'
#' @import ggplot2
#' @include util.R
#'
#' @param df.original
#' @param otherClones
#' @param lim.min
#' @param flip
#'
#' @return ggplot2 object
.cloneDistMarginal <- function(df.original, otherClones, lim.min, flip) {
    mask <- df.original$Clonotype %in% otherClones
    # clones that are shared between df.original and the other sample
    df.shared <- df.original[mask, ]
    # clones that are not shared with the other sample
    df.exclusive <- df.original[!mask, ]
    g <- ggplot() +
        stat_density(data = df.original, aes(x = prop, y = ..scaled..),
                     fill = "#808080", # grey
                     alpha = 0.4, adjust = 1, size = 0.1, color = "black") +
        stat_density(data = df.shared, aes(x = prop, y = ..scaled..),
                     fill = BLUEHEX,
                     alpha = 0.4, adjust = 1, size = 0.1, color = "black") +
        stat_density(data = df.exclusive, aes(x = prop, y = ..scaled..),
                     fill = "#e0ccff", # purplish blue
                     alpha = 0.4, adjust = 1, size = 0.1, color = "black") +
        scale_x_continuous(limits = c(lim.min, 0), expand = c(0, 0.25)) +
        theme_bw() +
        theme(legend.position = "none", axis.title.x = element_blank(),
              axis.text = element_blank(), axis.ticks = element_blank(),
              axis.title.y = element_blank(), panel.grid = element_blank())
    if (flip) {
        g <- g + coord_flip()
    }
    return(g)
}

#' Marginal histogram of clonotypes (blue for shared, grey for total). The y
#' axis is scaled by sqrt (but it doesn't really matter anyway, since we're
#' stripping away the y-ticks)
#'
#' @import ggplot2
#' @include util.R
#'
#' @param df.original
#' @param df.filtered
#' @param lim.min
#' @param flip
#'
#' @return ggplot2 object
.cloneDistHist <- function(df.original, otherClones, lim.min, flip) {
    df.filtered <- df.original[df.original$Clonotype %in% otherClones, ]
    g <- ggplot() +
        geom_histogram(data = df.original, aes(x = prop,
                                               y = sqrt(..count../sum(..count..))),
                       fill = "#808080", alpha = 0.4) +
        geom_histogram(data = df.filtered, aes(x = prop,
                                               y = sqrt(..count../sum(..count..))),
                       fill = BLUEHEX, alpha = 0.4) +
        scale_x_continuous(limits = c(lim.min, 0), expand = c(0, 0.25)) +
        theme_bw() +
        theme(legend.position = "none", axis.title.x = element_blank(),
              axis.text = element_blank(), axis.ticks = element_blank(),
              axis.title.y = element_blank(), panel.grid = element_blank())
    if (flip) {
        g <- g + coord_flip()
    }
    return(g)
}

#' Conduct all pairwise comparison analyses
#'
#' @include util.R
#'
#' @import ggplot2
#' @import BiocParallel
#'
#' @param dataframes
#' @param sampleNames
#' @param outputPath
#' @param save
#'
#' @return
.pairwiseComparison <- function(dataframes, sampleNames, outputPath, .save = T) {
    stopifnot(length(dataframes) == length(sampleNames))
    names(dataframes) <- sampleNames
    pairs <- .pairwiseDataFrames(dataframes)
    lapply(seq_along(pairs), function(i) {
        comparisonName <- unlist(strsplit(names(pairs)[i], split = "_vs_", fixed = T))

        # -----------------------  scatter plot ---------------------------
        p <- .scatterPlotComplex(pairs[[i]],
                                 dataframes[[comparisonName[1]]],
                                 dataframes[[comparisonName[2]]],
                                 comparisonName[1],
                                 comparisonName[2],
                                 "CDR3")
        saveName <- file.path(outputPath, paste0(names(pairs)[i], "_clone_scatter.png"))
        # square plot to get a evenly scaled scatter plot on both the x and y axis
        ggsave(saveName, plot = p, width = V_WIDTH, height = V_WIDTH)
        # too large to save, skip this!
        # .saveAs(.save, saveName, p)

        # -----------------------  something else ---------------------------
        message("Going to add more analyses here")
    })
}



#' Comprehensive clonotype analyses
#'
#' @include util.R
#'
#' @param diversityDirectories list type. List of directories to diversity dir
#' @param clonotypeOut string type. Output directory
#' @param sampleNames vector type. 1-1 with diversityDirectories
#' @param mashedNames string type. Prefix for ooutput files using "mashed-up"
#' @param .save logical type. Save ggplot object?
#'
#' @return Nothing
.clonotypeAnalysis <- function(diversityDirectories,
                              clonotypeOut, sampleNames,
                              mashedNames, .save = T) {
    message(paste("Starting clonotype analysis on samples",
                  paste(sampleNames, collapse = ",")))
    cdr3ClonesFile <- .listFilesInOrder(path = diversityDirectories,
                                        pattern = ".*_cdr3_clonotypes_.*_over\\.csv(\\.gz)?$")
    if (length(cdr3ClonesFile) > 1) {
        dataframes <- lapply(cdr3ClonesFile, function(fname) {
            df <- read.csv(fname, stringsAsFactors = F)
            df$prop <- df$Count / sum(df$Count)
            return(df[, c("Count", "Clonotype", "prop")])
        })
        # all vs all comparisons (scatter plot, similarity, etc ...)
        .pairwiseComparison(dataframes, sampleNames, clonotypeOut, .save = .save)

        # venn diagram
        .vennIntersection(dataframes, sampleNames,
                          file.path(clonotypeOut,
                                    paste0(mashedNames, "_cdr3_clonotypeIntersection.png")))
        # top 10 barplot
        g <- .topNDist(dataframes, sampleNames)
        saveName <- file.path(clonotypeOut,
                              paste0(mashedNames, "_top10Clonotypes.png"))
        ggsave(saveName, plot = g, width = V_WIDTH_L, height = V_HEIGHT_L)
        .saveAs(.save, saveName, g)
    }
}
