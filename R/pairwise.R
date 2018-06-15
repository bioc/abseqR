#' Title Clonotype table
#'
#' @import ggplot2
#' @import RColorBrewer
#'
#' @param dataframes list type. List of dataframes.
#' @param sampleNames vector type. vector of strings representing sample
#' names should have one-to-one correspondence with dataframes
#' @param top int type. Top N clonotypes to plot
#'
#' @return None
.topNDist <- function(dataframes, sampleNames, top = 10) {
    nsamples <- length(dataframes)

    if (nsamples != length(sampleNames)) {
        stop(paste("Expected equal number of sample names and dataframes, got",
                   length(sampleNames), "samples and", nsamples, "dataframes."))
    }

    message(paste("Plotting top", top, "clonotype distribution for samples",
                  paste(sampleNames, collapse = ", ")))

    # --- cleanup & pre-processing ---
    colNames <- c("Clonotype", "Count")
    for (i in seq_len(nsamples)) {
        df <- dataframes[[i]]
        # only need top N, also take the columns we're interested in only
        df <- head(df[colNames], top)
        # append sample name to distinguish data when merged later on
        df$sample <- rep(sampleNames[i], nrow(df))
        # normalize percentage to top N
        df$normPerc <- df$Count / sum(df$Count)
        dataframes[[i]] <- df
    }

    df.union <- do.call("rbind", dataframes)
    # --- done: cleanup & pre-processing ---

    # plot!
    # colour suggestion taken from
    # https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes/41230685
    # and modified
    palette <- c("dodgerblue2","#E31A1C", # red
                 "green4",
                 "#6A3D9A", # purple
                 "#FF7F00", # orange
                 "black","gold1",
                 "skyblue2","#FB9A99", # lt pink
                 "olivedrab1",
                 "#CAB2D6", # lt purple
                 "#FDBF6F", # lt orange
                 "gray70", "khaki2",
                 "maroon","orchid1","deeppink1","blue1","steelblue4",
                 "darkturquoise","green1","yellow4","yellow3",
                 "aquamarine", "darkorange4", "mediumpurple1", "dimgrey",
                 "darkseagreen1", "lightyellow", "coral2")
    if (length(unique(df.union$Clonotype)) > 30) {
        warning(paste0("Too many unique clonotypes are being plotted",
                       " - extrapolatingpalette for top10 clonotype dist plot."))
        getPalatte <- colorRampPalette(brewer.pal(8, 'Accent'))
        palette <- getPalatte(length(unique(df.union$Clonotype)))
    }
    g <- ggplot(df.union, aes(x = sample, y = normPerc)) +
        geom_bar(stat ='identity', aes(fill = Clonotype)) +
        theme(legend.position = "bottom", legend.box = "horizontal",
              legend.title = element_blank(),
              legend.text = element_text(size = 7)) +
        labs(title = paste("Top", top, "clonotype across each sample"),
             subtitle = paste("Colour coded clonotypes, distribution of each clonotype is scaled to top", top),
             x = "Sample",
             y = "Distribution")  +
        scale_fill_manual(values = palette)
    return(g)
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



#' Title Creates Venndiagram for clonotype intersection
#'
#' @import VennDiagram
#'
#' @param dataframes list type. List of sample dataframes. Only accepts 2 - 5
#' samples. Warning message will be generated for anything outside of the range
#' @param sampleNames vector type. 1-1 with dataframes
#' @param outFile string type. Filename to be saved as
#' @param top int type. Top N cutoff, defaults to ALL clones if not specified
#'
#' @return Nothing
.vennIntersection <- function(dataframes, sampleNames, outFile, top = Inf) {

    nsample <- length(dataframes)

    if (nsample != length(sampleNames)) {
        stop(paste("Expected equal number of sample names and dataframes, got",
                   length(sampleNames), "samples and", nsample, "dataframes."))
    }

    if (nsample >= 2 && nsample <= 5) {
        # output
        message(paste("Creating Venn diagram for samples",
                      paste(sampleNames, collapse = ", ")))
        png(file = outFile, width = 8, height = 7, units = "in", res = 300)

        # Get the top N clonotypes if specified, and only use the 2 columns
        # specified below
        colNames <- c("Clonotype", "Count")
        dataframes <- lapply(dataframes, function(df) {
            head(df[colNames], top)
        })

        # merge all dataframes
        df.union <- merge(dataframes[[1]], dataframes[[2]],
                          by = "Clonotype", all.y = TRUE, all.x = TRUE)
        colnames(df.union) <- c("Clonotype", sampleNames[1], sampleNames[2])
        df.union[is.na(df.union)] <- 0

        # check if there's more
        if (nsample > 2) {
            for (i in 3:nsample) {
                oldColNames <- colnames(df.union)
                df.union <- merge(df.union, dataframes[[i]], by = "Clonotype",
                                  all.y = TRUE, all.x = TRUE)
                colnames(df.union) <- c(oldColNames, sampleNames[i])
                df.union[is.na(df.union)] <- 0
            }
        }

        # plot!
        area1 <- sum(df.union[sampleNames[1]] > 0)
        area2 <- sum(df.union[sampleNames[2]] > 0)
        area12 <- sum(df.union[sampleNames[1]] > 0 & df.union[sampleNames[2]] > 0)
        if (nsample == 2) {
            draw.pairwise.venn(area1, area2, area12, category = c(sampleNames[1], sampleNames[2]),
                               lty = "blank",
                               col = "transparent",
                               cat.fontface = "bold",
                               cex = 1.7,
                               fill = 2:3,
                               alpha = 0.5,
                               filename = NULL)
        } else {
            area3 <- sum(df.union[sampleNames[3]] > 0)
            area13 <- sum(df.union[sampleNames[1]] > 0 & df.union[sampleNames[3]] > 0)
            area23 <- sum(df.union[sampleNames[2]] > 0 & df.union[sampleNames[3]] > 0)
            area123 <- sum(df.union[sampleNames[1]] > 0 & df.union[sampleNames[2]] > 0 & df.union[sampleNames[3]] > 0)
            if (nsample == 3) {
                draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, n12 = area12, n13 = area13, n23 = area23, n123 = area123,
                                 category = c(sampleNames[1], sampleNames[2], sampleNames[3]),
                                 lty = "blank",
                                 col = "transparent", cat.fontface = "bold", cex = 1.7,
                                 fill = 2:4, alpha = 0.5, filename = NULL)
            } else {
                area4 <- sum(df.union[sampleNames[4]] > 0)
                area14 <- sum(df.union[sampleNames[1]] > 0 & df.union[sampleNames[4]] > 0)
                area24 <- sum(df.union[sampleNames[2]] > 0 & df.union[sampleNames[4]] > 0)
                area34 <- sum(df.union[sampleNames[3]] > 0 & df.union[sampleNames[4]] > 0)
                area124 <- sum(df.union[sampleNames[1]] > 0 & df.union[sampleNames[2]] > 0 & df.union[sampleNames[4]] > 0)
                area134 <- sum(df.union[sampleNames[1]] > 0 & df.union[sampleNames[3]] > 0 & df.union[sampleNames[4]] > 0)
                area234 <- sum(df.union[sampleNames[2]] > 0 & df.union[sampleNames[3]] > 0 & df.union[sampleNames[4]] > 0)
                area1234 <- sum(df.union[sampleNames[1]] > 0 &
                                    df.union[sampleNames[2]] > 0 &
                                    df.union[sampleNames[3]] > 0 &
                                    df.union[sampleNames[4]] > 0)
                if (nsample == 4) {
                    draw.quad.venn(area1 = area1, area2 = area2, area3 = area3, area4 = area4,
                                   n12 = area12, n13 = area13, n14 = area14, n23 = area23, n24 = area24, n34 = area34,
                                   n123 = area123, n124 = area124, n134 = area134, n234 = area234, n1234 = area1234,
                                   category = c(sampleNames[1], sampleNames[2], sampleNames[3], sampleNames[4]),
                                   lty = "blank",
                                   col = "transparent", cat.fontface = "bold", cex = 1.7,
                                   fill = 2:5, alpha = 0.5, filename = NULL
                    )
                } else {
                    # quuntuple plot
                    area5 <- sum(df.union[sampleNames[5]] > 0)
                    area15 <- sum(df.union[sampleNames[1]] > 0 & df.union[sampleNames[5]] > 0)
                    area25 <- sum(df.union[sampleNames[2]] > 0 & df.union[sampleNames[5]] > 0)
                    area35 <- sum(df.union[sampleNames[3]] > 0 & df.union[sampleNames[5]] > 0)
                    area45 <- sum(df.union[sampleNames[4]] > 0 & df.union[sampleNames[5]] > 0)
                    area125 <- sum(df.union[sampleNames[1]] > 0 & df.union[sampleNames[2]] > 0 & df.union[sampleNames[5]] > 0)
                    area135 <- sum(df.union[sampleNames[1]] > 0 & df.union[sampleNames[3]] > 0 & df.union[sampleNames[5]] > 0)
                    area145 <- sum(df.union[sampleNames[1]] > 0 & df.union[sampleNames[4]] > 0 & df.union[sampleNames[5]] > 0)
                    area235 <- sum(df.union[sampleNames[2]] > 0 & df.union[sampleNames[3]] > 0 & df.union[sampleNames[5]] > 0)
                    area245 <- sum(df.union[sampleNames[2]] > 0 & df.union[sampleNames[4]] > 0 & df.union[sampleNames[5]] > 0)
                    area345 <- sum(df.union[sampleNames[3]] > 0 & df.union[sampleNames[4]] > 0 & df.union[sampleNames[5]] > 0)
                    area1235 <- sum(df.union[sampleNames[1]] > 0 &
                                        df.union[sampleNames[2]] > 0 &
                                        df.union[sampleNames[3]] > 0 &
                                        df.union[sampleNames[5]] > 0)
                    area1245 <- sum(df.union[sampleNames[1]] > 0 &
                                        df.union[sampleNames[2]] > 0 &
                                        df.union[sampleNames[4]] > 0 &
                                        df.union[sampleNames[5]] > 0)
                    area1345 <- sum(df.union[sampleNames[1]] > 0 &
                                        df.union[sampleNames[3]] > 0 &
                                        df.union[sampleNames[4]] > 0 &
                                        df.union[sampleNames[5]] > 0)
                    area2345 <- sum(df.union[sampleNames[2]] > 0 &
                                        df.union[sampleNames[3]] > 0 &
                                        df.union[sampleNames[4]] > 0 &
                                        df.union[sampleNames[5]] > 0)
                    area12345 <- sum(df.union[sampleNames[1]] > 0 &
                                         df.union[sampleNames[2]] > 0 &
                                         df.union[sampleNames[3]] > 0 &
                                         df.union[sampleNames[4]] > 0 &
                                         df.union[sampleNames[5]] > 0)
                    draw.quintuple.venn(
                        area1 = area1, area2 = area2, area3 = area3, area4 = area4, area5 = area5,
                        n12 = area12, n13 = area13, n14 = area14, n15 = area15, n23 = area23, n24 = area24, n25 = area25, n34 = area34,
                        n35 = area35, n45 = area45, n123 = area123, n124 = area124, n125 = area125, n134 = area134,
                        n135 = area135, n145 = area145, n234 = area234, n235 = area235, n245 = area245, n345 = area345, n1234 = area1234,
                        n1235 = area1235, n1245 = area1245, n1345 = area1345, n2345 = area2345, n12345 = area12345,
                        category = c(sampleNames[1], sampleNames[2], sampleNames[3], sampleNames[4], sampleNames[5]),
                        lty = "blank",
                        col = "transparent", cat.fontface = "bold", cex = 1.7,
                        fill = 2:6, alpha = 0.5, filename = NULL
                    )
                }
            }
        }
        dev.off()
    } else {
        warning(paste("Skipping venn diagram plot for",
                      paste(sampleNames, collapse = ", "),
                      "because they do not fall within the range of 2 <= x <= 5"))
    }
}


#' Conduct all pairwise comparison analyses
#'
#' @include util.R
#' @include statistics.R
#'
#' @import ggplot2
#' @import BiocParallel
#' @import ggcorrplot
#' @import reshape2
#'
#' @param dataframes
#' @param sampleNames
#' @param outputPath
#' @param .save
#'
#' @return nothing
.pairwiseComparison <- function(dataframes, sampleNames, outputPath, .save = TRUE) {

    stopifnot(length(dataframes) == length(sampleNames))
    # TODO: depending on memory consumption and CPU usage, can parallelize this.
    df <- do.call("rbind", lapply(head(seq_along(dataframes), n = -1), function(i) {
        do.call("rbind", lapply((i + 1):length(dataframes), function(j) {

            df.union <- merge(dataframes[[i]],
                              dataframes[[j]],
                              by = "Clonotype",
                              all.y = T,
                              all.x = T)

            # --- scatter plot ---
            p <- .scatterPlotComplex(df.union,
                                     dataframes[[i]],
                                     dataframes[[j]],
                                     sampleNames[i],
                                     sampleNames[j],
                                     "CDR3")
            saveName <- file.path(outputPath,
                                  paste0(sampleNames[i], "_vs_",
                                         sampleNames[j], "_clone_scatter.png"))
            # square plot to get a evenly scaled scatter plot on both the x and y axis
            ggsave(saveName, plot = p, width = V_WIDTH, height = V_WIDTH)
            # too large to save, skip this!
            # .saveAs(.save, saveName, p)
            # -- end scatter plot --

            # make all NAs 0 before conducting distance / similarity analyses
            df.union <- replace(df.union, is.na(df.union), 0)
            # --- pearson correlation of clonotype frequencies ---
            correlations <- suppressWarnings(.correlationTest(df.union))
            # -- end pearson correaltion --

            # --- Morisita Horn similarity index ---
            distance <- .distanceMeasure(df.union)
            # -- end Morisita Horn similarity index --


            data.frame(from = sampleNames[i],
                       to = sampleNames[j],
                       morisita.horn = distance$morisita.horn,
                       jaccard = distance$jaccard,
                       bray.curtis = distance$bray.curtis,
                       pearson = correlations$pearson,
                       pearson.p = correlations$pearson.p,
                       spearman = correlations$spearman,
                       spearman.p = correlations$spearman.p
                       )
        }))
    }))
    # what df should loop like now: length(sampleNames) Choose 2 rows
    # +------------------------------------------------------------------------------------------------+
    # |from | to | morisita.horn | jaccard | bray.curtis | pearson | pearson.p | spearman | spearman.p |
    # +------------------------------------------------------------------------------------------------+
    # | ...                                                                                            |
    # +------------------------------------------------------------------------------------------------+
    #acast(df, from ~ to, value.var = "a", fill = 0, fun.aggregate = sum)
    # step 1. save the dataframe as a TSV file
    write.table(file = file.path(outputPath, "indices.tsv"), df,
                quote = FALSE, row.names = FALSE, sep = "\t")

    # step 2. plot and save the correlograms
    lapply(c("pearson", "spearman"), function(method) {
        mat <- acast(df, from ~ to, value.var = method)
        p <- ggcorrplot(mat,
                        lab = TRUE,
                        ggtheme = theme_bw,
                        method = "circle",
                        title = paste(.capitalize(method), "correlation"))
        saveName <- file.path(outputPath, paste0(method, ".png"))
        ggsave(saveName, plot = p, width = V_WIDTH, height = V_HEIGHT)
        .saveAs(.save, saveName, p)
    })

}


#' Comprehensive clonotype analyses
#'
#' @include util.R
#'
#' @import ggplot2
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

