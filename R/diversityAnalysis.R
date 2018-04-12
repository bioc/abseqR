#' Title Creates a scatter plot
#'
#' @import ggplot2
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
    options(scipen = 999)
    theme_set(theme_bw())

    message(paste("Generating scatter plot for", name1, "and", name2))

    df.union <- merge(df1, df2, by = "Clonotype", all.y = TRUE, all.x = TRUE)

    # replace NaN with 0
    df.union[is.na(df.union)] <- 0

    # plot!
    gg <- ggplot(df.union, aes(x = Count.x, y = Count.y)) +
        geom_point() +
        labs(subtitle = paste(name2, " vs ", name1,
                            " plot based on clonotype counts"),
             y = name2,
             x = name1,
             title = paste("Scatter plot of ", cloneClass, " clonotype counts"))
    return(gg)
}


#' Title Creates pairing scatter plots
#'
#' Plots a scatterplot of (i, i+1) samples. EG: given samples 1-3, plots
#' a scatterplot of pairs: (1,2) and (2,3).
#'
#' @include util.R
#' @import ggplot2
#'
#' @param dataframes list type. List of dataframes
#' @param sampleNames vector type. 1-1 with dataframes
#' @param outputPath string type. Output directory
#' @param cloneClass string type.
#' What region was used to classify clonotypes - appears in title. For example,
#' CDR3 or V region
#' @param .save logical type. Save ggplot object?
#'
#' @return None
.scatterClones <- function(dataframes, sampleNames, outputPath, cloneClass,
                           .save = T) {
    nsamples <- length(dataframes)
    message(paste("Plotting pairwise scatter plot for",
                  paste(sampleNames, collapse = ", ")))
    # this scatter plot doesn't make sense for 1 sample only
    if (nsamples <= 1) {
        stop("Expected > 1 samples for scatter plot")
    }
    if (nsamples != length(sampleNames)) {
        stop(paste("Expected equal number of sample names and dataframes, got",
                   length(sampleNames), "samples and", nsamples, "dataframes."))
    }

    # only need these 2 columns from dataframe
    colnames <- c("Clonotype", "Count")
    for (i in seq_len(nsamples - 1)) {
        p <- .scatterPlot(dataframes[[i]][colnames],
                          dataframes[[i + 1]][colnames],
                          sampleNames[i],
                          sampleNames[i + 1], cloneClass)
        saveName <- file.path(outputPath, paste0(sampleNames[i], "_vs_",
                      sampleNames[i + 1], "_clone_scatter.png"))
        ggsave(saveName, plot = p, width = V_WIDTH, height = V_HEIGHT)
        .saveAs(.save, saveName, p)
    }
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



#' Title
#'
#' @import ggplot2
#' @include util.R
#'
#' @param files list type. List of strings to _cdr_v_duplication.csv pathname
#' @param sampleNames vector type. Vector of strings each representing
#' sample names
#' @param regions vector type.
#'  Which regions to include in the plot. Default = c("CDR3", "V")
#'
#' @return ggplot2 object
.plotDuplication <- function(files, sampleNames, regions = c("CDR3", "V")) {
    theme_set(theme_bw())
    nsamples <- length(files)

    if (nsamples != length(sampleNames)) {
        stop(paste("Expected equal number of sample names and dataframes, got",
                   length(sampleNames), "samples and", nsamples, "dataframes."))
    }

    message(paste("Creating duplication plot for samples",
                  paste(sampleNames, collapse = ", ")))

    # read xticks and xlimits from first 2 row
    trimwsNoQuotes <- function(x) {
        gsub("'", "", trimws(x))
    }

    # trimwsNoQuotes strips whitespaces AND single quotes
    fp <- file(files[[1]], "r")
    xticks <- strtoi(unlist(lapply(strsplit(readLines(fp, n = 1), ",")[[1]],
                                   trimws)))
    xlabels <- unlist(lapply(strsplit(readLines(fp, n = 1), ",")[[1]],
                             trimwsNoQuotes))
    close(fp)

    # read files into dataframes
    dataframes <- lapply(files, read.csv, skip = 2)

    # pre-processing & cleanup
    for (i in seq_len(nsamples)) {
        df <- dataframes[[i]]
        df$sample  <- rep(sampleNames[[i]], nrow(df))
        dataframes[[i]] <- df[df$region %in% regions, ]
    }

    # combine!
    df.union <- do.call("rbind", dataframes)

    g <- ggplot(df.union, aes(x = x, y = y))

    if (nsamples == 1) {
        g <- g + geom_line(aes(linetype = region, color = sample),
                           color = BLUEHEX, size = 0.65) +
            guides(color = FALSE) +
            #scale_size_manual(values = head(seq(1, length(regions), by = .5),
            #                                n = length(regions))) +
            scale_linetype_manual(values = .getLineTypes(regions))
    } else {
        g <- g + geom_line(aes(linetype = region, color = sample),
                           size = 0.65) +
            #scale_size_manual(values = head(seq(1, length(regions), by = .5),
            #                                n = length(regions))) +
            scale_linetype_manual(values = .getLineTypes(regions))
    }

    g <- g + scale_x_continuous(breaks = xticks, labels = xlabels) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = paste("Sequence duplication levels of",
                           paste(regions, collapse = ", "),
                           "in", paste(sampleNames, collapse = ", ")),
             x = "Duplication level",
             y = "Proportion of duplicated sequences")
    return(g)
}




#' Title Rarefaction plot
#'
#' @import ggplot2
#' @include util.R
#'
#' @param files list type. A list of files consisting of path to samples
#' @param sampleNames vector type. A vector of strings,
#'  each being the name of samples in files
#' @param regions vector type. A vector of strings,
#' regions to be included. Defaults to c("CDR3", "V")
#'
#' @return ggplot2 object
.plotRarefaction <- function(files, sampleNames, regions = c("CDR3", "V")) {
    theme_set(theme_bw())

    nsamples <- length(files)
    # sanity check
    if (length(sampleNames) != nsamples) {
        stop(paste("Expected equal number of sample names and dataframes, got",
                   length(sampleNames), "samples and", nsamples, "dataframes."))
    }

    message(paste("Creating rarefaction plot for samples",
            paste(sampleNames, collapse = ", ")))

    # find the minimum xtick value from all the samples to plot as the
    # max xtick value on the actual graph (i.e. the graph is truncated to the
    # smallest maximum xtick value from the pool of samples)
    fp <- file(files[[1]], "r")
    xticks <- strtoi(unlist(lapply(strsplit(readLines(fp, n = 1), ",")[[1]],
                                   trimws)))
    close(fp)

    # if there are more
    if (nsamples > 1) {
        for (i in 2:nsamples) {
            fp <- file(files[[i]], "r")
            candidate <- strtoi(
                unlist(
                    lapply(
                        strsplit(
                            readLines(fp, n = 1), ",")[[1]], trimws)))
            close(fp)
            if (tail(candidate, n = 1) < tail(xticks, n = 1)) {
                xticks <- candidate
            }
        }
    }

    # read files
    dataframes <- lapply(files, read.csv, skip = 1)

    # pre-processing & cleaning
    for (i in 1:nsamples) {
        df <- dataframes[[i]]
        df$sample <- rep(sampleNames[[i]], nrow(df))
        df <- df[df$region %in% regions, ]
        dataframes[[i]] <- .summarySE(df, measurevar = 'y',
                                      groupvars = c('x', 'region', 'sample'))
    }

    # merge
    df.union <- do.call("rbind", dataframes)

    # make compound column of region . sample for geom_ribbon
    df.union$compound <- paste(df.union$sample, df.union$region)

    g <- ggplot(df.union, aes(x = x, y = y))

    if (nsamples == 1) {
        g <- g + geom_line(aes(linetype = region, color = sample),
                           color = BLUEHEX, size = 0.65) +
            guides(color = FALSE) +
            #scale_size_manual(values = head(seq(1, length(regions), by = .5),
            #                                n = length(regions))) +
            scale_linetype_manual(values = .getLineTypes(regions))
    } else {
        g <- g +
            geom_line(aes(linetype = region, color = sample), size = 0.65) +
            #scale_size_manual(values = head(seq(1, length(regions), by = .5),
            #                                n = length(regions))) +
            scale_linetype_manual(values = .getLineTypes(regions))
    }

    g <- g + scale_x_continuous(breaks = xticks,
                                limits = c(head(xticks, n = 1),
                                           tail(xticks, n = 1))) +
        geom_ribbon(aes(ymin = y - ci, ymax = y + ci, fill = compound),
                    alpha = 0.1, show.legend = FALSE) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = paste("Rarefaction of", paste(regions, collapse = ", "),
                           "in", paste(sampleNames, collapse = ", ")),
             subtitle = "Mean number of deduplicated sequences with 95% confidence interval",
             x = 'Sample size', y = "Number of deduplicated sequences")
    return(g)
}





#' Title Plots recapture
#'
#' @import ggplot2
#' @include util.R
#'
#' @param files list type. List of _cdr_v_recapture.csv.gz files.
#' @param sampleNames vector type. A vector of strings each
#' representing the name of samples in files.
#' @param regions vector type. A vector of strings,
#' regions to be included in the plot. defaults to c("CDR3", "V")
#'
#' @return ggplot2 object
.plotRecapture <- function(files, sampleNames, regions = c("CDR3", "V")) {
    theme_set(theme_bw())
    nsamples <- length(files)
    # sanity check
    if (nsamples != length(sampleNames)) {
        stop(paste("Expected equal number of sample names and dataframes, got",
                   length(sampleNames), "samples and", nsamples, "dataframes."))
    }

    message(paste("Creating recapture plot for samples",
            paste(sampleNames, collapse = ", ")))

    # find the minimum xtick value from all the samples to plot as the
    # max xtick value on the actual graph (i.e. the graph is truncated to the
    # smallest maximum xtick value from the pool of samples)
    fp <- file(files[[1]], "r")
    xticks <- strtoi(unlist(lapply(strsplit(readLines(fp, n = 1), ",")[[1]],
                                   trimws)))
    close(fp)

    # if there are more
    if (nsamples > 1) {
        for (i in 2:nsamples) {
            fp <- file(files[[i]], "r")
            candidate <- strtoi(
                unlist(
                    lapply(
                        strsplit(
                            readLines(fp, n = 1), ",")[[1]], trimws)))
            close(fp)
            if (tail(candidate, n = 1) < tail(xticks, n = 1)) {
                xticks <- candidate
            }
        }
    }

    # read dataframes
    dataframes <- lapply(files, read.csv, skip = 1)

    # cleanup & pre-processing
    for (i in 1:nsamples) {
        df <- dataframes[[i]]
        # append sample name to a new column named sample
        df$sample <- rep(sampleNames[[i]], nrow(df))
        # only want selected regions - ignore others
        df <- df[df$region %in% regions,]
        # get mean, sd, se, and ci
        dataframes[[i]] <- .summarySE(df, measurevar = 'y',
                                      groupvars = c("x", "region", "sample"))
    }

    df.union <- do.call("rbind", dataframes)

    # make compound column for geom_ribbon (region . sample)
    df.union$compound <- paste(df.union$region, df.union$sample)

    # plot!
    p <- ggplot(df.union, aes(x = x, y = y))

    if (nsamples == 1) {
        p <- p + geom_line(aes(linetype = region, color = sample),
                           color = BLUEHEX, size = 0.65) +
            guides(color = FALSE) +
            #scale_size_manual(values = head(seq(1, length(regions), by = .5),
            #                                n = length(regions))) +
            scale_linetype_manual(values = .getLineTypes(regions))
    } else {
        p <- p +
            geom_line(aes(linetype = region, color = sample), size = 0.65) +
            #scale_size_manual(values = head(seq(1, length(regions), by = .5),
            #                                n = length(regions))) +
            scale_linetype_manual(values = .getLineTypes(regions))
    }

    p <- p + scale_x_continuous(breaks = xticks) +
        geom_ribbon(aes(ymin = y - ci, ymax = y + ci, fill = compound),
                    alpha = 0.1, show.legend = FALSE) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = paste("Percent recapture of",
                           paste(regions, collapse = ", "), "in",
                           paste(sampleNames, collapse = ", ")),
             subtitle = "Mean number of recaptured sequences with 95% confidence interval",
             x = "Sample size", y = "Percent Recapture")
    return(p)
}




#' Title Shows varying regions for a given clonotype defined by its CDR3
#'
#' @import reshape2
#' @import ggplot2
#'
#' @param path string type. Path to diversity folder
#' where <sampleName>_clonotype_diversity_region_analysis.csv.gz is located
#' @param sampleName string type
#' @param top int type. Top N number of clones to analyze
#'
#' @return ggplot2 object
.regionAnalysis <- function(path, sampleName, top = 15) {

    message(paste("Starting clonotype region analysis for", sampleName))

    df <- read.csv(paste0(path, sampleName,
                          "_clonotype_diversity_region_analysis.csv.gz"),
                   stringsAsFactors = FALSE)

    headers <- c("fr1", "cdr1", "fr2", "cdr2", "fr3", "fr4")
    # sort the df with decreasing counts of CDR3 occurance
    df <- df[with(df, order(-count)), ]

    # add new column to sum the "unique" regions
    df$sumcounts = rowSums(df[, headers])

    # grab only those whos V-domain will differ.
    #  This means sumcounts != 6 (> 6) where 6 = length(headers)
    df <- df[df$sumcounts > 6, ]

    # grab top N
    df <- head(df, top)

    # add new row as reference of "unique regions"
    de <- data.frame("REFERENCE", Inf, 1, 1, 1, 1, 1, 1, 6)
    names(de) <- names(df)
    df <- rbind(df, de)


    # reshape to multi-row
    df.mel <- melt(df, measure.vars = headers)
    df.mel[, "value"] = df.mel[, "value"] / df.mel[, "sumcounts"]


    # plot!
    g <- ggplot(df.mel, aes(x = cdr3, y = value, fill = variable,
                            label = sprintf("%0.2f%%",
                                            round(value * 100, digits = 2)))) +
        geom_bar(stat = 'identity') +
        theme(text = element_text(size = 10),
              axis.text.x = element_text(angle = 65, hjust = 1)) +
        #geom_text(position = position_stack(vjust = 0.5)) +
        #stat_summary(fun.y = sum, aes(label = sumcounts, group=cdr3),
        #             geom='text', vjust=-.2) +
        labs(title = paste(sampleName, "varying levels of FRs and CDRs of top",
                           top, "CDR3 clonotype" ),
             subtitle = "Counts of unique regions for a given CDR3",
             x = "CDR3",
             y = "Proportion") +
        guides(fill = guide_legend(title = "Region")) +
        scale_x_discrete(limits = c("REFERENCE", head(df, -1)$cdr3))
    return(g)
}

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
        df$round <- rep(sampleNames[i], nrow(df))
        # normalize percentage to top N
        df$Count <- df$Count / sum(df$Count)
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
    g <- ggplot(df.union, aes(x = round, y = Count)) +
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

#' Title Diversity analysis
#'
#' @import ggplot2
#' @include util.R
#' @include distributions.R
#'
#' @param diversityDirectories list type. List of directories to diversity dir
#' @param diversityOut string type. Output directory
#' @param sampleNames vector type. 1-1 with diversityDirectories
#' @param mashedNames string type. Prefix for output files using "mashed-up"
#' sample names
#' @param .save logical type. Save ggplot object?
#'
#' @return None
.diversityAnalysis <- function(diversityDirectories, diversityOut,
                              sampleNames, mashedNames, .save = T) {
    message(paste("Starting diversity analysis on samples",
            paste(sampleNames, collapse = ", ")))
    # clonotype plots
    cdr3ClonesFile <-
        .listFilesInOrder(path = diversityDirectories,
                          pattern = ".*_cdr3_clonotypes_.*_over\\.csv(\\.gz)?$")

    # these plots can only be done if there are > 1 samples
    if (length(cdr3ClonesFile) > 0 && length(sampleNames) > 1) {
        # plot scatter plot (CDR3 clonotypes)
        .scatterClones(
            lapply(cdr3ClonesFile, read.csv, stringsAsFactors = FALSE),
            sampleNames, diversityOut, "CDR3")

        # plot venn diagram (clonotypes)
        .vennIntersection(
            lapply(cdr3ClonesFile, read.csv, stringsAsFactors = FALSE),
            sampleNames,
            file.path(diversityOut, paste0(mashedNames,
                                           "_cdr3_clonotypeIntersection.png"))
            )

        # plot top N distribution (clonotypes)
        g <- .topNDist(
            lapply(cdr3ClonesFile, read.csv, stringsAsFactors = FALSE),
            sampleNames
            )
        saveName <- file.path(diversityOut, paste0(mashedNames,
                                              "_top10Clonotypes.png"))
        ggsave(saveName, plot = g, width = V_WIDTH_L, height = V_HEIGHT_L)
        .saveAs(.save, saveName, g)
    }

    # fr/cdr plots
    # plot duplication, rarefaction, recapture
    for (reg in c("cdr", "cdr_v", "fr")) {
        if (reg == "cdr") {
            includedRegions <- c("CDR1", "CDR2", "CDR3")
        } else if (reg == "cdr_v") {
            includedRegions <- c("CDR3", "V")
        } else {
            includedRegions <- c("FR1", "FR2", "FR3")
        }
        searchFiles <-
            .listFilesInOrder(path = diversityDirectories,
                              pattern =  paste0(".*_", reg,
                                                "_duplication\\.csv(\\.gz)?$"))
        if (length(searchFiles) > 0) {
            # duplication
            g <- .plotDuplication(searchFiles,
                                 sampleNames,
                                 includedRegions)
            saveName <- file.path(diversityOut, paste0(mashedNames, "_", reg, "_duplication.png"))
            ggsave(saveName, plot = g, width = V_WIDTH, height = V_HEIGHT)
            .saveAs(.save, saveName, g)
        } else {
            warning(paste("Could not find duplication files in",
                          paste(diversityDirectories, collapse = ", ")))
        }

        searchFiles <-
            .listFilesInOrder(path = diversityDirectories,
                              pattern = paste0(".*_", reg,
                                               "_rarefaction\\.csv(\\.gz)?$"))
        if (length(searchFiles) > 0) {
            # plot rarefaction
            g <- .plotRarefaction(searchFiles,
                                 sampleNames,
                                 includedRegions)
            saveName <- file.path(diversityOut, paste0(mashedNames, "_", reg, "_rarefaction.png"))
            ggsave(saveName, plot = g, width = V_WIDTH, height = V_HEIGHT)
            .saveAs(.save, saveName, g)
        } else {
            warning(paste("Could not find rarefaction files in",
                          paste(diversityDirectories, collapse = ", ")))
        }
        searchFiles <-
            .listFilesInOrder(path = diversityDirectories,
                              pattern = paste0(".*_", reg,
                                               "_recapture\\.csv(\\.gz)?$"))
        if (length(searchFiles) > 0) {
            # plot recapture
            g <- .plotRecapture(searchFiles,
                               sampleNames,
                               includedRegions)
            saveName <- file.path(diversityOut, paste0(mashedNames, "_", reg, "_recapture.png"))
            ggsave(saveName,
                   plot = g,
                   width = V_WIDTH, height = V_HEIGHT)
            .saveAs(.save, saveName, g)
        } else {
            warning(paste("Could not find recapture files in",
                          paste(diversityDirectories, collapse = ", ")))
        }
    }

    # MINISECTION:
    # ## SPECTRATYPES ###
    specOut <- file.path(diversityOut, "spectratypes")
    if (!file.exists(specOut)) {
        dir.create(specOut)
    }
    message(paste("Plotting spectratypes on samples",
            paste(sampleNames, collapse = ", ")))
    # CDR 1 - 3
    for (i in 1:3) {
        specFiles <-
            .listFilesInOrder(path = diversityDirectories,
                              pattern = paste0(".*_cdr", i,
                                               "_spectratype\\.csv(\\.gz)?$"))
        if (length(specFiles) > 0) {
            g <- .plotSpectratype(lapply(specFiles, read.csv,
                                         stringsAsFactors = FALSE),
                                  sampleNames,
                                  paste0("CDR", i))
            saveName <- file.path(specOut, paste0(mashedNames, "_cdr", i, "_spectratype.png"))
            ggsave(saveName, plot = g, width = V_WIDTH, height = V_HEIGHT)
            .saveAs(.save, saveName, g)
        } else {
            warning(paste0("Could not find CDR", i, " spectratype files in ",
                          paste(diversityDirectories, collapse = ", ")))
        }
    }

    # special case, no outliers plot for CDR3 only
    specFiles <-
        .listFilesInOrder(path = diversityDirectories,
                          pattern = ".*_cdr3_spectratype_no_outliers\\.csv(\\.gz)?$")

    if (length(specFiles) > 0) {
        g <- .plotSpectratype(lapply(specFiles, read.csv, stringsAsFactors = FALSE),
                             sampleNames,
                             "CDR3")
        saveName <- file.path(specOut, paste0(mashedNames, "_cdr3_spectratype_no_outliers.png"))
        ggsave(saveName, plot = g, width = V_WIDTH, height = V_HEIGHT)
        .saveAs(.save, saveName, g)
    } else {
            warning(paste("Could not find CDR3 spectratype (no outlier) files in",
                          paste(diversityDirectories, collapse = ", ")))
    }

    # FR 1 - 4
    for (i in 1:4) {
        specFiles <-
            .listFilesInOrder(path = diversityDirectories,
                              pattern = paste0(".*_fr", i,
                                               "_spectratype\\.csv(\\.gz)?$"))
        if (length(specFiles) > 0) {
            g <- .plotSpectratype(lapply(specFiles, read.csv, stringsAsFactors = FALSE),
                                 sampleNames,
                                 paste0("FR", i))
            saveName <- file.path(specOut, paste0(mashedNames, "_fr",
                                             i, "_spectratype.png"))
            ggsave(saveName, plot = g, width = V_WIDTH, height = V_HEIGHT)
            .saveAs(.save, saveName, g)
        } else {
            warning(paste0("Could not find FR", i, " spectratype files in ",
                          paste(diversityDirectories, collapse = ", ")))
        }
    }

    # entire V-domain
    specFiles <- .listFilesInOrder(path = diversityDirectories,
                                   pattern = ".*_v_spectratype\\.csv(\\.gz)?$")
    if (length(specFiles) > 0) {
        g <- .plotSpectratype(lapply(specFiles, read.csv, stringsAsFactors = FALSE),
                             sampleNames,
                             "V domain")
        saveName <- file.path(specOut, paste0(mashedNames, "_v_spectratype.png"))
        ggsave(saveName, plot = g, width = V_WIDTH, height = V_HEIGHT)
        .saveAs(.save, saveName, g)
    }

    # MINISECTION:
    # ## COMPOSITION LOGOS ###
    if (length(sampleNames) == 1) {
        compOut <- file.path(diversityOut, "composition_logos")
        if (!file.exists(compOut)) {
            dir.create(compOut)
        }
        message(paste("Plotting composition logos on samples",
                paste(sampleNames, collapse = ", ")))
        .aminoAcidPlot(compOut, sampleNames[1])
    }

    # we can plot region analysis if there's only one sample
    #if (length(sampleNames) == 1) {
    #  # default = top 15
    #  g <- .regionAnalysis(diversityOut, sampleNames[1])
    #  ggsave(paste0(diversityOut, mashedNames, "_region_analysis.png"),
    #         plot = g, width = V_WIDTH, height = V_HEIGHT)
    #}
}


#' Composition logo plot
#'
#' @param compositionDirectory string type.
#' @param sampleName string type.
#' @param regions logical type. vector of FR/CDR regions to plot
#' @param .save logical type. save ggplot object
#'
#' @return none
.aminoAcidPlot <- function(compositionDirectory, sampleName,
                           regions = c("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"),
                           .save = T) {
    for (region in regions) {
        dirName <- file.path(compositionDirectory, region)
        summaryPlot <- file.path(dirName, paste0(sampleName, "_cumulative_logo.csv"))
        df <- read.csv(summaryPlot)
        g1 <- .aminoAcidBar(df, scale = F, region)
        g2 <- .aminoAcidBar(df, scale = T, region)
        fname1 <- file.path(dirName, paste0(sampleName, "_cumulative_logo.png"))
        fname2 <- file.path(dirName, paste0(sampleName, "_cumulative_logo_scaled.png"))
        ggsave(fname1, plot = g1, width = V_WIDTH, height = V_HEIGHT)
        ggsave(fname2, plot = g2, width = V_WIDTH, height = V_HEIGHT)
        .saveAs(.save, fname1, g1)
        .saveAs(.save, fname2, g2)

        germlineSpecific <- list.files(path = dirName, pattern = paste0(sampleName, "_.+_cumulative_logo\\.csv(\\.gz)?$"), full.names = T)
        lapply(germlineSpecific, function(gLogoFile) {
            germName <-
                sub(paste0(dirName, .Platform$file.sep),
                    "",
                    gsub(
                        paste0(sampleName, "_(.*)_cumulative_logo\\.csv$"),
                        "\\1",
                        gLogoFile
                    ))
            df <- read.csv(gLogoFile)
            g1 <- .aminoAcidBar(df, scale = F, region, germ = germName)
            g2 <- .aminoAcidBar(df, scale = T, region, germ = germName)
            fname1 <- file.path(dirName, paste0(sampleName, "_", germName, "_cumulative_logo.png"))
            fname2 <- file.path(dirName, paste0(sampleName, "_", germName, "_cumulative_logo_scaled.png"))
            ggsave(fname1, plot = g1, width = V_WIDTH, height = V_HEIGHT)
            ggsave(fname2, plot = g2, width = V_WIDTH, height = V_HEIGHT)
            .saveAs(.save, fname1, g1)
            .saveAs(.save, fname2, g2)
        })
    }
}

.aminoAcidBar <- function(df, scale, region, germ = "") {
    group.colors <-
        c(
            # oranges
            G = "#ff9900",
            A = "#ffc266",
            S = "#ffe0b3",
            T = "#fff5e6",
            # greens
            C = "#145214",
            V = "#248f24",
            I = "#2eb82e",
            L = "#47d147",
            P = "#70db70",
            F = "#99e699",
            Y = "#c2f0c2",
            M = "#d6f5d6",
            W = "#ebfaeb",
            # purples
            N = "#cc00cc",
            Q = "#ff99ff",
            H = "#ffe6ff",
            # reds
            D = "#ff1a1a",
            E = "#ffb3b3",
            # blues
            K = "#6699ff",
            R = "#e6eeff"
        )
    df.agg <- aggregate(count ~ position, df, sum)

    # get the max counts for each position - then xlabel will contain
    # the amino acid character for that position - break ties on first occurance
    df.max <- merge(aggregate(count ~ position, df, max), df)
    df.max <- df.max[!duplicated(df.max[c(1,2)]), ]
    xlabels <- df.max[with(df.max, order(position)), ]$aa

    total <- max(df.agg$count)
    if (scale) {
        df$proportion <- df$count / total
        subs <- "Scaled to proportion"
    } else {
        df.tmp <- merge(df, df.agg, by = "position")
        df.tmp <- df.tmp[with(df.tmp, order(position)), ]
        # if not scaled, divide within its own position rather than
        # overall (i.e. the max)
        df$proportion <- df.tmp$count.x / df.tmp$count.y
        subs <- ""
    }

    df$aa <-
        factor(
            df$aa,
            levels = c(
                "G",
                "A",
                "S",
                "T",
                "C",
                "V",
                "I",
                "L",
                "P",
                "F",
                "Y",
                "M",
                "W",
                "N",
                "Q",
                "H",
                "D",
                "E",
                "K",
                "R"
            )
        )
    g <- ggplot(df, aes(x = position, y = proportion)) +
        geom_bar(stat = "identity", aes(fill = aa)) +
        labs(title = paste0(germ, " ", region, " (", total, ")"),
             subtitle = subs, x = "amino acid", y = "proportion") +
        scale_x_continuous(breaks = df.agg$position, labels = xlabels) +
        scale_fill_manual(values = group.colors, drop = F) +
        theme(legend.title = element_blank(),
              legend.text = element_text(size = 5))
    return(g)
}
