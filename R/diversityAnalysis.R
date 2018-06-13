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

#' Composition logo plot
#'
#' @param compositionDirectory string type.
#' @param outdir string type.
#' @param sampleName string type.
#' @param regions logical type. vector of FR/CDR regions to plot
#' @param .save logical type. save ggplot object
#'
#' @import stringr
#'
#' @return none
.aminoAcidPlot <- function(compositionDirectory, outdir, sampleName,
                           regions = c("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"),
                           .save = T) {
    for (region in regions) {
        # todo: change dirname such that compositionDirectory and outdir
        # is differentiated!
        dirName <- file.path(compositionDirectory, region)
        outputPath <- file.path(outdir, region)
        if (!dir.exists(outputPath)) {
            dir.create(outputPath)
        }

        summaryPlot <- file.path(dirName, paste0(sampleName, "_cumulative_logo.csv"))
        df <- read.csv(summaryPlot)
        g1 <- .aminoAcidBar(df, scale = F, region)
        g2 <- .aminoAcidBar(df, scale = T, region)
        fname1 <- file.path(outputPath, paste0(sampleName, "_cumulative_logo.png"))
        fname2 <- file.path(outputPath, paste0(sampleName, "_cumulative_logo_scaled.png"))
        ggsave(fname1, plot = g1, width = V_WIDTH, height = V_HEIGHT)
        ggsave(fname2, plot = g2, width = V_WIDTH, height = V_HEIGHT)
        .saveAs(.save, fname1, g1)
        .saveAs(.save, fname2, g2)

        germlineSpecific <-
            list.files(path = dirName,
                       pattern = paste0(sampleName,
                                        "_.+_cumulative_logo\\.csv(\\.gz)?$"),
                       full.names = T)

        lapply(germlineSpecific, function(gLogoFile) {
            germName <- sub("_cumulative_logo\\.csv(\\.gz)?$", "",
                            stringr::str_extract(gLogoFile, "IG[HKL][VDJ].*"))
            df <- read.csv(gLogoFile)
            g1 <- .aminoAcidBar(df, scale = F, region, germ = germName)
            g2 <- .aminoAcidBar(df, scale = T, region, germ = germName)
            fname1 <- file.path(outputPath, paste0(sampleName, "_", germName, "_cumulative_logo.png"))
            fname2 <- file.path(outputPath, paste0(sampleName, "_", germName, "_cumulative_logo_scaled.png"))
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
            G = "#e65c00",
            A = "#ff751a",
            S = "#ff8533",
            T = "#ff944d",
            # greens
            C = "#003300",
            V = "#145214",
            I = "#006622",
            L = "#1f7a1f",
            P = "#009933",
            F = "#29a329",
            Y = "#00b33c",
            M = "#2eb82e",
            W = "#33cc33",
            # purples
            N = "#330066",
            Q = "#4d0099",
            H = "#6600cc",
            # reds
            D = "#990000",
            E = "#b30000",
            # blues
            K = "#000099",
            R = "#0000cc"
        )
    df.agg <- aggregate(count ~ position, df, sum)

    # get the max counts for each position - then xlabel will contain
    # the amino acid character for that position - break ties on first occurance
    df.max <- merge(aggregate(count ~ position, df, max), df)
    df.max <- df.max[!duplicated(df.max[c(1,2)]), ]
    xlabels <- lapply(df.max[with(df.max, order(position)), ]$aa, as.character)

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

##  dataframe is something like:
##  +--------------------------+
##  | Clonotype | Count | prop |
##  +--------------------------+
##  |           |       |      |
##  |           |       |      |
##  |           |       |      |
##  |  .......  |  ...  | .... |
##  +--------------------------+
#' Reports abundance-based (Lower bound) diversity estimates using the Vegan package
#'
#' @import vegan
#'
#' @param df clonotype dataframe. Vegan format:
#' +---------------------------+
#' | S.1| S.2| S.3 | S.4 | ... |   (each species should have its own column)
#' +---------------------------+
#' | v1 |v2  | v3  | ....      |   (each species' count values are placed in the corresponding column)
#' +---------------------------+
#'
#' @return dataframe with the format:
#' +----------------------------------------------------------------+
#' | S.obs | S.chao1 | se.chao1 | S.ACE | se.ACE | s.jack1 | s.jack2|
#' +----------------------------------------------------------------+
#' | v1    |  v2 ....                                               |
#' +----------------------------------------------------------------+
.reportLBE <- function(df) {
    f1.f2 <- unlist(lapply(1:2, function(i) {
        sum(df[1, ] == i)
    }))

    lbe <- estimateR(df)
    s.obs <- lbe[1]
    s.jack1 <- s.obs + f1.f2[1]
    s.jack2 <- s.obs + 2 * f1.f2[1] - f1.f2[2]

    lbe <- rbind(as.data.frame(lbe), S.jack1 = s.jack1, S.jack2 = s.jack2)

    df.out <- as.data.frame(t(lbe[, 1]))
    names(df.out) <- rownames(lbe)
    return(df.out)
}

#' Calculates the "standard" diversity indices
#'
#' @import vegan
#'
#' @param df clonotype dataframe. Vegan format:
#' +---------------------------+
#' | S.1| S.2| S.3 | S.4 | ... |   (each species should have its own column)
#' +---------------------------+
#' | v1 |v2  | v3  | ....      |   (each species' count values are placed in the corresponding column)
#' +---------------------------+
#'
#' @return dataframe with the column headers:
#' shannon , simpson.con , simpson.inv , simpson.gini , renyi.0 ,
#' renyi.1 , renyi.2 , renyi.Inf , hill.0 , hill.1 , hill.2 , hill.Inf
#'
#' renyi.0 => species richness
#' renyi.1 => shannon entropy
#' renyi.2 => inv.gini
#' renyi.Inf => min.entropy
#'
#' finally:
#'    hill_a = exp(renyi_a)
#'
.calculateDInd <- function(df) {
    renyi.scales = c(0, 1, 2, Inf)
    renyi <- renyi(df, scales = renyi.scales)
    hill <- exp(renyi)
    shannon <- vegan::diversity(df, index = "shannon")
    simpson <- vegan::diversity(df, index = "simpson")
    n.species <- ncol(df)

    as.data.frame(cbind(
          shannon = shannon,
          shannon.norm = shannon / log(n.species),
          simpson.gini = simpson,
          simpson.inv = vegan::diversity(df, index = "invsimpson"),
          simpson.con = 1 - simpson,
          renyi.0 = renyi['0'],
          renyi.1 = renyi['1'],
          renyi.2 = renyi['2'],
          renyi.Inf = renyi['Inf'],
          hill.0 = hill['0'],
          hill.1 = hill['1'],
          hill.2 = hill['2'],
          hill.Inf = hill['Inf']
          ))
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
        compDir <- file.path(diversityDirectories[[1]], "composition_logos")
        compOut <- file.path(diversityOut, "composition_logos")
        if (!file.exists(compOut)) {
            dir.create(compOut)
        }
        message(paste("Plotting composition logos on samples",
                paste(sampleNames, collapse = ", ")))
        .aminoAcidPlot(compDir, compOut, sampleNames[1])
    }

    # we can plot region analysis if there's only one sample
    #if (length(sampleNames) == 1) {
    #  # default = top 15
    #  g <- .regionAnalysis(diversityOut, sampleNames[1])
    #  ggsave(paste0(diversityOut, mashedNames, "_region_analysis.png"),
    #         plot = g, width = V_WIDTH, height = V_HEIGHT)
    #}

    ####################################################
    #                  UNSEEN SPECIES                  #
    ####################################################
    cdr3ClonesFile <- .listFilesInOrder(path = diversityDirectories,
                                        pattern = ".*_cdr3_clonotypes_.*_over\\.csv(\\.gz)?$")
    lb.fname <- "lower_bound_estimate.tsv"
    ind.fname <- "diversity_indices.tsv"

    if (length(cdr3ClonesFile)) {
        # dataframes is in vegan input format, the clonotypes are now column headers
        # instead of column values
        dataframes <- lapply(cdr3ClonesFile, function(fname) {
            df <- read.csv(fname, stringsAsFactors = F)
            d.trans <- as.data.frame(t(df[, "Count"]))
            names(d.trans) <- df$Clonotype
            return(d.trans)
        })
        lbeOut <- file.path(diversityOut, lb.fname)
        indOut <- file.path(diversityOut, ind.fname)

        if (length(dataframes) == 1) {
            ##### Lower bound estimates of unseen species
            message(paste("Calculating unseen species lower bound estimates for",
                          sampleNames[1]))
            df.lbe <- .reportLBE(dataframes[[1]])
            write.table(df.lbe, file = lbeOut, sep = "\t", quote = F, row.names = F)

            ##### Diversity indices
            message(paste("Calculating standard diversity indices for",
                          sampleNames[1]))
            df.ind <- .calculateDInd(dataframes[[1]])
            write.table(df.ind, file = indOut, sep = "\t", quote = F, row.names = F)

        } else {
            ##### Lower bound estimates of unseen species
            # don't bother recalculating .reportLBE unless the file can't be found
            # .reportLBE might not be calculated if the user manually combined repertoires
            # before actually calling abseqR::report() on the individual repertoires

            # try to find all lower bound estimate files in all the sample directories, if
            # even one of them does not exist, we'll have to generate the tsv file
            lbesFiles <- .listFilesInOrder(path = diversityDirectories,
                                               pattern = paste0(lb.fname, "\\.tsv(\\.gz)?$"))
            if (length(lbesFiles) != length(sampleNames)) {
                message(paste("Calculating unseen species lower bound estimates for",
                              paste(sampleNames, collapse = ", ")))
                df.lbes <- lapply(dataframes, .reportLBE)
                stopifnot(length(df.lbes) == length(sampleNames) &&
                              length(sampleNames) == length(diversityDirectories))

                # write out to individual sample directories
                lapply(seq_along(df.lbes), function(i) {
                    write.table(df.lbes[[i]],
                                file = file.path(diversityDirectories[[i]], lb.fname),
                                sep = "\t",
                                quote = F,
                                row.names = F)
                })
            } else {
                # expect 1-1 relationship
                message(paste("Loading precomputed lower bound estimates from",
                              paste(sampleNames, collapse = ", ")))
                stopifnot(length(diversityDirectories) == length(sampleNames))
                lbesFiles <- .listFilesInOrder(path = diversityDirectories,
                                               pattern = paste0(lb.fname, "\\.tsv(\\.gz)?$"))
                stopifnot(length(lbesFiles) == length(sampleNames))
                df.lbes <- lapply(lbesFiles, read.table, header = T)
            }
            dfs <- do.call("rbind", lapply(seq_along(df.lbes), function(i) {
                df <- df.lbes[[i]]
                df$sample <- sampleNames[i]
                return(df)
            }))
            write.table(dfs, file = lbeOut, sep = "\t", quote = F, row.names = F)

            ##### Diversity indices
            # same logic as above
            indFiles <- .listFilesInOrder(path = diversityDirectories,
                                          pattern = paste0(ind.fname, "\\.tsv(\\.gz)?$"))
            # if we can't recover all the required indices tsv file, re-generate
            # all of them
            if (length(indFiles) != length(sampleNames)) {
                message(paste("Calculating standard diversity indices for",
                              paste(sampleNames, collapse = ", ")))

                df.inds <- lapply(dataframes, .calculateDInd)
                stopifnot(length(df.inds) == length(sampleNames) &&
                              length(sampleNames) == length(diversityDirectories))

                # write out to individual sample directories
                lapply(seq_along(df.inds), function(i) {
                    write.table(df.inds[[i]],
                                file = file.path(diversityDirectories[[i]],
                                                 ind.fname),
                                sep = "\t",
                                quote = F,
                                row.names = F)
                })
            } else {
                message(paste("Loading precomputed diversity indices from",
                              paste(sampleNames, collapse = ", ")))
                stopifnot(length(diversityDirectories) == length(sampleNames))
                indFiles <- .listFilesInOrder(path = diversityDirectories,
                                              pattern = paste0(ind.fname, "\\.tsv(\\.gz)?$"))
                stopifnot(length(indFiles) == length(sampleNames))
                df.inds <- lapply(indFiles, read.table, header = T)
            }
            dfs <- do.call("rbind", lapply(seq_along(df.inds), function(i) {
                df <- df.inds[[i]]
                df$sample <- sampleNames[i]
                return(df)
            }))
            write.table(dfs, file = indOut, sep = "\t", quote = F, row.names = F)
        }

    } else {
        warning(paste(paste(sampleNames, collapse = ", "),
                      "is missing CDR3 clonotype counts file,",
                      "skipping LBE and IND analysis"))
    }
}
