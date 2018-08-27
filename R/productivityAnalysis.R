#' Summary of productivity
#'
#' @import ggplot2
#'
#' @param dataframes list type. List of sample dataframes
#' @param sampleNames vector type. 1-1 with dataframes
#'
#' @return ggplot2 object
.productivityPlot <- function(dataframes, sampleNames) {
    nsamples <- length(dataframes)
    if (length(sampleNames) != nsamples) {
        stop(paste("Expected equal number of sample names and dataframes, got",
                length(sampleNames), "samples and", nsamples, "dataframes."))
    }

    message(paste("Plotting productivity summary plot for",
                   paste(sampleNames, collapse = ", ")))

    #   ---- clean & pre-processing ----
    unusedCols <- c("X")
    # label each row with sample name and drop unused columns
    for (i in seq_len(nsamples)) {
        df <- dataframes[[i]]
        df$sample <- sampleNames[i]
        df[df == "Both"] <- "Stop & Out"
        dataframes[[i]] <- df[, !(names(df) %in% unusedCols)]
    }
    #   ---- done: clean & pre-processing ----

    # merge all samples
    df.union <- do.call("rbind", dataframes)

    # plot!
    g <- ggplot(df.union, aes(sample, Percentage)) +
        geom_bar(stat = "identity", aes(fill = Reason), width = 0.5) +
        facet_wrap(~Productivity) +
        labs(
            title = "Productivity proportions",
            subtitle = "Percentage of unproductive reads due to stop codons and frameshifts",
            y = "Percentage (%)"
        ) +
        coord_cartesian(ylim = c(0, 100)) +
        geom_text(position = position_stack(vjust = 0.5),
                  aes(label = sprintf("%0.2f%%", Percentage),
                      group = Reason)) +
        theme(axis.text.x = element_text(angle = 55, hjust = 1),
              legend.title = element_blank(), axis.title.x = element_blank())
    return(g)
}

#' Plots a distribution plot for different productivity analysis files
#'
#' @description A wrapper for plotDist
#'
#' @include util.R
#' @include distributions.R
#' @import ggplot2
#'
#' @param productivityDirectories vector type.
#' directories where all productivity csv
#' files lives (usually <samplename>/productivity/)
#' @param sampleNames vector type.
#' @param title string type.
#' @param reg string type. Regular expression to find the right files
#' for this particular distribution plot
#' @param outputFileName string type.
#' Vector of file names to save in the order of \code{regions}
#' @param region string type. Most of the dist plots are
#' regional based. use "" if no regions are involved
#' @param .save logical type. Save Rdata?
#'
#' @return None
.prodDistPlot <-
    function(productivityDirectories,
             sampleNames, title, reg, outputFileName,
             region, .save = TRUE) {
    fs <- .listFilesInOrder(path = productivityDirectories,
                            pattern = paste0(".*", region, reg))
    if (length(fs)) {
        dataframes <- lapply(fs, read.csv,
                             stringsAsFactors = FALSE, skip = 1)
        plotTitle <- paste(title, toupper(region), "in",
                           paste(sampleNames, collapse = ", "))
        subtitle <-
            paste("Total is" , paste(lapply(fs, function(x) {
                as.integer(.getTotal(x))
            }), collapse = ", "))
        g <- .plotDist(dataframes,
                       sampleNames,
                       plotTitle,
                       .checkVert(fs[[1]]),
                       subs = subtitle)
        ggsave(outputFileName, plot = g, width = V_WIDTH, height = V_HEIGHT)
        .saveAs(.save, outputFileName, g)
    }
}


#' Conducts productivty analysis
#'
#' @include util.R
#' @include distributions.R
#' @import ggplot2
#'
#' @param productivityDirectories list type. List of directories
#' @param prodOut string type. Output directory
#' @param sampleNames vector type. 1-1 with productivity directories
#' @param combinedNames string type. Title friendly "combined" sample names
#' @param mashedNames string type. File friendly "mashed-up" sample names
#' @param .save logical type. Save Rdata
#'
#' @return None
.productivityAnalysis <- function(productivityDirectories,
                                  prodOut,
                                  sampleNames,
                                  combinedNames,
                                  mashedNames,
                                  .save = TRUE) {
    message("Commencing productivity analysis")

    # summary productivity file
    prodFiles <- .listFilesInOrder(path = productivityDirectories,
                                   pattern = ".*_productivity\\.csv(\\.gz)?$")


    if (length(prodFiles) > 0) {
        g <- .productivityPlot(lapply(prodFiles, read.csv,
                                      stringsAsFactors = FALSE),
                               sampleNames)
        saveName <- file.path(prodOut, paste0(mashedNames, "_productivity.png"))
        ggsave(saveName, plot = g, width = V_WIDTH, height = V_HEIGHT)
        .saveAs(.save, saveName, g)
    }


    # sub-productivity files
    regions <- c("cdr1",
                 "cdr2",
                 "cdr3",
                 "fr1",
                 "fr2",
                 "fr3",
                 "igv",
                 "igd",
                 "igj")

    # gaps_dist plots only ---------------------------------------------------
    message(paste("Plotting indel distributions for", combinedNames))
    lapply(regions, function(region) {
        .prodDistPlot(
            productivityDirectories,
            sampleNames,
            "Gaps in",
            "_gaps_dist\\.csv(\\.gz)?$",
            file.path(prodOut, paste0(mashedNames, "_",
                                      region,
                                      "_gaps_dist.png")),
            region)
    })

    # gaps_out_of_frame plots only (no igv, igd, ihj plots for this) ---------
    subregions <- head(regions, n = 6)
    message(paste(
        "Plotting indel out of frame distributions for",
        combinedNames
    ))

    lapply(subregions, function(region) {
        .prodDistPlot(
            productivityDirectories,
            sampleNames,
            "Gaps in",
            "_gaps_dist_out_of_frame\\.csv(\\.gz)?$",
            file.path(prodOut, paste0(mashedNames, "_", region,
                                      "_gaps_dist_out_of_frame.png")),
            region)
    })

    # mismatch dist only -----------------------------------------------------
    message(paste("Plotting mismatches distributions for", combinedNames))
    lapply(regions, function(region) {
        .prodDistPlot(
            productivityDirectories,
            sampleNames,
            "Mismatches in",
            "_mismatches_dist\\.csv(\\.gz)?$",
            file.path(prodOut, paste0(mashedNames,
                                      "_", region, "_mismatches_dist.png")),
            region)
    })

    # stop codon dist plot ---------------------------------------------------
    message(paste(
        "Plotting stop codon in-frame distributions for",
        combinedNames
    ))
    .prodDistPlot(
        productivityDirectories,
        sampleNames,
        "Stop codon in In-frame Clones",
        "_stopcodon_dist_in_frame\\.csv(\\.gz)?$",
        file.path(prodOut, paste0(mashedNames, "_stopcodon_dist_in_frame.png")),
        ""
    ) # no regions

    # vjframe plot -----------------------------------------------------------
    message(paste("Plotting vjframe distributions for", combinedNames))
    .prodDistPlot(
        productivityDirectories,
        sampleNames,
        "V-D-J Rearrangement",
        "_vjframe_dist\\.csv(\\.gz)?$",
        file.path(prodOut, paste0(mashedNames, "_vjframe_dist.png")),
        ""
    ) # no regions

    # 3 special cases for IGV region -----------------------------------------
    # igv - inframe_unproductive, out of frame, productive

    # in frame unprod ----------
    message(paste(
        "Plotting inframe unproductive distributions for",
        combinedNames
    ))

    .prodDistPlot(
        productivityDirectories,
        sampleNames,
        "Abundance of In-frame Unproductive Clones in",
        "_dist_inframe_unproductive\\.csv(\\.gz)?$",
        file.path(prodOut, paste0(mashedNames,
                                  "_igv_dist_inframe_unproductive.png")),
        "igv"
    )

    # out of frame unprod --------
    message(paste(
        "Plotting out of frame clones distributions for",
        combinedNames
    ))
    .prodDistPlot(
        productivityDirectories,
        sampleNames,
        "Abundance of Out-Of-Frame Clones in",
        "_dist_out_of_frame\\.csv(\\.gz)?$",
        file.path(prodOut, paste0(mashedNames, "_igv_dist_out_of_frame.png")),
        "igv"
    )

    # productive ----------------
    message(paste(
        "Plotting out IGV productive distributions for",
        combinedNames
    ))
    .prodDistPlot(
        productivityDirectories,
        sampleNames,
        "Abundance of Productive Clones in",
        "_dist_productive\\.csv(\\.gz)?$",
        file.path(prodOut, paste0(mashedNames, "_igv_dist_productive.png")),
        "igv"
    )

    # stop codon in FR/CDR region proportion plot ----------------------------
    message(paste(
        "Plotting stop codon region distributions for",
        combinedNames
    ))
    lapply(c("inframe", "outframe"), function(framestat) {
        stopCodonRegionFiles <-
            .listFilesInOrder(
                path = productivityDirectories,
                pattern = paste0(".*_stopcodon_region_",
                                 framestat, "\\.csv(\\.gz)?$")
            )
        if (length(stopCodonRegionFiles) > 0) {
            vert = .checkVert(stopCodonRegionFiles[[1]])
            if (framestat == 'inframe') {
                titlestatus <- "In-frame"
            } else {
                titlestatus <- "Out-of-frame"
            }
            subtitle <- paste("Total is",
                              paste(
                                  lapply(stopCodonRegionFiles, function(x) {
                                      as.integer(.getTotal(x))
                                  }),
                                  collapse = ", "
                              ))

            stopcodonRegion <- .plotDist(
                lapply(stopCodonRegionFiles, read.csv, skip = 1),
                sampleNames,
                paste(
                    "Stop codon in FRs and CDRs of",
                    titlestatus,
                    "sequences in",
                    combinedNames
                ),
                vert,
                subs = subtitle,
                sorted = FALSE
            )
            saveName <- file.path(prodOut,
                                  paste0(
                                      mashedNames,
                                      "_stopcodon_region_",
                                      framestat,
                                      ".png"
                                  ))
            if (vert) {
                ggsave(saveName, plot = stopcodonRegion,
                       width = V_WIDTH, height = V_HEIGHT)
            } else {
                ggsave(saveName, plot = stopcodonRegion,
                       width = H_WIDTH, height = H_HEIGHT)
            }
            .saveAs(.save, saveName, stopcodonRegion)
        }
    })
}
