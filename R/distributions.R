.plotDist <- function(dataframes, sampleNames, plotTitle, vert = T,
                      xlabel = "", ylabel = "", perc = T, subtitle = "",
                      sorted = T, cutoff = 15) {
    require(ggplot2)
    frames <- length(dataframes)
    # sanity check
    if (length(sampleNames) != frames) {
        stop(paste("Expected equal number of sample names and dataframes, got",
                   length(sampleNames), "samples and", frames, "dataframes."))
    }

    # If there was a cutoff, caps will display the right message,
    #  by default it's empty
    caps <- ""

    # if it was percentage, put a % sign in labels, otherwise just normal values
    if (perc) {
        placeholder <- "%0.2f%%"
    } else {
        placeholder <- "%0.2f"
    }

    # ------------  cleaning & transforming data   ----------- #
    for (i in seq_len(frames)) {

        df <- dataframes[[i]]

        # for each dataframe, append a new column 'sample'
        # to indicate where this df came from
        df$sample <- rep(sampleNames[i], nrow(df))

        # sort the values by descending order (strictly speaking,
        # not necessary because it's already sorted
        # - but better be safe than sorry)
        if (sorted) {
            if (vert) {
                dataframes[[i]] <- df[with(df, order(-y)), ]
            } else {
                dataframes[[i]] <- df[with(df, order(-x)), ]
            }
        } else {
            dataframes[[i]] <- df
        }

        # canonicalize NaN, NA, and N/As to NA.
        dataframes[[i]][is.na(dataframes[[i]])] <- "NA"
        dataframes[[i]][dataframes[[i]] == "N/A"] <- "NA"
    }

    # keep a copy of the original list of dataframes without cutoff applied
    originals <- dataframes

    # make sure only the top cutoff is considered
    for (i in seq_len(frames)) {
        df <- dataframes[[i]]
        if (nrow(df) > cutoff) {
            caps <- "Cutoff at top 15"
            dataframes[[i]] <- head(df, cutoff)
        }
    }

    # --------- complete: cleaning & transforming data -------- #


    # merge dataframes into one
    df.union <- do.call("rbind", dataframes)


    # need to compensate for topN in each dataframe
    # 1. gather unique x/y
    # 2. if x/y is in original table but not in top N table,
    #    append row to df.union
    # 3. otherwise, do nothing

    if (frames > 1) {
        if (vert) {
            topNx <- unique(df.union$x)
            for (x_ in topNx) {
                for (i in seq_len(frames)) {
                    if (x_ %in% originals[[i]]$x
                        && !(x_ %in% dataframes[[i]]$x)) {
                        newRow <- originals[[i]][originals[[i]]$x == x_, ]
                        df.union <- rbind(df.union, newRow)
                    }
                }
            }
        } else {
            topNy <- unique(df.union$y)
            for (y_ in topNy) {
                for (i in seq_len(frames)) {
                    if (y_ %in% originals[[i]]$y
                        && !(y_ %in% dataframes[[i]]$y)) {
                        newRow <- originals[[i]][originals[[i]]$y == y_, ]
                        df.union <- rbind(df.union, newRow)
                    }
                }
            }
        }
    }

    if (missing(ylabel)) {
        ylabel <- "Proportion"
        if (perc) {
            ylabel <- paste(ylabel, "(%)")
        }
    }
    if (!vert) {
        if (sorted) {
            g <- ggplot(df.union,
                        aes(x = reorder(y, x), y = x,
                            label = sprintf(placeholder, x)))
            + coord_flip()
        } else {
            df.union$y <- factor(df.union$y, levels = unique(df.union$y))
            g <- ggplot(df.union,
                        aes(x = y, y = x,
                            label = sprintf(placeholder, x)))
            + coord_flip()
        }

        if (frames == 1) {
            # single sample -> blue colour plot
            g <- g + geom_text(hjust = 0.50,
                               vjust = -0.5, size = 3, angle = -90) +
                geom_bar(stat = "identity", aes(fill = sample),
                         position = "dodge",
                         fill = BLUEHEX, show.legend = FALSE)
        } else {
            # multiple samples -> multi-coloured plot
            g <- g + geom_bar(stat = "identity",
                              aes(fill = sample), position = "dodge")
        }

    } else {
        if (sorted) {
            g <- ggplot(df.union,
                        aes(x = reorder(x, -y), y = y,
                            label = sprintf(placeholder, y)))
        } else {
            df.union$x <- factor(df.union$x, levels = unique(df.union$x))
            g <- ggplot(df.union, aes(x = x, y = y,
                                      label = sprintf(placeholder, y)))
        }

        if (frames == 1) {
            # single sample -> blue colour plot
            g <- g + geom_text(vjust = -0.5, size = 3) +
                geom_bar(stat = 'identity', aes(fill = sample),
                         position = 'dodge',
                         fill = BLUEHEX, show.legend = FALSE)
        } else {
            # multiple samples -> multi-coloured plot
            g <- g + geom_bar(stat = 'identity',
                              aes(fill = sample), position = 'dodge')
        }
    }
    g <- g + theme(text = element_text(size = 10)) +
        labs(title = plotTitle,
             subtitle = subtitle,
             x = xlabel,
             y = ylabel,
             caption = caps)
    return(g)
}
