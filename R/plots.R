

#' A custom ggplot2 theme builder.
#'
#' Can be added to new `ggplot2` plots using `p <- p + plot_theme()`. Further
#' customization is possible by applying a new [ggplot2::theme] to the returned plot
#' (see `?ggplot2::theme()` for details).
#'
#' @param font_size Font size (default: 10)
#' @param font_family Font family (default: "sans")
#' @param line_width Line width (default: 0.5)
#' @param legend_position Legend position (default: "none")
#' @param top_margin Top margin (default: 1 cm)
#' @param bottom_margin Bottom margin (default: 0.25 cm)
#' @param right_margin Right margin (default: 0.1 cm)
#' @param left_margin Left margin (default: 0.25 cm)
#' @param margin_units Margin units (default: "cm")
#' @param device_font Loads font libraries from the users computer. See the [extrafont::loadfonts] documentation for details
#' @returns A ggplot2 theme object
#' @import ggplot2
plot_theme <- function(font_size=10,
                       font_family="sans",
                       line_width=0.5,
                       legend_position="none",
                       top_margin = 1.0,
                       bottom_margin = 0.25,
                       right_margin = 0.1,
                       left_margin = 0.25,
                       margin_units = "cm",
                       device_font = "win"
){
  # Try loading device fonts
  tryCatch(extrafont::loadfonts(device = device_font), error = function(e) invisible(NULL))

  # hjust [0, 0.5, 1] = [left, center, right]
  th <- theme_bw(base_size = font_size, base_family = font_family) +
    theme(plot.margin = margin(t=top_margin, r=right_margin, b=bottom_margin, l=left_margin, unit=margin_units), # Plot margins
          panel.border = element_rect(color = "black", fill = NA, linewidth = line_width), # Plot border
          strip.background = element_blank(), # Strip above/beside plot where facet plot labels are placed
          axis.ticks.x = element_line(colour = "black", linewidth = 0.5*line_width), # x axis ticks
          axis.ticks.y = element_line(colour = "black", linewidth = 0.5*line_width), # y axis ticks
          axis.line = element_line(colour = "black", linewidth = 0.5*line_width), # Axis lines
          plot.title = element_text(size=font_size, face = "bold", color = "black", hjust = 0.5), # Plot title
          axis.text.x = element_text(size=font_size, face = "plain", color = "black"), # x axis labels
          axis.text.y = element_text(size=font_size, face = "plain", color = "black"), # y axis labels
          axis.title = element_text(size=font_size, face = "plain", color = "black"), # Axis titles
          strip.text.x = element_text(size = font_size, face = "plain", colour = "black", angle = 0, hjust = 0.5), # Facet plot x axis labels
          strip.text.y = element_text(size = font_size, face = "plain", colour = "black", angle = -90, hjust = 0.5), # Facet plot y axis labels
          legend.title = element_text(size = font_size, face = "plain", colour = "black"), # Legend title font
          legend.text = element_text(size = font_size, face = "plain", colour = "black"), # Legend text font
          legend.background = element_rect(fill = NA), # Legend border and fill
          # Legend position. Options are "none", "top", "bottom", "left", "right" or a vector of coordinates
          legend.position = legend_position
    )
  return(th)
}

#' Plot one sample and return. Argument `stats` must be one row.
#' Helper function for [plotLDA()]. No input validation done here.
#'
#' @param df The data.frame from results$summary returned from [getMLE()]
#' @param stats The data.frame from results$statistics returned from [getMLE()]
#' @param dilutions Column name in df and stats
#' @param groupby Column name in df and stats or NULL
#' @param signf_dig Number of significant digits to use for statistics label
#' @returns A ggplot
#' @import ggplot2
plot_one <- function(df,
                     stats,
                     dilutions,
                     groupby = NULL,
                     signf_dig = 3){
  # Label for plot
  plot_label <- paste0("MLE: ", signif(100*stats$BC_MLE, signf_dig), "%\n",
                       100*stats$CONF_LEVEL, "% CI [",
                       signif(100*stats$Exact_CI_min, signf_dig), "-",
                       signif(100*stats$Exact_CI_max, signf_dig), "%]\n",
                       "GoF: p = ", signif(stats$Exact_PGOF, 4)) #, "\n",
                       #"Cutoff: >", cutoff, " cells")
  if (is.null(groupby)){
    plot_title <- "Limiting Dilution"
  }else{
    plot_title <- paste0("Limiting Dilution: ", stats[[groupby]])
  }

  # Plot LDA results
  p <- ggplot(df, aes(x = !!sym(dilutions), y = pct_neg)) +
    geom_point() +
    labs(x="Dilution", y="Negative wells (%)") +
    annotate("text", x=max(df[[dilutions]]), y=100, label=plot_label, hjust=1, vjust=1) +
    ylim(0,100) +
    ggtitle(plot_title)

  return(p)
}

#' Plot results of LDA
#'
#' Plots dilutions along the x-axis and % negative wells along the y-axis.
#' Adds a label including the bias-corrected MLE, exact CI, exact GoF p-value,
#' and cutoff threshold.
#'
#' @param results Output from [getMLE()]
#' @param dilutions Column name of dilutions. Should be as provided to [getMLE()]
#' @param groupby Optional character name of column containing sample or group names. Should be as provided to [getMLE()]
#' @param signf_dig Number of significant digits to include in statistics (default: 2)
#' @param nrow Number of rows to plot. Only used when argument 'groupby' provided.
#' @param ncol Number of columns to plot. Only used when argument 'groupby' provided
#' @param ... Additional arguments passed to [plot_theme()]
#' @returns A single ggplot or grid of ggplots
#' @importFrom gridExtra grid.arrange
#' @export
plotLDA <- function(results,
                    dilutions,
                    groupby = NULL,
                    signf_dig = 2,
                    nrow = NULL,
                    ncol = NULL,
                    ...){
  # Validate input
  if (!inherits(results, "list") | !all(names(results) %in% c("summary","statistics"))){
    stop("Argument 'results' should be the output of 'getMLE()'")
  }
  # Get data to plot
  df <- results$summary
  stats <- results$statistics

  # Check columns exist
  if (!(dilutions %in% colnames(df))){
    stop(sprintf("%s must be the column name of the dilutions", dilutions))
  }else if (!is.null(groupby)){
    if (!(groupby %in% colnames(df))) stop(sprintf("%s must be the column name of the groups", groupby))
  }
  # Single plot
  if (is.null(groupby)){
    p <- plot_one(df,
                  stats,
                  dilutions,
                  groupby = groupby,
                  signf_dig = signf_dig) + plot_theme(...)
  # Multiple plots
  }else{
    p <- lapply(unique(df[[groupby]]), function(group){
      # Select group to plot
      df_sub <- df[df[[groupby]] == group,]
      stats_sub <- stats[stats[[groupby]] == group,]
      plot_one(df_sub,
               stats_sub,
               dilutions,
               groupby = groupby,
               signf_dig = signf_dig) + plot_theme(...)
    })

    # Determine grid dimensions
    n_plots <- length(p)

    if (!is.null(ncol) & !is.null(nrow)) {
      grid_ncol <- ncol
      grid_nrow <- nrow
    } else if (!is.null(ncol)) {
      grid_ncol <- ncol
      grid_nrow <- ceiling(n_plots / ncol)
    } else if (!is.null(nrow)) {
      grid_nrow <- nrow
      grid_ncol <- ceiling(n_plots / nrow)
    } else {
      # Auto-calculate reasonable grid dimensions
      grid_ncol <- ceiling(sqrt(n_plots))
      grid_nrow <- ceiling(n_plots / grid_ncol)
    }
    # Arrange plots in grid
    p <- do.call(grid.arrange, c(p,
                                 ncol = grid_ncol,
                                 nrow = grid_nrow))
  }
  return (p)
}

#' Plot MLE with confidence intervals
#'
#' Creates a plot showing the bias-corrected MLE and the upper and lower
#' exact confidence interval for each sample in the column specified by `groupby`.
#'
#' @param results Output from [getMLE()]
#' @param groupby Character name of column containing sample or group names. Should be as provided to [getMLE()]
#' @param mle_colour Line colour of the MLE value (default: "black")
#' @param ci_colour Line colour of the CI box (default: "gray")
#' @param ... Additional arguments passed to [plot_theme()]
#' @returns A ggplot object
#' @import ggplot2
#' @export
plotMLE <- function(results,
                    groupby,
                    width = 0.5,
                    mle_colour = "black",
                    ci_colour = "gray",
                    ...){
  # Validate input
  if (!inherits(results, "list") | !all(names(results) %in% c("summary","statistics"))){
    stop("Argument 'results' should be the output of 'getMLE()'")
  }
  # Get data to plot
  stats <- results$statistics

  # Make percentages
  stats$BC_MLE <- 100 * stats$BC_MLE
  stats$Exact_CI_min <- 100 * stats$Exact_CI_min
  stats$Exact_CI_max <- 100 * stats$Exact_CI_max

  # Make plot
  p <- ggplot(stats, aes(x = !!sym(groupby), y = BC_MLE)) +
    geom_crossbar(aes(ymin = Exact_CI_min, ymax = Exact_CI_max), width = width, middle.colour = mle_colour, box.colour = ci_colour) +
    labs(x="", y="MLE (%)") +
    plot_theme(...)

  return (p)
}
