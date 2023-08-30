#!/usr/bin/env Rscript

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

suppressMessages(library("optparse"))
suppressMessages(library("alluvial"))
suppressMessages(library("data.table"))
suppressMessages(library("svglite"))
suppressMessages(library("stringr"))

option_list <- list(
  make_option(
    c("--cliques"),
    type = "character",
    default = NULL,
    metavar = "path",
    help = "Path to a bedGraph with clique sizes."
  ),
  make_option(
    c("--highlight_label"),
    type = "character",
    default = NULL,
    metavar = "string",
    help = "Condition to highlight."
  ),
  make_option(
    c("--highlight_color"),
    type = "character",
    default = "orange",
    metavar = "character",
    help = "Color used for highlighting."
  ),
  make_option(
    c("--base_color"),
    type = "character",
    default = "grey",
    metavar = "character",
    help = "Default color."
  ),
  make_option(
    c("-o", "--outprefix"),
    type = "character",
    default = NULL,
    metavar = "path",
    help = "Output prefix."
  )
)

generate_links <-
  function(df) {
    cols <- setdiff(names(df), c("chrom", "start", "end", "size"))
    df <- df[, cols]

    return(df)
  }

generate_frequencies <-
  function(df) {
    return(df$size)
  }

generate_layer_mask <-
  function(df, label) {
    if (is.null(label)) {
      return(rep(TRUE, length(df$size)))
    }

    col <- names(df)[[1]]
    return(df[, col] != label)
  }

generate_colors <-
  function(mask, base_color, highlight_color) {
    return(ifelse(mask, highlight_color, base_color))
  }

plot_alluvial <-
  function(counts, outprefix, highlight_label, highlight_color, base_color, alpha = 0.8, blocks = FALSE) {
    links <- generate_links(counts)
    freq <- generate_frequencies(counts)
    layer_mask <- generate_layer_mask(counts, highlight_label)
    colors <- generate_colors(layer_mask, highlight_color, base_color)

    svglite(
      file = paste(outprefix, "svg", sep = "."),
      width = 10,
      height = 7
    )
    alluvial(
      links,
      freq = freq,
      layer = layer_mask,
      col = colors,
      border = colors,
      alpha = alpha,
      blocks = blocks
    )
    dev.off()

    return(invisible())
  }

# Parse CLI options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

missing_options <- list()

if (is.null(opt$cliques)) {
  missing_options <- append(missing_options, "--cliques")
}
if (is.null(opt$outprefix)) {
  missing_options <- append(missing_options, "--outprefix")
}

if (length(missing_options) != 0) {
  num_missing_options <- length(missing_options)
  missing_options <- paste(missing_options, "\n - ")
  stop(
    stringr::str_interp(
      "The following ${num_missing_options} mandatory option(s) are missing:\n${missing_options}"
    ),
    .call = FALSE
  )
}

counts <- as.data.frame(fread(opt$cliques))
plot_alluvial(counts, opt$outprefix, opt$highlight_label, opt$highlight_color, opt$base_color)
