#!/usr/bin/env Rscript

suppressMessages(library("optparse"))
suppressMessages(library("alluvial"))
suppressMessages(library("data.table"))
suppressMessages(library("svglite"))
suppressMessages(library("stringr"))

option_list <- list(
  make_option(
    "--compartments",
    type = "character",
    default = NULL,
    metavar = "path",
    help = "Path to a subcompartment file bedGraph produced by dcHiC."
  ),
  make_option(
    "--highlight_label",
    type = "character",
    default = NULL,
    metavar = "string",
    help = "Condition to highlight."
  ),
  make_option(
    "--only_show_highlighted",
    type = "logical",
    action = "store_true",
    default = FALSE,
    help = "Hide non-higlighted ribbons."
  ),
  make_option(
    "--highlight_color",
    type = "character",
    default = "orange",
    metavar = "character",
    help = "Color used for highlighting."
  ),
  make_option(
    "--base_color",
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
  ),
  make_option(
    "--width",
    type = "double",
    default = 10,
    metavar = "double",
    help = "Plot width."
  ),
  make_option(
    "--height",
    type = "double",
    default = 7,
    metavar = "double",
    help = "Plot height."
  )
)

generate_links <-
  function(df) {
    cols <- grep(".*.state$", names(df), value = TRUE)
    df <- df[, cols]
    names(df) <- lapply(names(df), sub, pattern = ".state$", replacement = "")

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
  function(counts, outprefix, highlight_label, highlight_color, base_color, only_show_higlighted, width, height, blocks = TRUE, alpha = 0.8) {
    if (base_color == "white") {
      only_show_higlighted <- TRUE
      alpha <- 1.0
    }

    links <- generate_links(counts)
    freq <- generate_frequencies(counts)
    layer_mask <- generate_layer_mask(counts, highlight_label)
    colors <- generate_colors(layer_mask, highlight_color, base_color)

    if (only_show_higlighted) {
      mask <- layer_mask
    } else {
      mask <- replicate(length(layer_mask), FALSE)
    }

    svglite(
      file = paste(outprefix, "svg", sep = "."),
      width = width,
      height = height
    )
    alluvial(
      links,
      freq = freq,
      layer = layer_mask,
      col = colors,
      hide = mask,
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

if (is.null(opt$compartments)) {
  missing_options <- append(missing_options, "--compartments")
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

counts <- as.data.frame(fread(opt$compartments))

plot_alluvial(counts,
              opt$outprefix,
              opt$highlight_label,
              opt$highlight_color,
              opt$base_color,
              opt$only_show_highlighted,
              opt$width,
              opt$height)
