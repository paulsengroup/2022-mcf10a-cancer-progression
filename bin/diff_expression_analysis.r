#!/usr/bin/env Rscript

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

suppressMessages(library("optparse"))
suppressMessages(library("DESeq2"))
suppressMessages(library("stringr"))
suppressMessages(library("fs"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("pheatmap"))
suppressMessages(library("parallel"))
suppressMessages(library("BiocParallel"))

option_list <- list(
  make_option(
    c("-i", "--count_matrix"),
    type = "character",
    default = NULL,
    metavar = "path",
    help = "Path to one of the count matrices produce dy nf-core/rnaseq"
  ),
  make_option(
    c("-c", "--contrast"),
    type = "character",
    default = NULL,
    metavar = "string",
    help = "Sample name to use as contrast."
  ),
  make_option(
    c("--min_counts"),
    type = "numeric",
    default = 10,
    metavar = "numeric",
    help = "Cutoff used to remove low count genes."
  ),
  make_option(
    c("-o", "--outprefix"),
    type = "character",
    default = NULL,
    metavar = "path",
    help = "Output prefix."
  ),
  make_option(
    c("-r", "--sample_suffix"),
    type = "character",
    default = "_REP\\d$",
    metavar = "string",
    help = "Regex used to strip sample suffixes (e.g. repl number)"
  ),
  make_option(
    c("-p", "--cpu_cores"),
    type = "integer",
    default = 1,
    metavar = "integer",
    help = "Number of cores."
  )
)

construct_dds_from_count_matrix <-
  function(path_to_count_matrix,
           contrast,
           sample_suffix,
           row_names = "id",
           round = TRUE,
           seq_types = "paired-end",
           min_counts = 10) {
    num_conditions <-
      length(read.table(path_to_count_matrix,
                        nrows = 1)) - 1

    colClasses <-
      c("character", replicate(num_conditions, "numeric"))
    counts <-
      as.matrix(
        read.table(
          path_to_count_matrix,
          header = TRUE,
          row.names = row_names,
          colClasses = colClasses
        )
      )

    print(head(counts))

    if (round) {
      counts <- round(counts)
    }


    cols <- colnames(counts)
    conditions <- gsub(opt$sample_suffix, "", cols)
    if (is.character(seq_types)) {
      seq_types <- replicate(length(cols), seq_types)
    }

    # Generate sample df
    coldata <-
      data.frame(condition = factor(conditions),
                 type = factor(seq_types))

    # Construct dds
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = coldata,
                                  design = ~condition)

    keep <- rowSums(counts(dds)) >= min_counts
    dds <- dds[keep,]

    dds$condition <- relevel(dds$condition, ref = contrast)

    return(dds)
  }


write_results_to_disk <- function(dds, contrast, outprefix) {
  for (cond in colData(dds)$condition) {
    if (cond != contrast) {
      res <-
        results(dds,
                contrast = c("condition", cond, contrast),
                parallel = TRUE)
      write.table(
        as.data.frame(res),
        file = stringr::str_interp("${outprefix}${contrast}_vs_${cond}.tsv")
        ,
        sep = "\t"
      )
    }
  }

  return(invisible())
}


plot_sample_to_sample_dist <- function(dds, contrast, outprefix) {
  plot_file <-
    stringr::str_interp("${outprefix}${contrast}_sample_to_sample_dist_heatmap.pdf")

  vsd <- vst(dds, blind = FALSE)
  sampleDists <- dist(t(assay(vsd)))

  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <-
    paste(vsd$condition, vsd$type, sep = "-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)


  hmap <- pheatmap(
    sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    col = colors
  )

  pdf(
    file = plot_file,
    onefile = TRUE,
    width = 10,
    height = 7
  )
  plot(hmap$gtable)
  dev.off()

  return(invisible())
}

plot_count_matrix <- function(dds, contrast, outprefix) {
  plot_file <-
    stringr::str_interp("${outprefix}${contrast}_count_matrix_heatmap.pdf")
  vsd <- vst(dds, blind = FALSE)
  select <- order(rowMeans(counts(dds, normalized = TRUE)),
                  decreasing = TRUE)[1:20]
  df <- as.data.frame(colData(dds)[, c("condition", "type")])

  hmap <-
    pheatmap(
      assay(vsd)[select,],
      cluster_rows = FALSE,
      show_rownames = FALSE,
      cluster_cols = FALSE,
      annotation_col = df
    )

  pdf(
    file = plot_file,
    onefile = TRUE,
    width = 10,
    height = 10
  )
  plot(hmap$gtable)
  dev.off()

  return(invisible())
}

# Parse CLI options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

missing_options <- list()

if (is.null(opt$count_matrix)) {
  missing_options <- append(missing_options, "--count_matrix")
}
if (is.null(opt$contrast)) {
  missing_options <- append(missing_options, "--contrast")
}
if (is.null(opt$outprefix)) {
  missing_options <- append(missing_options, "--outprefix")
}

opt$cpu_cores <- min(opt$cpu_cores, detectCores())
register(MulticoreParam(opt$cpu_cores))
options(Ncpus = opt$cpu_cores)
options(mc.cores = opt$cpu_cores)

if (!dir.exists(path_dir(opt$outprefix))) {
  dir.create(path_dir(opt$outprefix))
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


dds <-
  construct_dds_from_count_matrix(opt$count_matrix, opt$contrast, opt$sample_suffix)

mask <- rowSums(counts(dds)) >= opt$min_counts
dds <- dds[mask,]

# Run DESeq
dds <- DESeq(dds, parallel = TRUE)

# Write results to disk
write_results_to_disk(dds, opt$contrast, opt$outprefix)


# Make plots
plot_sample_to_sample_dist(dds, opt$contrast, opt$outprefix)
plot_count_matrix(dds, opt$contrast, opt$outprefix)
