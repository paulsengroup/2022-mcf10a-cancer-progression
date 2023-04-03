#!/usr/bin/env Rscript

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

suppressMessages(library("argparse"))
suppressMessages(library("clusterProfiler"))
suppressMessages(library("fs"))
suppressMessages(library("ggplot2"))
suppressMessages(library("org.Hs.eg.db"))
suppressMessages(library("stringr"))


import_de_genes <- function(path_to_tsv,
                            src_id) {
  df <- read.table(path_to_tsv, stringsAsFactors = F)

  if (startsWith(src_id, "ENSEMBL")) {
    # Trim version from ENSEMBL.* IDs
    ids <- gsub("\\..*$", "", row.names(df))
  } else {
    ids <- row.names(df)
  }

  df["id"] <- ids
  df <- df[!is.na(df$id),]

  return(df)
}


log_transform <- function(v) {
  v1 <- log10(v)

  # Deal with -inf/+inf
  min_ <- min(v1[is.finite(v1)])
  max_ <- max(v1[is.finite(v1)])
  v1[is.na(v1)] <- min_ - 1
  v1 <- pmax(pmin(v1, max_), min_)

  return(v1)
}


generate_gene_list <- function(df,
                               lfc_threshold,
                               pval,
                               col_name,
                               log_transform = FALSE) {
  df1 <- df[abs(df["log2FoldChange"]) >= lfc_threshold & df["padj"] <= pval,]

  if (log_transform) {
    gene_list <- log_transform(df1[, col_name])
    universe <- log_transform(df[, col_name])
  } else {
    gene_list <- df1[, col_name]
    universe <- df[, col_name]
  }

  names(gene_list) <- as.character(df1[, "id"])
  names(universe) <- as.character(df[, "id"])

  gene_list <- sort(gene_list, decreasing = TRUE)
  universe <- sort(universe, decreasing = TRUE)

  return(list(de_genes = gene_list, universe = universe))
}

enrich_go_terms <- function(genes,
                            key_type,
                            ont,
                            pval,
                            qval) {
  ego <- enrichGO(gene = names(genes$de_genes),
                  universe = names(genes$universe),
                  OrgDb = org.Hs.eg.db,
                  keyType = key_type,
                  ont = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff = pval,
                  qvalueCutoff = qval,
                  readable = TRUE)

  return(ego)
}

gene_set_enrichment_go <- function(genes,
                                   key_type,
                                   ont,
                                   pval,
                                   min_gs_size,
                                   max_gs_size) {
  stopifnot(min_gs_size < max_gs_size)
  ego <- gseGO(geneList = genes$de_genes,
               OrgDb = org.Hs.eg.db,
               keyType = key_type,
               ont = ont,
               minGSSize = min_gs_size,
               maxGSSize = max_gs_size,
               pvalueCutoff = pval,
               eps = 0)

  return(ego)
}

cluster_go_terms <- function(gene_lists,
                             labels,
                             key_type,
                             fx) {
  cluster <- list()
  for (genes in gene_lists) {
    if (fx == "gseGO") {
      x <- genes$de_genes
    } else {
      x <- list(names(genes$de_genes))
    }
    cluster <- append(cluster, x)
  }
  names(cluster) <- labels

  ck <- compareCluster(geneCluster = cluster,
                       OrgDb = org.Hs.eg.db,
                       keyType = key_type,
                       fun = fx)
  ck <- setReadable(ck,
                    OrgDb = org.Hs.eg.db,
                    keyType = key_type)

  return(ck)
}

plot_go_network <- function(ego,
                            label,
                            outprefix,
                            outsuffix,
                            layout = "sugiyama") {
  # available layouts: https://rdrr.io/cran/ggraph/man/layout_tbl_graph_igraph.html
  plot_file <- stringr::str_interp("${outprefix}_${label}_go_network_${outsuffix}.pdf")

  plt <- goplot(ego, layout = layout)
  ggsave(plot_file,
         plot = plt)

  return(invisible())
}

plot_cluster <- function(cluster,
                         outprefix,
                         outsuffix,
                         num_categories = 10,
                         font_size = 8) {
  plot_file <- stringr::str_interp("${outprefix}_go_cluster_dotplot_${outsuffix}.pdf")

  plt <- dotplot(cluster, showCategory = num_categories, font.size = font_size)
  ggsave(plot_file,
         plot = plt)

  plot_file <- stringr::str_interp("${outprefix}_go_cluster_network_${outsuffix}.pdf")

  plt <- cnetplot(cluster)
  ggsave(plot_file,
         plot = plt)

  return(invisible())
}


run_go_ora <- function(de_gene_dfs,
                       ontology,
                       outprefix,
                       identifier_type,
                       lfc_thresh,
                       pval_thresh_de,
                       pval_thresh_go,
                       qval_thresh_go) {

  gene_lists <- lapply(de_gene_dfs,
                       generate_gene_list,
                       lfc_threshold = lfc_thresh,
                       pval = pval_thresh_de,
                       col_name = "log2FoldChange")

  labels <- names(de_gene_dfs)

  for (i in 1:length(de_gene_dfs)) {
    genes <- gene_lists[[i]]
    label <- labels[[i]]

    ego <- enrich_go_terms(genes,
                           key_type = identifier_type,
                           ont = ontology,
                           pval = pval_thresh_go,
                           qval = qval_thresh_go)

    plot_go_network(ego = ego,
                    label = label,
                    outprefix = outprefix,
                    outsuffix = "ora")
  }

  go_clusters <- cluster_go_terms(gene_lists,
                                  labels,
                                  key_type = identifier_type,
                                  fx = "enrichGO")

  plot_cluster(go_clusters,
               outprefix = outprefix,
               outsuffix = "ora")

  return(invisible())
}

run_go_gse <- function(de_gene_dfs,
                       ontology,
                       outprefix,
                       identifier_type,
                       lfc_thresh,
                       pval_thresh_de,
                       pval_thresh_go,
                       min_gs_size,
                       max_gs_size) {

  gene_lists <- lapply(de_gene_dfs,
                       generate_gene_list,
                       lfc_threshold = lfc_thresh,
                       pval = pval_thresh_de,
                       col_name = "log2FoldChange")

  labels <- names(de_gene_dfs)

  for (i in 1:length(de_gene_dfs)) {
    genes <- gene_lists[[i]]
    label <- labels[[i]]

    ego <- gene_set_enrichment_go(genes,
                                  key_type = identifier_type,
                                  ont = ontology,
                                  pval = pval_thresh_go,
                                  min_gs_size = min_gs_size,
                                  max_gs_size = max_gs_size)

    plot_go_network(ego = ego,
                    label = label,
                    outprefix = outprefix,
                    outsuffix = "gse")
  }

  go_clusters <- cluster_go_terms(gene_lists,
                                  labels,
                                  key_type = identifier_type,
                                  fx = "gseGO")

  plot_cluster(go_clusters,
               outprefix = outprefix,
               outsuffix = "gse")

  return(invisible())
}


# Parse CLI options
parser <- ArgumentParser()
parser$add_argument("tsvs",
                    type = "character",
                    nargs = "+",
                    help = "Path to two or more TSVs with the list of differentially expressed genes.")
parser$add_argument("--ora",
                    action="store_true",
                    default=FALSE,
                    help="Perform analysis using enrichGO().")
parser$add_argument("--gse",
                    action="store_true",
                    default=FALSE,
                    help="Perform analysis using gseGO().")
parser$add_argument("--outprefix",
                    type = "character",
                    required = TRUE)
parser$add_argument("--identifier-type",
                    type = "character",
                    default = "ENSEMBL",
                    help = "Gene/transcript identifier type. Can be any type supported by org.Hs.eg.db.")
parser$add_argument("--labels",
                    type = "character",
                    help = "Comma-separated list of labels to use for plotting.")
parser$add_argument("--ontology",
                    type = "character",
                    required = TRUE,
                    choices = c("BP", "CC", "MF"),
                    help = "Gene ontology type.")
parser$add_argument("--lfc-thresh",
                    type = "double",
                    default = 2.5,
                    help = "Log2FoldChange threshold.")
parser$add_argument("--pval-thresh-de",
                    type = "double",
                    default = 0.01,
                    help = "Adjusted p-value threshold used to select differentially expressed genes.")
parser$add_argument("--pval-thresh-go",
                    type = "double",
                    default = 0.05,
                    help = "Adjusted p-value threshold used to select significantly enriched GO terms.")
parser$add_argument("--qval-thresh-go",
                    type = "double",
                    default = 0.05,
                    help = "Adjusted p-value threshold used to select significantly enriched GO terms.")
parser$add_argument("--min-gs-size",
                    type = "integer",
                    default = 10,
                    help = "Lower bound for the size of gene sets to be analyzed.")
parser$add_argument("--max-gs-size",
                    type = "integer",
                    default = 500,
                    help = "Upper bound for the size of gene sets to be analyzed.")
parser$add_argument("--seed",
                    type = "integer",
                    default = 891233477)

opt <- parser$parse_args()

if (opt$ora + opt$gse != 1) {
  stop("Please specify one of --ora --gse.")
}

if (!is.null(opt$labels)) {
  labels <- strsplit(opt$labels, ",")
} else {
  labels <- lapply(opt$tsvs, path_file)
  labels <- lapply(labels, path_ext_remove)
}

for (thresh in c(opt$pval_thresh_go, opt$pval_thresh_de)) {
  if (abs(thresh) > 1.0) {
    stop(
      stringr::str_interp(
        "Invalid p-value threshold: expected a number between 0.0 and 1.0, found ${thresh}"
      ),
      .call = FALSE
    )
  }
}

if (length(labels) != length(opt$tsvs)) {
  num_labels <- length(labels)
  num_tsvs <- length(opt$tsvs)

  stop(
    stringr::str_interp(
      "Expected ${num_tsv} labels, found ${num_labels}"
    ),
    .call = FALSE
  )
}

set.seed(opt$seed)

if (!dir.exists(path_dir(opt$outprefix))) {
  dir.create(path_dir(opt$outprefix))
}

# Import data
de_gene_dfs <- lapply(opt$tsvs,
                      import_de_genes,
                      src_id = opt$identifier_type)
names(de_gene_dfs) <- labels


### Run GO over-representation analysis
if (opt$ora) {
  run_go_ora(de_gene_dfs,
             opt$ontology,
             opt$outprefix,
             opt$identifier_type,
             opt$lfc_thresh,
             opt$pval_thresh_de,
             opt$pval_thresh_go,
             opt$qval_thresh_go)
}


### Run GO gene set enrichment analysis
if (opt$gse) {
  run_go_gse(de_gene_dfs,
             opt$ontology,
             opt$outprefix,
             opt$identifier_type,
             opt$lfc_thresh,
             opt$pval_thresh_de,
             opt$pval_thresh_go,
             opt$min_gs_size,
             opt$max_gs_size)
}
