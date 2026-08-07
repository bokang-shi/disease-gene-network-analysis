#!/usr/bin/env Rscript

# Apply unsupervised PRPS/RUV-III normalization to an integrated TCGA+GTEx
# SummarizedExperiment object.
#
# Input requirements:
# - an assay containing raw-ish RNA-seq counts, default: counts
# - colData columns: condition and unwanted_source
#
# Use --samples-for-prps normal for TCGA-normal/GTEx-normal validation, or
# --samples-for-prps all for exploratory full tumor+normal correction before
# downstream DEG testing.

parse_args <- function(args) {
  cfg <- list(
    input_rds = "",
    output_rds = "",
    assay_name = "counts",
    k = "1,2,3",
    samples_for_prps = "normal",
    hvg_n = 5000L,
    min_sample_for_ps = 3L,
    max_sample_for_ps = 10L,
    max_prps_sets = 3L,
    approach = "cca"
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (i == length(args)) stop(sprintf("Missing value for argument %s", key))
    value <- args[[i + 1L]]
    if (key == "--input-rds") cfg$input_rds <- value
    else if (key == "--output-rds") cfg$output_rds <- value
    else if (key == "--assay-name") cfg$assay_name <- value
    else if (key == "--k") cfg$k <- value
    else if (key == "--samples-for-prps") cfg$samples_for_prps <- value
    else if (key == "--hvg-n") cfg$hvg_n <- as.integer(value)
    else if (key == "--min-sample-for-ps") cfg$min_sample_for_ps <- as.integer(value)
    else if (key == "--max-sample-for-ps") cfg$max_sample_for_ps <- as.integer(value)
    else if (key == "--max-prps-sets") cfg$max_prps_sets <- as.integer(value)
    else if (key == "--approach") cfg$approach <- value
    else stop(sprintf("Unknown argument: %s", key))
    i <- i + 2L
  }
  if (!nzchar(cfg$input_rds)) stop("--input-rds is required")
  if (!nzchar(cfg$output_rds)) {
    cfg$output_rds <- sub("\\.rds$", "_ruvprps_unsupervised.rds", cfg$input_rds, ignore.case = TRUE)
  }
  cfg$k <- as.integer(trimws(strsplit(cfg$k, ",")[[1L]]))
  cfg
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))

for (pkg in c("SummarizedExperiment", "RUVprps", "matrixStats")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop(sprintf("Package %s is required", pkg))
}

se_in <- readRDS(cfg$input_rds)
cd <- as.data.frame(SummarizedExperiment::colData(se_in), stringsAsFactors = FALSE)
if (!all(c("condition", "unwanted_source") %in% names(cd))) {
  stop("Input object must contain condition and unwanted_source in colData")
}

counts <- SummarizedExperiment::assay(se_in, cfg$assay_name)
ruv_assay_name <- "RawCount"
hvg <- rownames(counts)
if (!is.na(cfg$hvg_n) && cfg$hvg_n > 0L && nrow(counts) > cfg$hvg_n) {
  vars <- matrixStats::rowVars(log2(as.matrix(counts) + 1))
  hvg <- rownames(counts)[order(vars, decreasing = TRUE)[seq_len(cfg$hvg_n)]]
}

if (cfg$samples_for_prps == "normal") {
  samples_to_use <- cd$condition == "normal"
} else if (cfg$samples_for_prps == "all") {
  samples_to_use <- "all"
} else {
  stop("--samples-for-prps must be normal or all")
}

source_check <- if (identical(samples_to_use, "all")) cd$unwanted_source else cd$unwanted_source[samples_to_use]
if (length(unique(source_check)) < 2L) {
  stop("Selected PRPS samples must include at least two unwanted_source levels")
}

if (cfg$approach == "cca" && !identical(samples_to_use, "all")) {
  counts <- counts[, samples_to_use, drop = FALSE]
  cd <- cd[samples_to_use, , drop = FALSE]
  samples_to_use <- "all"
}
if (cfg$approach == "cca" && length(hvg) < nrow(counts)) {
  counts <- counts[hvg, , drop = FALSE]
}

se <- RUVprps::prepareSeObj(
  data = setNames(list(counts), ruv_assay_name),
  sample.annotation = cd,
  raw.count.assay.name = ruv_assay_name,
  calculate.library.size = TRUE,
  create.gene.annotation = TRUE,
  add.housekeeping.genes = FALSE,
  verbose = TRUE
)

se <- RUVprps::createPrPsUnSupervised(
  se.obj = se,
  assay.name = ruv_assay_name,
  uv.variables = "unwanted_source",
  approach = cfg$approach,
  samples.to.use = samples_to_use,
  hvg = hvg,
  min.sample.for.ps = cfg$min_sample_for_ps,
  max.sample.for.ps = cfg$max_sample_for_ps,
  filter.prps.sets = TRUE,
  max.prps.sets = cfg$max_prps_sets,
  min.batches.to.cover = "all",
  cover.all.batches = TRUE,
  check.prps.connectedness = FALSE,
  normalization = NULL,
  apply.log = FALSE,
  apply.log.for.prps = FALSE,
  prps.group.name = "tcga_gtex_unsupervised",
  prps.sets.name = "tcga_gtex_unsupervised_prps",
  plot.output = FALSE,
  save.se.obj = TRUE,
  verbose = TRUE
)

se <- RUVprps::findNcgUnSupervised(
  se.obj = se,
  assay.name = ruv_assay_name,
  uv.variables = "unwanted_source",
  form = ~ unwanted_source,
  nb.ncg = 0.1,
  samples.to.use = samples_to_use,
  ncg.selection.method = "non.overlap",
  normalization = NULL,
  apply.log = FALSE,
  create.ncg.rank.plot = FALSE,
  plot.ncg.rank = FALSE,
  plot.ncg.assessment = FALSE,
  ncg.group.name = "tcga_gtex_unsupervised",
  ncg.set.name = "tcga_gtex_unsupervised_ncg",
  save.se.obj = TRUE,
  verbose = TRUE
)

prps_set_name <- ifelse(cfg$approach == "cca", paste("unwanted_source", "CcaPca", ruv_assay_name, sep = "|"), "tcga_gtex_unsupervised_prps")
md <- S4Vectors::metadata(se)
prps_node <- md[["PRPS"]][["un.supervised"]][["tcga_gtex_unsupervised"]][[prps_set_name]][["prps.data"]]
if (is.list(prps_node) && !is.null(dim(prps_node)) && length(prps_node) >= 1L) {
  prps_node <- as.matrix(prps_node[[1L]])
} else if (is.list(prps_node)) {
  prps_node <- lapply(prps_node, function(x) {
    if (is.matrix(x) && ncol(x) == nrow(se) && nrow(x) != nrow(se)) return(t(x))
    x
  })
} else if (is.matrix(prps_node) && ncol(prps_node) == nrow(se) && nrow(prps_node) != nrow(se)) {
  prps_node <- t(prps_node)
}
md[["PRPS"]][["un.supervised"]][["tcga_gtex_unsupervised"]][[prps_set_name]][["prps.data"]] <- prps_node
S4Vectors::metadata(se) <- md

se <- RUVprps::RUVIIIprps(
  se.obj = se,
  assay.name = ruv_assay_name,
  control.sample.types = "prps",
  prps.type = "un.supervised",
  prps.group.names = "tcga_gtex_unsupervised",
  prps.set.names = prps_set_name,
  ncg.type = "un.supervised",
  ncg.group.names = "tcga_gtex_unsupervised",
  ncg.set.names = "tcga_gtex_unsupervised_ncg",
  k = cfg$k,
  apply.log = TRUE,
  data.to.log = "assay",
  save.se.obj = TRUE,
  verbose = TRUE
)

dir.create(dirname(cfg$output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(se, cfg$output_rds)
message("Wrote unsupervised RUVprps-normalized object: ", normalizePath(cfg$output_rds, winslash = "\\", mustWork = FALSE))
