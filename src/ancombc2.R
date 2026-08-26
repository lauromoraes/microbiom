#!/usr/bin/env Rscript

# ==============================================================================
# ANCOM-BC2 differential abundance analysis from QIIME 2 exported TSV files
#
# Required inputs:
#   feature-table.tsv   BIOM/QIIME 2 count table exported as TSV
#   taxonomy.tsv        QIIME 2 taxonomy exported as TSV
#   metadata.tsv        QIIME 2 sample metadata
#
# Main outputs:
#   tables/             ANCOM-BC2 result tables and QC summaries
#   figures/            QC and differential-abundance plots (PNG + PDF)
#   objects/            RDS objects
#   logs/               parameters and sessionInfo()
#
# IMPORTANT:
#   feature-table.tsv must contain NON-RAREFIED integer counts.
# ==============================================================================


# ==============================================================================
# 1. COMMAND-LINE INTERFACE
# ==============================================================================

usage <- function(exit_status = 0) {
  cat(
"
ANCOM-BC2 analysis for QIIME 2 exported TSV files

USAGE

  Rscript ancombc2_qiime2.R [options]

REQUIRED OPTIONS

  --input-dir PATH
      Directory containing the input TSV files.

  --output-dir PATH
      Directory where results will be written.

  --group-var NAME
      Main exposure/group variable in metadata.tsv.

  --reference-level VALUE
      Reference level of --group-var.

OPTIONAL INPUT FILE NAMES

  --feature-table FILE       Default: feature-table.tsv
  --taxonomy FILE            Default: taxonomy.tsv
  --metadata FILE            Default: metadata.tsv

MODEL OPTIONS

  --covariates LIST
      Comma-separated metadata variables used as adjustment covariates.
      Example: --covariates age,bmi,hrt-use
      Default: none

  --tax-level LEVEL
      Taxonomic level used for aggregation:
      Domain, Phylum, Class, Order, Family, Genus, Species, or ASV.
      Default: Genus

  --prevalence-cutoff NUM
      Minimum prevalence used by ANCOM-BC2 (0 to <1).
      Default: 0.10

  --library-cutoff INT
      Minimum library size used by ANCOM-BC2.
      Default: 1000

  --alpha NUM
      Significance level.
      Default: 0.05

  --p-adjust METHOD
      P-value adjustment method accepted by p.adjust().
      Examples: holm, BH, bonferroni.
      Default: BH

  --s0-perc NUM
      Percentile used for variance regularization.
      Default: 0.05

  --cores INT
      Number of parallel workers.
      Default: 4

  --pseudo-sens BOOL
      Run pseudo-count sensitivity analysis.
      Default: true

  --structural-zero BOOL
      Detect structural zeros.
      Default: true

  --negative-lower-bound BOOL
      Apply negative lower-bound criterion for structural zeros.
      Default: true

MULTI-GROUP OPTIONS

  --global BOOL
      Run global test. Useful for >=3 groups.
      Default: false

  --pairwise BOOL
      Run all pairwise directional comparisons.
      Default: false

  --dunnet BOOL
      Run Dunnett-type comparisons against the reference.
      Default: false

  --trend BOOL
      Run ordered-group trend test.
      Default: false

FIGURE OPTIONS

  --top-taxa INT
      Maximum number of taxa in effect-size plots/heatmaps.
      Default: 25

  --plot-width NUM
      Figure width in inches.
      Default: 9

  --plot-height NUM
      Base figure height in inches.
      Default: 6

OTHER OPTIONS

  --seed INT
      Random seed.
      Default: 123

  --help
      Show this help message.

EXAMPLE

  Rscript ancombc2_qiime2.R \\
    --input-dir /input \\
    --output-dir /output \\
    --group-var menopausal-status \\
    --reference-level premenopause \\
    --covariates age,hrt-use \\
    --tax-level Genus \\
    --prevalence-cutoff 0.10 \\
    --library-cutoff 1000 \\
    --pseudo-sens true \\
    --cores 8

"
  )
  quit(status = exit_status)
}


parse_bool <- function(x, option_name) {
  value <- tolower(trimws(x))

  if (value %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (value %in% c("false", "f", "0", "no", "n")) return(FALSE)

  stop(
    option_name, " must be true or false; received: ", x,
    call. = FALSE
  )
}


parse_cli <- function(args) {
  defaults <- list(
    input_dir = NULL,
    output_dir = NULL,
    feature_table = "feature-table.tsv",
    taxonomy = "taxonomy.tsv",
    metadata = "metadata.tsv",
    group_var = NULL,
    reference_level = NULL,
    covariates = character(0),
    tax_level = "Genus",
    prevalence_cutoff = 0.10,
    library_cutoff = 1000,
    alpha = 0.05,
    p_adjust = "BH",
    s0_perc = 0.05,
    cores = 4L,
    pseudo_sens = TRUE,
    structural_zero = TRUE,
    negative_lower_bound = TRUE,
    global = FALSE,
    pairwise = FALSE,
    dunnet = FALSE,
    trend = FALSE,
    top_taxa = 25L,
    plot_width = 9,
    plot_height = 6,
    seed = 123L
  )

  if (length(args) == 0 || "--help" %in% args || "-h" %in% args) {
    usage(0)
  }

  option_map <- list(
    "--input-dir" = "input_dir",
    "--output-dir" = "output_dir",
    "--feature-table" = "feature_table",
    "--taxonomy" = "taxonomy",
    "--metadata" = "metadata",
    "--group-var" = "group_var",
    "--reference-level" = "reference_level",
    "--covariates" = "covariates",
    "--tax-level" = "tax_level",
    "--prevalence-cutoff" = "prevalence_cutoff",
    "--library-cutoff" = "library_cutoff",
    "--alpha" = "alpha",
    "--p-adjust" = "p_adjust",
    "--s0-perc" = "s0_perc",
    "--cores" = "cores",
    "--pseudo-sens" = "pseudo_sens",
    "--structural-zero" = "structural_zero",
    "--negative-lower-bound" = "negative_lower_bound",
    "--global" = "global",
    "--pairwise" = "pairwise",
    "--dunnet" = "dunnet",
    "--trend" = "trend",
    "--top-taxa" = "top_taxa",
    "--plot-width" = "plot_width",
    "--plot-height" = "plot_height",
    "--seed" = "seed"
  )

  out <- defaults
  i <- 1

  while (i <= length(args)) {
    arg <- args[[i]]

    # Accept both "--option value" and "--option=value".
    if (grepl("^--[^=]+=", arg)) {
      option <- sub("=.*$", "", arg)
      value <- sub("^[^=]+=", "", arg)
    } else {
      option <- arg

      if (!option %in% names(option_map)) {
        stop("Unknown option: ", option, "\nUse --help for usage.", call. = FALSE)
      }

      if (i == length(args)) {
        stop("Missing value for option: ", option, call. = FALSE)
      }

      value <- args[[i + 1]]
      i <- i + 1
    }

    if (!option %in% names(option_map)) {
      stop("Unknown option: ", option, "\nUse --help for usage.", call. = FALSE)
    }

    key <- option_map[[option]]

    if (key == "covariates") {
      value <- trimws(value)
      out[[key]] <- if (value == "" || tolower(value) %in% c("none", "null")) {
        character(0)
      } else {
        trimws(strsplit(value, ",", fixed = TRUE)[[1]])
      }
    } else if (key %in% c(
      "pseudo_sens", "structural_zero", "negative_lower_bound",
      "global", "pairwise", "dunnet", "trend"
    )) {
      out[[key]] <- parse_bool(value, option)
    } else if (key %in% c(
      "prevalence_cutoff", "alpha", "s0_perc",
      "plot_width", "plot_height"
    )) {
      out[[key]] <- as.numeric(value)
    } else if (key %in% c(
      "library_cutoff", "cores", "top_taxa", "seed"
    )) {
      out[[key]] <- as.integer(value)
    } else {
      out[[key]] <- value
    }

    i <- i + 1
  }

  required <- c("input_dir", "output_dir", "group_var", "reference_level")
  missing_required <- required[
    vapply(required, function(x) {
      is.null(out[[x]]) || length(out[[x]]) == 0 || !nzchar(out[[x]])
    }, logical(1))
  ]

  if (length(missing_required) > 0) {
    stop(
      "Missing required option(s): ",
      paste(paste0("--", gsub("_", "-", missing_required)), collapse = ", "),
      "\nUse --help for usage.",
      call. = FALSE
    )
  }

  if (!is.finite(out$prevalence_cutoff) ||
      out$prevalence_cutoff < 0 ||
      out$prevalence_cutoff >= 1) {
    stop("--prevalence-cutoff must be >= 0 and < 1.", call. = FALSE)
  }

  if (!is.finite(out$alpha) || out$alpha <= 0 || out$alpha >= 1) {
    stop("--alpha must be > 0 and < 1.", call. = FALSE)
  }

  if (!is.finite(out$s0_perc) || out$s0_perc < 0 || out$s0_perc > 1) {
    stop("--s0-perc must be between 0 and 1.", call. = FALSE)
  }

  if (is.na(out$library_cutoff) || out$library_cutoff < 0) {
    stop("--library-cutoff must be >= 0.", call. = FALSE)
  }

  if (is.na(out$cores) || out$cores < 1) {
    stop("--cores must be >= 1.", call. = FALSE)
  }

  if (is.na(out$top_taxa) || out$top_taxa < 1) {
    stop("--top-taxa must be >= 1.", call. = FALSE)
  }

  valid_tax_levels <- c(
    "Domain", "Phylum", "Class", "Order",
    "Family", "Genus", "Species", "ASV"
  )

  # Case-insensitive normalization.
  tax_match <- valid_tax_levels[
    tolower(valid_tax_levels) == tolower(out$tax_level)
  ]

  if (length(tax_match) != 1) {
    stop(
      "--tax-level must be one of: ",
      paste(valid_tax_levels, collapse = ", "),
      call. = FALSE
    )
  }

  out$tax_level <- tax_match
  out
}


opt <- parse_cli(commandArgs(trailingOnly = TRUE))


# ==============================================================================
# 2. PACKAGES
# ==============================================================================

required_packages <- c(
  "ANCOMBC",
  "phyloseq",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running the script.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(ANCOMBC)
  library(phyloseq)
  library(ggplot2)
})


# ==============================================================================
# 3. PATHS AND OUTPUT DIRECTORIES
# ==============================================================================

input_dir <- normalizePath(opt$input_dir, mustWork = TRUE)

feature_table_file <- file.path(input_dir, opt$feature_table)
taxonomy_file <- file.path(input_dir, opt$taxonomy)
metadata_file <- file.path(input_dir, opt$metadata)

output_dir <- opt$output_dir
tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
objects_dir <- file.path(output_dir, "objects")
logs_dir <- file.path(output_dir, "logs")

for (d in c(output_dir, tables_dir, figures_dir, objects_dir, logs_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# ==============================================================================
# 4. HELPERS
# ==============================================================================

section <- function(text) {
  cat("\n", paste(rep("=", 78), collapse = ""), "\n", sep = "")
  cat(text, "\n")
  cat(paste(rep("=", 78), collapse = ""), "\n", sep = "")
}


check_file <- function(path) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path, call. = FALSE)
  }
}


safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}


save_plot_both <- function(
  plot,
  basename,
  width = opt$plot_width,
  height = opt$plot_height
) {
  png_path <- file.path(figures_dir, paste0(basename, ".png"))
  pdf_path <- file.path(figures_dir, paste0(basename, ".pdf"))

  ggsave(
    filename = png_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  ggsave(
    filename = pdf_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )
}


write_tsv <- function(x, path, row_names = FALSE) {
  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = row_names,
    col.names = TRUE,
    na = "NA"
  )
}


write_if_present <- function(x, filename) {
  if (!is.null(x)) {
    write_tsv(x, file.path(tables_dir, filename))
  }
}


# ==============================================================================
# 5. READ QIIME 2 TSV FILES
# ==============================================================================

read_feature_table <- function(path) {
  first_line <- readLines(path, n = 1, warn = FALSE)

  # BIOM TSV usually starts with "# Constructed from biom file".
  skip_n <- if (
    length(first_line) == 1 &&
    grepl("^# Constructed from biom file", first_line)
  ) 1 else 0

  tab <- read.delim(
    path,
    sep = "\t",
    header = TRUE,
    skip = skip_n,
    check.names = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )

  if (ncol(tab) < 2) {
    stop("feature-table.tsv is not a valid feature x sample table.", call. = FALSE)
  }

  feature_ids <- as.character(tab[[1]])

  if (anyNA(feature_ids) || any(feature_ids == "")) {
    stop("Missing feature IDs in feature-table.tsv.", call. = FALSE)
  }

  if (anyDuplicated(feature_ids)) {
    stop("Duplicated feature IDs in feature-table.tsv.", call. = FALSE)
  }

  tab[[1]] <- NULL

  count_matrix <- as.matrix(tab)
  suppressWarnings(storage.mode(count_matrix) <- "numeric")
  rownames(count_matrix) <- feature_ids

  if (anyNA(count_matrix)) {
    stop(
      "NA or non-numeric values detected in feature-table.tsv.",
      call. = FALSE
    )
  }

  if (any(count_matrix < 0)) {
    stop("Negative counts detected in feature-table.tsv.", call. = FALSE)
  }

  integer_like <- abs(count_matrix - round(count_matrix)) <=
    .Machine$double.eps^0.5

  if (!all(integer_like)) {
    warning(
      "Non-integer values were detected in feature-table.tsv. ",
      "ANCOM-BC2 should receive raw count data, not relative abundances ",
      "or transformed values."
    )
  }

  count_matrix
}


read_metadata <- function(path) {
  meta <- read.delim(
    path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN")
  )

  if (nrow(meta) == 0 || ncol(meta) < 2) {
    stop("metadata.tsv appears empty or invalid.", call. = FALSE)
  }

  id_candidates <- c(
    "#SampleID", "#Sample ID", "sample-id", "sample_id",
    "SampleID", "sampleid", "ID", "id"
  )

  id_col <- id_candidates[id_candidates %in% names(meta)][1]

  if (is.na(id_col)) {
    id_col <- names(meta)[1]
    warning(
      "No standard QIIME 2 sample-ID column found. ",
      "Using the first metadata column: ", id_col
    )
  }

  # QIIME 2 metadata can contain a "#q2:types" row.
  q2_type_row <- as.character(meta[[id_col]]) == "#q2:types"
  q2_type_row[is.na(q2_type_row)] <- FALSE

  if (any(q2_type_row)) {
    meta <- meta[!q2_type_row, , drop = FALSE]
  }

  sample_ids <- trimws(as.character(meta[[id_col]]))

  if (anyNA(sample_ids) || any(sample_ids == "")) {
    stop("Missing sample IDs in metadata.tsv.", call. = FALSE)
  }

  if (anyDuplicated(sample_ids)) {
    stop("Duplicated sample IDs in metadata.tsv.", call. = FALSE)
  }

  meta[[id_col]] <- NULL
  rownames(meta) <- sample_ids

  original_names <- names(meta)
  r_names <- make.names(original_names, unique = TRUE)
  names(meta) <- r_names

  name_map <- data.frame(
    original_name = original_names,
    r_name = r_names,
    stringsAsFactors = FALSE
  )

  attr(meta, "name_map") <- name_map
  meta
}


read_taxonomy <- function(path) {
  tax <- read.delim(
    path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )

  feature_candidates <- c(
    "Feature ID", "Feature.ID", "#OTU ID", "OTU ID",
    "feature-id", "feature_id"
  )

  taxonomy_candidates <- c(
    "Taxon", "Taxonomy", "taxonomy", "taxon"
  )

  feature_col <- feature_candidates[feature_candidates %in% names(tax)][1]
  tax_col <- taxonomy_candidates[taxonomy_candidates %in% names(tax)][1]

  if (is.na(feature_col)) {
    feature_col <- names(tax)[1]
    warning(
      "Taxonomy feature-ID column not recognized. ",
      "Using the first column: ", feature_col
    )
  }

  if (is.na(tax_col)) {
    stop(
      "Could not identify Taxon/Taxonomy column in taxonomy.tsv.",
      call. = FALSE
    )
  }

  feature_ids <- as.character(tax[[feature_col]])
  tax_strings <- as.character(tax[[tax_col]])

  if (anyDuplicated(feature_ids)) {
    stop("Duplicated feature IDs in taxonomy.tsv.", call. = FALSE)
  }

  ranks <- c(
    "Domain", "Phylum", "Class", "Order",
    "Family", "Genus", "Species"
  )

  tax_matrix <- matrix(
    NA_character_,
    nrow = length(feature_ids),
    ncol = length(ranks),
    dimnames = list(feature_ids, ranks)
  )

  normalize_taxon <- function(x) {
    x <- trimws(x)

    # Typical SILVA/Greengenes prefixes: d__, k__, p__, c__, o__, f__, g__, s__
    x <- sub("^[A-Za-z]__", "", x)
    x <- sub("^[A-Za-z]_[0-9]+__", "", x)

    if (x %in% c(
      "", "NA", "NaN", "Unassigned",
      "uncultured", "unidentified"
    )) {
      return(NA_character_)
    }

    x
  }

  split_tax <- strsplit(tax_strings, ";", fixed = TRUE)

  for (i in seq_along(split_tax)) {
    values <- vapply(split_tax[[i]], normalize_taxon, character(1))
    n <- min(length(values), length(ranks))

    if (n > 0) {
      tax_matrix[i, seq_len(n)] <- values[seq_len(n)]
    }
  }

  tax_matrix
}


resolve_metadata_name <- function(name, name_map) {
  hit_original <- name_map$r_name[name_map$original_name == name]

  if (length(hit_original) == 1) {
    return(hit_original)
  }

  hit_r <- name_map$r_name[name_map$r_name == name]

  if (length(hit_r) == 1) {
    return(hit_r)
  }

  stop(
    "Metadata variable not found: ", name,
    "\nAvailable variables: ",
    paste(name_map$original_name, collapse = ", "),
    call. = FALSE
  )
}


# ==============================================================================
# 6. LOAD AND ALIGN DATA
# ==============================================================================

section("Reading input files")

check_file(feature_table_file)
check_file(taxonomy_file)
check_file(metadata_file)

counts <- read_feature_table(feature_table_file)
taxonomy <- read_taxonomy(taxonomy_file)
metadata <- read_metadata(metadata_file)

name_map <- attr(metadata, "name_map")
write_tsv(name_map, file.path(tables_dir, "metadata-column-name-map.tsv"))

cat("Feature table: ", nrow(counts), " features x ",
    ncol(counts), " samples\n", sep = "")
cat("Taxonomy:      ", nrow(taxonomy), " features\n", sep = "")
cat("Metadata:      ", nrow(metadata), " samples x ",
    ncol(metadata), " variables\n", sep = "")


section("Matching feature IDs")

common_features <- intersect(rownames(counts), rownames(taxonomy))

if (length(common_features) == 0) {
  stop(
    "No shared feature IDs between feature-table.tsv and taxonomy.tsv.",
    call. = FALSE
  )
}

features_without_taxonomy <- setdiff(rownames(counts), rownames(taxonomy))
taxonomy_without_counts <- setdiff(rownames(taxonomy), rownames(counts))

if (length(features_without_taxonomy) > 0) {
  warning(
    length(features_without_taxonomy),
    " count-table features have no taxonomy and will be removed."
  )
}

counts <- counts[common_features, , drop = FALSE]
taxonomy <- taxonomy[common_features, , drop = FALSE]


section("Matching sample IDs")

common_samples <- intersect(colnames(counts), rownames(metadata))

if (length(common_samples) == 0) {
  stop(
    "No shared sample IDs between feature-table.tsv and metadata.tsv.",
    call. = FALSE
  )
}

samples_only_counts <- setdiff(colnames(counts), rownames(metadata))
samples_only_metadata <- setdiff(rownames(metadata), colnames(counts))

if (length(samples_only_counts) > 0) {
  warning(
    length(samples_only_counts),
    " sample(s) occur only in feature-table.tsv and will be removed."
  )
}

if (length(samples_only_metadata) > 0) {
  warning(
    length(samples_only_metadata),
    " sample(s) occur only in metadata.tsv and will be removed."
  )
}

# Preserve count-table sample order.
counts <- counts[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

stopifnot(identical(colnames(counts), rownames(metadata)))


# ==============================================================================
# 7. PREPARE THE STATISTICAL MODEL
# ==============================================================================

section("Preparing the statistical model")

group_var <- resolve_metadata_name(opt$group_var, name_map)

covariates <- if (length(opt$covariates) > 0) {
  vapply(
    opt$covariates,
    resolve_metadata_name,
    character(1),
    name_map = name_map
  )
} else {
  character(0)
}

model_variables <- unique(c(group_var, covariates))

# Remove samples with missing data in model variables.
complete_idx <- complete.cases(metadata[, model_variables, drop = FALSE])

if (!all(complete_idx)) {
  warning(
    sum(!complete_idx),
    " sample(s) have missing values in model variables and will be removed."
  )

  metadata <- metadata[complete_idx, , drop = FALSE]
  counts <- counts[, rownames(metadata), drop = FALSE]
}

metadata[[group_var]] <- factor(metadata[[group_var]])

if (!opt$reference_level %in% levels(metadata[[group_var]])) {
  stop(
    "Reference level '", opt$reference_level,
    "' is not present in ", opt$group_var,
    ". Available levels: ",
    paste(levels(metadata[[group_var]]), collapse = ", "),
    call. = FALSE
  )
}

metadata[[group_var]] <- stats::relevel(
  metadata[[group_var]],
  ref = opt$reference_level
)

metadata[[group_var]] <- droplevels(metadata[[group_var]])

if (nlevels(metadata[[group_var]]) < 2) {
  stop(
    "The group variable contains fewer than two levels after filtering.",
    call. = FALSE
  )
}

fix_formula <- paste(c(group_var, covariates), collapse = " + ")

cat("Group variable (input): ", opt$group_var, "\n", sep = "")
cat("Group variable (R):     ", group_var, "\n", sep = "")
cat("Levels:                 ",
    paste(levels(metadata[[group_var]]), collapse = ", "), "\n", sep = "")
cat("Reference:              ", levels(metadata[[group_var]])[1], "\n", sep = "")
cat("Fixed-effects formula:  ", fix_formula, "\n", sep = "")

if (nlevels(metadata[[group_var]]) < 3 &&
    any(c(opt$global, opt$pairwise, opt$dunnet, opt$trend))) {
  warning(
    "Multi-group tests were requested but ", opt$group_var,
    " has only ", nlevels(metadata[[group_var]]), " levels. ",
    "These procedures are generally intended for multi-group analyses."
  )
}


# ==============================================================================
# 8. BUILD PHYLOSEQ
# ==============================================================================

section("Building phyloseq object")

ps <- phyloseq(
  otu_table(counts, taxa_are_rows = TRUE),
  tax_table(taxonomy),
  sample_data(metadata)
)

print(ps)


# ==============================================================================
# 9. INPUT QC TABLES AND FIGURES
# ==============================================================================

section("Generating input QC")

library_sizes <- sample_sums(ps)

qc_samples <- data.frame(
  sample_id = names(library_sizes),
  library_size = as.numeric(library_sizes),
  group = as.character(sample_data(ps)[[group_var]]),
  stringsAsFactors = FALSE
)

write_tsv(
  qc_samples,
  file.path(tables_dir, "sample-library-sizes.tsv")
)

cat("Library-size summary:\n")
print(summary(library_sizes))

if (all(library_sizes < opt$library_cutoff)) {
  stop(
    "All samples are below --library-cutoff=",
    opt$library_cutoff,
    ". Review the threshold.",
    call. = FALSE
  )
}

# 9.1 Library size by sample
qc_samples$sample_id <- factor(
  qc_samples$sample_id,
  levels = qc_samples$sample_id[order(qc_samples$library_size)]
)

p_library <- ggplot(
  qc_samples,
  aes(x = sample_id, y = library_size, fill = group)
) +
  geom_col() +
  geom_hline(
    yintercept = opt$library_cutoff,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  coord_flip() +
  labs(
    title = "Sequencing depth by sample",
    subtitle = paste0(
      "Dashed line: library cutoff = ", opt$library_cutoff
    ),
    x = "Sample",
    y = "Library size (reads)",
    fill = opt$group_var
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank()
  )

library_height <- max(
  opt$plot_height,
  min(18, 3 + 0.18 * nrow(qc_samples))
)

save_plot_both(
  p_library,
  "qc-library-size-by-sample",
  height = library_height
)


# 9.2 Library-size distribution by group
p_library_group <- ggplot(
  qc_samples,
  aes(x = group, y = library_size, fill = group)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.7, size = 2) +
  geom_hline(
    yintercept = opt$library_cutoff,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  labs(
    title = "Sequencing depth by group",
    x = opt$group_var,
    y = "Library size (reads)"
  ) +
  theme_bw() +
  theme(legend.position = "none")

save_plot_both(
  p_library_group,
  "qc-library-size-by-group"
)


# 9.3 Feature prevalence
otu_mat <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)

feature_prevalence <- rowMeans(otu_mat > 0)

qc_prevalence <- data.frame(
  feature_id = names(feature_prevalence),
  prevalence = as.numeric(feature_prevalence),
  pass_cutoff = feature_prevalence >= opt$prevalence_cutoff,
  stringsAsFactors = FALSE
)

write_tsv(
  qc_prevalence,
  file.path(tables_dir, "feature-prevalence.tsv")
)

p_prevalence <- ggplot(
  qc_prevalence,
  aes(x = prevalence)
) +
  geom_histogram(
    bins = 30,
    boundary = 0,
    closed = "left"
  ) +
  geom_vline(
    xintercept = opt$prevalence_cutoff,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  labs(
    title = "Feature prevalence before ANCOM-BC2 filtering",
    subtitle = paste0(
      "Dashed line: prevalence cutoff = ",
      opt$prevalence_cutoff
    ),
    x = "Proportion of samples with non-zero counts",
    y = "Number of features"
  ) +
  theme_bw()

save_plot_both(
  p_prevalence,
  "qc-feature-prevalence"
)


# 9.4 Zero fraction by sample
zero_fraction <- colMeans(otu_mat == 0)

qc_zeros <- data.frame(
  sample_id = names(zero_fraction),
  zero_fraction = as.numeric(zero_fraction),
  group = as.character(sample_data(ps)[[group_var]]),
  stringsAsFactors = FALSE
)

write_tsv(
  qc_zeros,
  file.path(tables_dir, "sample-zero-fraction.tsv")
)

p_zeros <- ggplot(
  qc_zeros,
  aes(x = group, y = zero_fraction, fill = group)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
  labs(
    title = "Fraction of zero counts by sample",
    x = opt$group_var,
    y = "Zero-count fraction"
  ) +
  theme_bw() +
  theme(legend.position = "none")

save_plot_both(
  p_zeros,
  "qc-zero-fraction-by-group"
)


# ==============================================================================
# 10. RUN ANCOM-BC2
# ==============================================================================

section("Running ANCOM-BC2")

set.seed(opt$seed)

tax_level_arg <- if (opt$tax_level == "ASV") NULL else opt$tax_level

if (!is.null(tax_level_arg) && !tax_level_arg %in% rank_names(ps)) {
  stop(
    "Taxonomic level '", tax_level_arg,
    "' is unavailable. Available ranks: ",
    paste(rank_names(ps), collapse = ", "),
    call. = FALSE
  )
}

output <- ancombc2(
  data = ps,
  assay_name = "counts",
  tax_level = tax_level_arg,
  fix_formula = fix_formula,
  rand_formula = NULL,
  p_adj_method = opt$p_adjust,
  pseudo_sens = opt$pseudo_sens,
  prv_cut = opt$prevalence_cutoff,
  lib_cut = opt$library_cutoff,
  s0_perc = opt$s0_perc,
  group = group_var,
  struc_zero = opt$structural_zero,
  neg_lb = opt$negative_lower_bound,
  alpha = opt$alpha,
  n_cl = opt$cores,
  verbose = TRUE,
  global = opt$global,
  pairwise = opt$pairwise,
  dunnet = opt$dunnet,
  trend = opt$trend
)


# ==============================================================================
# 11. EXPORT ANCOM-BC2 TABLES
# ==============================================================================

section("Exporting ANCOM-BC2 tables")

write_if_present(output$res, "ancombc2-primary-results.tsv")
write_if_present(output$res_global, "ancombc2-global-results.tsv")
write_if_present(output$res_pair, "ancombc2-pairwise-results.tsv")
write_if_present(output$res_dunn, "ancombc2-dunnett-results.tsv")
write_if_present(output$res_trend, "ancombc2-trend-results.tsv")
write_if_present(output$zero_ind, "ancombc2-structural-zeros.tsv")

if (!is.null(output$bias_correct_log_table)) {
  bias_table <- data.frame(
    taxon = rownames(output$bias_correct_log_table),
    output$bias_correct_log_table,
    check.names = FALSE
  )

  write_tsv(
    bias_table,
    file.path(tables_dir, "ancombc2-bias-corrected-log-abundance.tsv")
  )
}


# ==============================================================================
# 12. PRIMARY-RESULT PLOTS
# ==============================================================================

section("Generating ANCOM-BC2 result figures")

res <- output$res

if (is.null(res) || nrow(res) == 0) {
  warning("ANCOM-BC2 returned no primary result rows; result plots skipped.")
} else {

  # ANCOM-BC2 uses the same coefficient suffix for:
  # lfc_, se_, W_, p_, q_, diff_, and passed_ss_.
  lfc_cols <- grep("^lfc_", names(res), value = TRUE)

  # The main group contrasts generally contain the sanitized group variable
  # in their coefficient name. If this is not detectable, retain all
  # non-intercept fixed-effect LFC columns.
  group_lfc_cols <- lfc_cols[
    grepl(group_var, lfc_cols, fixed = TRUE)
  ]

  group_lfc_cols <- setdiff(
    group_lfc_cols,
    c("lfc_(Intercept)", "lfc_Intercept")
  )

  if (length(group_lfc_cols) == 0) {
    warning(
      "Could not identify group-specific LFC columns by name. ",
      "Using all non-intercept LFC columns for result figures."
    )

    group_lfc_cols <- setdiff(
      lfc_cols,
      c("lfc_(Intercept)", "lfc_Intercept")
    )
  }

  all_plot_rows <- list()

  for (lfc_col in group_lfc_cols) {

    suffix <- sub("^lfc_", "", lfc_col)
    q_col <- paste0("q_", suffix)
    p_col <- paste0("p_", suffix)
    diff_col <- paste0("diff_", suffix)
    ss_col <- paste0("passed_ss_", suffix)

    if (!q_col %in% names(res)) {
      warning(
        "Skipping coefficient '", suffix,
        "' because ", q_col, " is absent."
      )
      next
    }

    plot_df <- data.frame(
      taxon = as.character(res$taxon),
      lfc = suppressWarnings(as.numeric(res[[lfc_col]])),
      q = suppressWarnings(as.numeric(res[[q_col]])),
      stringsAsFactors = FALSE
    )

    plot_df$p <- if (p_col %in% names(res)) {
      suppressWarnings(as.numeric(res[[p_col]]))
    } else {
      NA_real_
    }

    plot_df$diff <- if (diff_col %in% names(res)) {
      as.logical(res[[diff_col]])
    } else {
      !is.na(plot_df$q) & plot_df$q <= opt$alpha
    }

    plot_df$passed_ss <- if (ss_col %in% names(res)) {
      as.logical(res[[ss_col]])
    } else {
      NA
    }

    plot_df$significant <- !is.na(plot_df$q) &
      plot_df$q <= opt$alpha &
      !is.na(plot_df$diff) &
      plot_df$diff

    if (opt$pseudo_sens && ss_col %in% names(res)) {
      plot_df$robust_significant <- plot_df$significant &
        !is.na(plot_df$passed_ss) &
        plot_df$passed_ss
    } else {
      plot_df$robust_significant <- plot_df$significant
    }

    plot_df$minus_log10_q <- -log10(
      pmax(plot_df$q, .Machine$double.xmin)
    )

    plot_df$contrast <- suffix

    all_plot_rows[[suffix]] <- plot_df

    contrast_file <- safe_filename(suffix)

    # --------------------------------------------------------------------------
    # 12.1 Volcano plot
    # --------------------------------------------------------------------------

    volcano_df <- plot_df[
      is.finite(plot_df$lfc) & is.finite(plot_df$minus_log10_q),
      ,
      drop = FALSE
    ]

    if (nrow(volcano_df) > 0) {
      volcano_df$status <- ifelse(
        volcano_df$robust_significant,
        "Significant",
        "Not significant"
      )

      p_volcano <- ggplot(
        volcano_df,
        aes(
          x = lfc,
          y = minus_log10_q,
          shape = status
        )
      ) +
        geom_point(alpha = 0.75, size = 2.3) +
        geom_vline(
          xintercept = 0,
          linetype = "dotted",
          linewidth = 0.5
        ) +
        geom_hline(
          yintercept = -log10(opt$alpha),
          linetype = "dashed",
          linewidth = 0.5
        ) +
        labs(
          title = paste0("ANCOM-BC2 volcano plot: ", suffix),
          subtitle = paste0(
            "Significance threshold: adjusted p-value ≤ ",
            opt$alpha,
            if (opt$pseudo_sens && ss_col %in% names(res)) {
              " and passed sensitivity analysis"
            } else {
              ""
            }
          ),
          x = "Log fold change",
          y = expression(-log[10]("adjusted p-value")),
          shape = NULL
        ) +
        theme_bw() +
        theme(legend.position = "bottom")

      save_plot_both(
        p_volcano,
        paste0("ancombc2-volcano-", contrast_file)
      )
    }


    # --------------------------------------------------------------------------
    # 12.2 Effect-size plot
    #
    # Prefer robust significant taxa. If none exist, show taxa with the
    # smallest q-values so the figure remains diagnostically useful.
    # --------------------------------------------------------------------------

    effect_df <- plot_df[
      is.finite(plot_df$lfc) & !is.na(plot_df$q),
      ,
      drop = FALSE
    ]

    sig_effect_df <- effect_df[
      effect_df$robust_significant,
      ,
      drop = FALSE
    ]

    if (nrow(sig_effect_df) > 0) {
      effect_df_show <- sig_effect_df[
        order(abs(sig_effect_df$lfc), decreasing = TRUE),
        ,
        drop = FALSE
      ]
      effect_subtitle <- paste0(
        "Differentially abundant taxa; q ≤ ", opt$alpha
      )
    } else {
      effect_df_show <- effect_df[
        order(effect_df$q, -abs(effect_df$lfc)),
        ,
        drop = FALSE
      ]
      effect_subtitle <- paste0(
        "No robust significant taxa; showing taxa with the smallest q-values"
      )
    }

    effect_df_show <- head(effect_df_show, opt$top_taxa)

    if (nrow(effect_df_show) > 0) {
      effect_df_show$taxon <- factor(
        effect_df_show$taxon,
        levels = effect_df_show$taxon[
          order(effect_df_show$lfc)
        ]
      )

      effect_df_show$direction <- ifelse(
        effect_df_show$lfc >= 0,
        "Positive LFC",
        "Negative LFC"
      )

      p_effect <- ggplot(
        effect_df_show,
        aes(x = taxon, y = lfc, fill = direction)
      ) +
        geom_col() +
        geom_hline(
          yintercept = 0,
          linewidth = 0.5
        ) +
        coord_flip() +
        labs(
          title = paste0("ANCOM-BC2 effect sizes: ", suffix),
          subtitle = effect_subtitle,
          x = NULL,
          y = "Log fold change",
          fill = NULL
        ) +
        theme_bw() +
        theme(
          legend.position = "bottom",
          panel.grid.major.y = element_blank()
        )

      effect_height <- max(
        opt$plot_height,
        min(14, 3 + 0.3 * nrow(effect_df_show))
      )

      save_plot_both(
        p_effect,
        paste0("ancombc2-effect-size-", contrast_file),
        height = effect_height
      )
    }


    # --------------------------------------------------------------------------
    # 12.3 q-value vs effect magnitude
    # --------------------------------------------------------------------------

    q_effect_df <- effect_df[
      is.finite(effect_df$lfc) & is.finite(effect_df$q),
      ,
      drop = FALSE
    ]

    if (nrow(q_effect_df) > 0) {
      q_effect_df$status <- ifelse(
        q_effect_df$robust_significant,
        "Significant",
        "Not significant"
      )

      p_q_lfc <- ggplot(
        q_effect_df,
        aes(
          x = abs(lfc),
          y = q,
          shape = status
        )
      ) +
        geom_point(alpha = 0.75, size = 2.3) +
        geom_hline(
          yintercept = opt$alpha,
          linetype = "dashed",
          linewidth = 0.5
        ) +
        scale_y_continuous(
          trans = "reverse",
          limits = c(1, 0)
        ) +
        labs(
          title = paste0("Adjusted p-value vs effect magnitude: ", suffix),
          x = "Absolute log fold change",
          y = "Adjusted p-value",
          shape = NULL
        ) +
        theme_bw() +
        theme(legend.position = "bottom")

      save_plot_both(
        p_q_lfc,
        paste0("ancombc2-qvalue-vs-effect-", contrast_file)
      )
    }
  }


  # ---------------------------------------------------------------------------
  # 12.4 Heatmap-style plot across multiple group contrasts
  # ---------------------------------------------------------------------------

  if (length(all_plot_rows) > 1) {
    long_df <- do.call(
      rbind,
      lapply(names(all_plot_rows), function(nm) {
        x <- all_plot_rows[[nm]]
        x[, c(
          "taxon", "lfc", "q", "robust_significant", "contrast"
        )]
      })
    )

    # Select taxa with significant findings first; otherwise select strongest
    # maximum absolute effects.
    taxon_summary <- aggregate(
      cbind(
        max_abs_lfc = abs(long_df$lfc),
        any_sig = as.numeric(long_df$robust_significant)
      ),
      by = list(taxon = long_df$taxon),
      FUN = max,
      na.rm = TRUE
    )

    taxon_summary <- taxon_summary[
      order(
        -taxon_summary$any_sig,
        -taxon_summary$max_abs_lfc
      ),
      ,
      drop = FALSE
    ]

    selected_taxa <- head(
      taxon_summary$taxon,
      opt$top_taxa
    )

    heat_df <- long_df[
      long_df$taxon %in% selected_taxa,
      ,
      drop = FALSE
    ]

    if (nrow(heat_df) > 0) {
      heat_df$taxon <- factor(
        heat_df$taxon,
        levels = rev(selected_taxa)
      )

      p_heat <- ggplot(
        heat_df,
        aes(x = contrast, y = taxon, fill = lfc)
      ) +
        geom_tile() +
        geom_point(
          data = heat_df[heat_df$robust_significant, , drop = FALSE],
          shape = 8,
          size = 2.8
        ) +
        scale_fill_gradient2(
          midpoint = 0
        ) +
        labs(
          title = "ANCOM-BC2 log fold changes across group contrasts",
          subtitle = paste0(
            "Asterisks indicate robust significant results (q ≤ ",
            opt$alpha, ")"
          ),
          x = "Contrast",
          y = NULL,
          fill = "LFC"
        ) +
        theme_bw() +
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1
          ),
          panel.grid = element_blank()
        )

      heat_height <- max(
        opt$plot_height,
        min(14, 3 + 0.3 * length(selected_taxa))
      )

      save_plot_both(
        p_heat,
        "ancombc2-lfc-heatmap",
        height = heat_height
      )
    }
  }
}


# ==============================================================================
# 13. GLOBAL / PAIRWISE / DUNNET RESULT HEATMAPS
# ==============================================================================

plot_lfc_result_table <- function(result_table, plot_prefix, plot_title) {
  if (is.null(result_table) || nrow(result_table) == 0) {
    return(invisible(NULL))
  }

  lfc_cols <- grep("^lfc_", names(result_table), value = TRUE)

  if (length(lfc_cols) == 0) {
    return(invisible(NULL))
  }

  rows <- lapply(lfc_cols, function(lfc_col) {
    suffix <- sub("^lfc_", "", lfc_col)
    diff_col <- paste0("diff_", suffix)
    ss_col <- paste0("passed_ss_", suffix)

    data.frame(
      taxon = as.character(result_table$taxon),
      contrast = suffix,
      lfc = suppressWarnings(as.numeric(result_table[[lfc_col]])),
      diff = if (diff_col %in% names(result_table)) {
        as.logical(result_table[[diff_col]])
      } else {
        FALSE
      },
      passed_ss = if (ss_col %in% names(result_table)) {
        as.logical(result_table[[ss_col]])
      } else {
        NA
      },
      stringsAsFactors = FALSE
    )
  })

  long_df <- do.call(rbind, rows)
  long_df <- long_df[is.finite(long_df$lfc), , drop = FALSE]

  if (nrow(long_df) == 0) return(invisible(NULL))

  long_df$robust_significant <- long_df$diff

  if (opt$pseudo_sens && any(!is.na(long_df$passed_ss))) {
    long_df$robust_significant <- long_df$diff &
      (is.na(long_df$passed_ss) | long_df$passed_ss)
  }

  taxon_summary <- aggregate(
    cbind(
      max_abs_lfc = abs(long_df$lfc),
      any_sig = as.numeric(long_df$robust_significant)
    ),
    by = list(taxon = long_df$taxon),
    FUN = max,
    na.rm = TRUE
  )

  taxon_summary <- taxon_summary[
    order(-taxon_summary$any_sig, -taxon_summary$max_abs_lfc),
    ,
    drop = FALSE
  ]

  selected_taxa <- head(taxon_summary$taxon, opt$top_taxa)

  long_df <- long_df[
    long_df$taxon %in% selected_taxa,
    ,
    drop = FALSE
  ]

  long_df$taxon <- factor(
    long_df$taxon,
    levels = rev(selected_taxa)
  )

  p <- ggplot(
    long_df,
    aes(x = contrast, y = taxon, fill = lfc)
  ) +
    geom_tile() +
    geom_point(
      data = long_df[long_df$robust_significant, , drop = FALSE],
      shape = 8,
      size = 2.7
    ) +
    scale_fill_gradient2(midpoint = 0) +
    labs(
      title = plot_title,
      subtitle = "Asterisks indicate differential-abundance calls",
      x = "Contrast",
      y = NULL,
      fill = "LFC"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )

  height <- max(
    opt$plot_height,
    min(14, 3 + 0.3 * length(selected_taxa))
  )

  save_plot_both(
    p,
    plot_prefix,
    height = height
  )

  invisible(NULL)
}


if (opt$pairwise) {
  plot_lfc_result_table(
    output$res_pair,
    "ancombc2-pairwise-lfc-heatmap",
    "ANCOM-BC2 pairwise log fold changes"
  )
}

if (opt$dunnet) {
  plot_lfc_result_table(
    output$res_dunn,
    "ancombc2-dunnett-lfc-heatmap",
    "ANCOM-BC2 Dunnett-type log fold changes"
  )
}


# ==============================================================================
# 14. STRUCTURAL-ZERO FIGURE
# ==============================================================================

if (!is.null(output$zero_ind) && nrow(output$zero_ind) > 0) {
  zero_df <- output$zero_ind

  # zero_ind generally contains the taxon column plus logical group columns.
  group_zero_cols <- setdiff(names(zero_df), "taxon")

  if (length(group_zero_cols) > 0) {
    zero_long <- do.call(
      rbind,
      lapply(group_zero_cols, function(col) {
        data.frame(
          taxon = as.character(zero_df$taxon),
          group = col,
          structural_zero = as.logical(zero_df[[col]]),
          stringsAsFactors = FALSE
        )
      })
    )

    zero_hits <- zero_long[
      !is.na(zero_long$structural_zero) &
      zero_long$structural_zero,
      ,
      drop = FALSE
    ]

    if (nrow(zero_hits) > 0) {
      p_zero <- ggplot(
        zero_hits,
        aes(x = group, y = taxon)
      ) +
        geom_point(shape = 4, size = 3, stroke = 1) +
        labs(
          title = "Structural zeros detected by ANCOM-BC2",
          x = opt$group_var,
          y = NULL
        ) +
        theme_bw()

      zero_height <- max(
        opt$plot_height,
        min(14, 3 + 0.3 * length(unique(zero_hits$taxon)))
      )

      save_plot_both(
        p_zero,
        "ancombc2-structural-zeros",
        height = zero_height
      )
    }
  }
}


# ==============================================================================
# 15. SAVE OBJECTS, PARAMETERS, AND SESSION INFO
# ==============================================================================

section("Saving reproducibility information")

saveRDS(
  output,
  file.path(objects_dir, "ancombc2-output.rds")
)

saveRDS(
  ps,
  file.path(objects_dir, "phyloseq-input.rds")
)

parameter_table <- data.frame(
  parameter = c(
    "input_dir",
    "output_dir",
    "feature_table",
    "taxonomy",
    "metadata",
    "group_var_original",
    "group_var_R",
    "reference_level",
    "covariates_original",
    "covariates_R",
    "fix_formula",
    "tax_level",
    "prevalence_cutoff",
    "library_cutoff",
    "alpha",
    "p_adjust",
    "s0_perc",
    "cores",
    "pseudo_sens",
    "structural_zero",
    "negative_lower_bound",
    "global",
    "pairwise",
    "dunnet",
    "trend",
    "top_taxa",
    "seed"
  ),
  value = c(
    input_dir,
    output_dir,
    opt$feature_table,
    opt$taxonomy,
    opt$metadata,
    opt$group_var,
    group_var,
    opt$reference_level,
    paste(opt$covariates, collapse = ","),
    paste(covariates, collapse = ","),
    fix_formula,
    opt$tax_level,
    opt$prevalence_cutoff,
    opt$library_cutoff,
    opt$alpha,
    opt$p_adjust,
    opt$s0_perc,
    opt$cores,
    opt$pseudo_sens,
    opt$structural_zero,
    opt$negative_lower_bound,
    opt$global,
    opt$pairwise,
    opt$dunnet,
    opt$trend,
    opt$top_taxa,
    opt$seed
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  parameter_table,
  file.path(logs_dir, "analysis-parameters.tsv")
)

sink(file.path(logs_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()


# ==============================================================================
# 16. FINISH
# ==============================================================================

section("ANCOM-BC2 analysis completed")

cat("Input directory:  ", input_dir, "\n", sep = "")
cat("Output directory: ", normalizePath(output_dir), "\n", sep = "")
cat("ANCOMBC version:  ", as.character(packageVersion("ANCOMBC")), "\n", sep = "")
cat("phyloseq version: ", as.character(packageVersion("phyloseq")), "\n", sep = "")
cat("ggplot2 version:  ", as.character(packageVersion("ggplot2")), "\n", sep = "")

cat("\nGenerated directories:\n")
cat("  ", tables_dir, "\n", sep = "")
cat("  ", figures_dir, "\n", sep = "")
cat("  ", objects_dir, "\n", sep = "")
cat("  ", logs_dir, "\n", sep = "")
