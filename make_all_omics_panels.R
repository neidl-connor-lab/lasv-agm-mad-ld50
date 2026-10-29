#!/usr/bin/env Rscript
# =============================================================================
# LASV NanoString / LEGENDplex figures -- Figures 4, 5, S1 and S2
#
#   Rscript make_all_panels.R [data-directory]
#
# Self-contained: shared setup and helpers first, then one section per figure.
# Needs data.xlsx and the ipa/ directory (Figure 4F) in the data directory.
#
# Everything is written into analysis/:
# Data comes from a single workbook with sheets:
#   animals                   one row per animal (ID, outcome)
#   nanostring.samplesheet    transcriptomics metadata
#   nanostring.thresholded    transcriptomics counts
#   legendplex.samplesheet    proteomics metadata
#   legendplex.concentrations proteomics concentrations
#   virology.blood            blood viral load by DPI (Figure 4B)
# IPA results for Figure 4F are read from the ipa/ directory.
# =============================================================================


## ---- configuration ----------------------------------------------------------
DATA.XLSX    <- "data.xlsx"
ANALYSIS.DIR <- "analysis"
FIG.DIR      <- "analysis"   # figures and panels/ sit beside the DE tables
IPA.DIR      <- "ipa"       # IPA exports used by Figure 4F
BASE.SIZE    <- 9       # font size for journal-scale panels
TAG.SIZE     <- 13      # panel letter size
FONT.FAMILY  <- "Arial"

# Data directory resolution order:
#   1. first command-line argument   Rscript make_all_panels.R /path/to/data
#   2. the LASV_DATA environment variable
#   3. the directory this script lives in
script.dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) return(dirname(normalizePath(f[1])))
  getwd()
}
if (!exists("DATA.DIR") || is.null(DATA.DIR)) {
  cli <- commandArgs(trailingOnly = TRUE)
  DATA.DIR <- if (length(cli) && nzchar(cli[1])) {
    cli[1]
  } else if (nzchar(Sys.getenv("LASV_DATA"))) {
    Sys.getenv("LASV_DATA")
  } else {
    script.dir()
  }
}


## ---- packages ---------------------------------------------------------------
.required <- c("tidyverse", "readxl", "edgeR", "limma", "ggpubr", "ggrepel",
               "patchwork")
.missing <- .required[!vapply(.required, requireNamespace, logical(1),
                              quietly = TRUE)]
if (length(.missing)) {
  .bioc <- intersect(.missing, c("edgeR", "limma"))
  .cran <- setdiff(.missing, .bioc)
  .msg <- c("Missing packages: ", paste(.missing, collapse = ", "), "\n",
            "R version : ", R.version.string, "\n",
            ".libPaths():\n", paste0("  ", .libPaths(), collapse = "\n"), "\n\n")
  if (length(.cran)) .msg <- c(.msg, '  install.packages(c(',
                               paste0('"', .cran, '"', collapse = ", "), '))\n')
  if (length(.bioc)) .msg <- c(.msg, '  BiocManager::install(c(',
                               paste0('"', .bioc, '"', collapse = ", "), '))\n')
  stop(.msg, call. = FALSE)
}

# patchwork < 1.2.0 cannot lay out ggplot2 >= 3.5 plots (the guide and panel
# gtables changed shape); the failure surfaces deep inside ggsave().
.pw <- utils::packageVersion("patchwork")
.gg <- utils::packageVersion("ggplot2")
if (.pw < "1.2.0" && .gg >= "3.5.0")
  stop("patchwork ", .pw, " cannot lay out ggplot2 ", .gg, " plots.\n",
       "  install.packages('patchwork')   # want >= 1.2.0", call. = FALSE)

options(lifecycle_verbosity = "quiet")   # silence ggpubr's internal deprecations
suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(patchwork)
})

if (!dir.exists(DATA.DIR))
  stop("Data directory does not exist: ", DATA.DIR, "\n",
       "  Pass one explicitly:  Rscript <script>.R /path/to/data", call. = FALSE)
setwd(DATA.DIR)
dir.create(ANALYSIS.DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(FIG.DIR, "panels"), showWarnings = FALSE, recursive = TRUE)
message("Data directory: ", getwd())
message("Using ggplot2 ", .gg, ", patchwork ", .pw)


## ---- fonts and theme --------------------------------------------------------
.font.candidates <- c(FONT.FAMILY, "Liberation Sans", "Arimo",
                      "Helvetica", "Nimbus Sans", "FreeSans")
if (requireNamespace("systemfonts", quietly = TRUE)) {
  .fdir <- path.expand("~/.fonts")
  if (dir.exists(.fdir) && "add_fonts" %in% getNamespaceExports("systemfonts"))
    try(systemfonts::add_fonts(.fdir), silent = TRUE)
  .have <- tolower(unique(systemfonts::system_fonts()$family))
  .hit  <- .font.candidates[tolower(.font.candidates) %in% .have]
  if (length(.hit)) {
    if (!identical(.hit[1], .font.candidates[1]))
      message("Note: '", .font.candidates[1], "' not installed - using '",
              .hit[1], "' (metric-compatible with Arial).")
    FONT.FAMILY <- .hit[1]
  } else {
    FONT.FAMILY <- "sans"
    message("Note: no Arial-metric font found; using generic 'sans'.")
  }
}
message("Font family: ", FONT.FAMILY)

theme_set(ggpubr::theme_pubr(base_size = BASE.SIZE, base_family = FONT.FAMILY) +
            theme(legend.position  = "right",
                  strip.background = element_blank(),
                  text     = element_text(family = FONT.FAMILY),
                  plot.tag = element_text(size = TAG.SIZE, face = "bold",
                                          family = FONT.FAMILY)))
update_geom_defaults("text",  list(family = FONT.FAMILY))
update_geom_defaults("label", list(family = FONT.FAMILY))


## ---- constants --------------------------------------------------------------
cols.stage <- c(Baseline  = "#bcbcbc",
                Early     = "#009444",
                Middle    = "#e69138",
                Late      = "#e06666",
                Recovered = "blue")
cols.reg       <- c(Up = "#e41a1c", Down = "#377eb8")
shapes.group   <- c("Fast" = 23, "Slow" = 22, "Non-responder" = 21)
shapes.outcome <- c("Survived" = 21, "Succumbed" = 24)
cols.outcome   <- c("Survived" = "#2c7fb8", "Succumbed" = "#c9314e")

# Day palette, shared by every PCA. Baseline is white
cols.day <- c("#FFC145",   # 3-4   gold
              "#F4791F",   # 7     orange
              "#E5383B",   # 10    red
              "#B5179E",   # 14-17 magenta
              "#7209B7",   # 21    violet
              "#4361EE",   # 28    bright blue
              "#10162F")   # 35    near-black navy
# Bin edges (left-closed): 0 | 3,4 | 7 | 10 | 14,15,17 | 21 | 28 | 35
day.edges    <- c(1.5, 5.5, 8.5, 12, 18.5, 24.5, 31.5)
breaks.dates <- c(0, 7, 14, 21, 28, 35)

# Responder groups (IDs as in the 'animals' sheet)
fast.animals <- c("Sym-1", "Sym-3", "Su-4")
slow.animals <- c("Sym-2", "Su-1", "Su-2", "Su-3")

# NanoString samples that failed QC, as nanostring.samplesheet IDs. They are
# present in the thresholded sheet and excluded at load time.
#   sample81 Sym-1 d3    sample82 Sym-1 d28   sample83 Asym-3 d3
#   sample84 Asym-3 d35  sample85 Sym-2 d3    sample86 Sym-2 d28
#   sample87 Su-4 d3     sample88 Su-1 d3     sample89 Asym-5 d3
#   sample92 Sym-3 d3    sample93 Sym-3 d35   sample94 Asym-2 d3
#   sample95 Su-3 d3     sample96 Asym-4 d3
qc.fail <- c(sprintf("sample%02d", 81:89), sprintf("sample%02d", 92:96))

# Duplicated LEGENDplex analytes to drop (the poorer-performing panel)
dup.analytes <- c("IFNg_cc", "IL1b_cc", "IL6_cc", "IL6_th", "IP10_inf",
                  "MCP1_cc", "TNFa_cc", "IL8_cc")


## ---- workbook access --------------------------------------------------------
.norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))

xlsx.path <- function() {
  p <- DATA.XLSX
  if (!file.exists(p))
    stop("Workbook not found: ", normalizePath(p, mustWork = FALSE), "\n",
         "  getwd(): ", getwd(), "\n",
         "  files here: ",
         paste(head(list.files(pattern = "\\.xlsx?$"), 10), collapse = " | "),
         "\n  Set DATA.XLSX, or pass the data directory as the first argument.",
         call. = FALSE)
  p
}

sheet.names <- function() readxl::excel_sheets(xlsx.path())

# Loose sheet lookup: "nanostring.samplesheet" also matches
# "NanoString_SampleSheet", "nanostring samplesheet", etc.
find.sheet <- function(candidates) {
  sh <- sheet.names()
  for (cand in candidates) {
    hit <- sh[.norm(sh) == .norm(cand)]
    if (length(hit)) return(hit[1])
  }
  NULL
}

# Read a table from the workbook
read.any <- function(sheets, files = character(), required = TRUE,
                     what = sheets[1]) {
  sh <- find.sheet(sheets)
  if (!is.null(sh)) {
    message("  reading sheet '", sh, "'")
    return(as.data.frame(readxl::read_excel(xlsx.path(), sheet = sh,
                                            .name_repair = "unique_quiet")))
  }
  for (f in files) {
    if (file.exists(f)) {
      message("  reading file '", f, "'")
      return(if (grepl("\\.tsv$", f, ignore.case = TRUE))
               read.delim(f) else read.csv(f))
    }
  }
  if (!required) return(NULL)
  stop("Could not find '", what, "'.\n",
       "  Sheets in ", DATA.XLSX, ": ",
       paste(sheet.names(), collapse = " | "), "\n",
       if (length(files))
         paste0("  Also looked for file(s): ", paste(files, collapse = ", "), "\n"),
       call. = FALSE)
}

# Resolve a column by any of several candidate names
pick.col <- function(df, candidates, what, required = TRUE) {
  nm <- names(df)
  for (cand in candidates) {
    hit <- nm[.norm(nm) == .norm(cand)]
    if (length(hit)) return(hit[1])
  }
  if (!required) return(NA_character_)
  stop("Could not find a '", what, "' column.\n",
       "  Tried: ", paste(candidates, collapse = ", "), "\n",
       "  Present: ", paste(nm, collapse = " | "), call. = FALSE)
}


## ---- generic helpers --------------------------------------------------------
# Concentrations can be zero (below limit of detection)
# Pseudocount applied only when needed.
log2.safe <- function(mat) {
  mat <- as.matrix(mat)
  nonpos <- mat <= 0
  if (any(nonpos, na.rm = TRUE)) {
    pseudo <- min(mat[mat > 0], na.rm = TRUE) / 2
    message("Note: ", sum(nonpos, na.rm = TRUE), " non-positive value(s); ",
            "adding pseudocount ", signif(pseudo, 3), " before log2.")
    mat <- mat + pseudo
  }
  log2(mat)
}

format.counts <- function(counts.mat, n = 3, threshold = 20) {
  counts.mat[rowSums(counts.mat > threshold, na.rm = TRUE) >= n, , drop = FALSE]
}

cpm.log <- function(counts.mat) {
  counts.mat %>% DGEList() %>% calcNormFactors() %>% cpm(log = TRUE)
}

run.pca <- function(mat, meta.mat) {
  pca <- prcomp(t(mat))
  pcs <- summary(pca)$importance["Proportion of Variance", 1:2]
  pcs <- paste0(names(pcs), " (", round(100 * pcs), "%)")
  contribs <- pca$rotation %>%
    as.data.frame() %>% rownames_to_column("Gene") %>% select(Gene, PC1, PC2)
  scores <- pca$x %>%
    as.data.frame() %>% rownames_to_column("Sample.ID") %>%
    select(Sample.ID, PC1, PC2) %>% left_join(meta.mat, by = "Sample.ID")
  list(pca = scores, contribs = contribs, pcs = pcs)
}

# Animal blocking
use.block <- function(block, n) {
  if (is.null(block)) return(FALSE)
  if (length(block) != n)
    stop("block has ", length(block), " entries for ", n, " samples.",
         call. = FALSE)
  any(duplicated(block))
}
report.correlation <- function(rho)
  message("  within-animal correlation (consensus): ", round(rho, 3))

# counts -> voom -> limma. Blocked fits use the "voom twice" pipeline
run.de.counts <- function(counts.mat, model.mat, contrasts.obj, block = NULL) {
  dge <- counts.mat %>% DGEList() %>% calcNormFactors()
  if (!use.block(block, ncol(dge))) {
    fit <- dge %>% voom(model.mat) %>% lmFit()
  } else {
    v   <- voom(dge, model.mat)
    rho <- duplicateCorrelation(v, model.mat, block = block)$consensus.correlation
    v   <- voom(dge, model.mat, block = block, correlation = rho)
    rho <- duplicateCorrelation(v, model.mat, block = block)$consensus.correlation
    report.correlation(rho)
    fit <- lmFit(v, model.mat, block = block, correlation = rho)
  }
  fit %>% contrasts.fit(contrasts.obj) %>% eBayes(robust = TRUE)
}

# already-log concentrations -> limma
run.de.conc <- function(logconc.mat, model.mat, contrasts.obj, block = NULL) {
  if (!use.block(block, ncol(logconc.mat))) {
    fit <- lmFit(logconc.mat, model.mat)
  } else {
    rho <- duplicateCorrelation(logconc.mat, model.mat,
                                block = block)$consensus.correlation
    report.correlation(rho)
    fit <- lmFit(logconc.mat, model.mat, block = block, correlation = rho)
  }
  fit %>% contrasts.fit(contrasts.obj) %>% eBayes()
}

format.results <- function(coef, ebayes.output, psig = 0.05, fsig = 1) {
  degs <- limma::topTable(ebayes.output, coef, number = Inf) %>%
    rownames_to_column("Gene") %>%
    dplyr::rename(lfc = logFC, padj = adj.P.Val) %>%
    mutate(Contrast    = coef,
           # facet key: "Early - Baseline" -> "Early"
           Stage       = str_extract(coef, "^[A-Za-z]+"),
           Significant = (padj < psig) & (abs(lfc) > fsig)) %>%
    select(Gene, lfc, padj, Contrast, Stage, Significant)
  degs$Regulation <- "None"
  degs$Regulation[degs$Significant] <- "Up"
  degs$Regulation[degs$Significant & degs$lfc < 0] <- "Down"
  degs$Regulation <- factor(degs$Regulation, levels = c("Up", "None", "Down"))
  degs
}

# Bin DPI into discrete groups, labelling each by the days actually observed
make.day.bins <- function(day, edges = day.edges, baseline.label = "Baseline") {
  b <- cut(day, breaks = c(-Inf, edges, Inf), right = FALSE)
  labs <- tapply(day, b, function(v) {
    v <- sort(unique(v[!is.na(v)]))
    if (!length(v)) NA_character_
    else if (length(v) == 1) as.character(v)
    else paste0(min(v), "-", max(v))
  })
  labs <- as.character(labs)
  out <- factor(labs[as.integer(b)], levels = labs[!is.na(labs)])
  if (!is.null(baseline.label))
    levels(out)[levels(out) == "0"] <- baseline.label
  out
}

day.palette <- function(bins) {
  n  <- nlevels(bins)
  hb <- "Baseline" %in% levels(bins)
  np <- n - as.integer(hb)
  post <- if (np == length(cols.day)) {
    cols.day
  } else {
    grDevices::colorRampPalette(cols.day)(max(np, 1))
  }
  pal <- c(if (hb) "#ffffff", post)[seq_len(n)]
  setNames(pal, levels(bins))
}

# Standard PCA panel
pca.panel <- function(scores, pcs, shape.var, shape.values, day.var = "DPI") {
  scores$Day.bin <- make.day.bins(scores[[day.var]])
  pal <- day.palette(scores$Day.bin)
  message("  day bins: ", paste(levels(scores$Day.bin), collapse = ", "))
  ggplot(scores, aes(PC1, PC2, fill = .data$Day.bin,
                     shape = .data[[shape.var]])) +
    geom_hline(yintercept = 0, linetype = 3, col = "lightgrey") +
    geom_vline(xintercept = 0, linetype = 3, col = "lightgrey") +
    geom_point(size = 2.6, colour = "grey20", stroke = 0.4) +
    scale_shape_manual(values = shape.values) +
    scale_fill_manual(values = pal, drop = FALSE) +
    labs(x = pcs[1], y = pcs[2], shape = NULL, fill = "DPI") +
    guides(fill  = guide_legend(order = 1, ncol = 1,
                                override.aes = list(shape = 21, size = 2.6,
                                                    colour = "grey20",
                                                    stroke = 0.4)),
           shape = guide_legend(order = 2, override.aes =
                                  list(fill = "grey70", size = 2.6))) +
    theme(legend.position   = "right",
          legend.title      = element_text(size = BASE.SIZE - 1),
          legend.text       = element_text(size = BASE.SIZE - 1),
          legend.key.height = unit(3.4, "mm"),
          legend.margin     = margin(0, 0, 0, 0),
          legend.spacing.y  = unit(1.5, "mm"))
}

volcano.plot <- function(res, facet = TRUE, n.up = 8, n.down = 0,
                         xlim = c(-3, 8), xbreaks = c(-3, 0, 3, 6),
                         ylim = c(NA, 10), labels.stage = NULL) {
  pick <- function(reg, n) {
    if (n <= 0) return(res[0, ])
    df <- filter(res, Regulation == reg)
    if (facet) df <- group_by(df, Stage)
    # slice_min for Down: slice_max(lfc) returns the LEAST downregulated
    df <- if (reg == "Up") slice_max(df, order_by = lfc, n = n, with_ties = FALSE)
          else             slice_min(df, order_by = lfc, n = n, with_ties = FALSE)
    ungroup(df)
  }
  topgenes <- bind_rows(pick("Up", n.up), pick("Down", n.down))

  p <- res %>%
    arrange(Significant) %>%
    ggplot(aes(lfc, -log10(padj))) +
    geom_point(aes(size = Regulation, col = Regulation), alpha = 0.6) +
    scale_color_manual(values = cols.reg, na.value = "lightgrey") +
    scale_size_manual(values = c(Up = 1, Down = 1), na.value = 0.5) +
    ggrepel::geom_text_repel(data = topgenes, aes(label = Gene),
                             size = 2, fontface = "bold", force = 50,
                             nudge_y = 1, min.segment.length = 2,
                             max.overlaps = Inf) +
    scale_x_continuous(limits = xlim, breaks = xbreaks) +
    scale_y_continuous(limits = ylim) +
    labs(x = "Fold change (log2)", y = "padj (-log10)") +
    theme(legend.position = "none")
  if (facet) {
    p <- p + if (is.null(labels.stage)) facet_wrap(~ Stage)
             else facet_wrap(~ Stage, labeller = labeller(Stage = labels.stage))
  }
  # fixed limits DROP out-of-range points rather than squashing them
  off <- sum(res$lfc < xlim[1] | res$lfc > xlim[2] |
               (-log10(res$padj) > ylim[2]), na.rm = TRUE)
  if (off)
    warning(off, " point(s) fall outside the fixed axis limits and are not ",
            "drawn. Widen xlim/ylim if that matters.", call. = FALSE)
  p
}

# Shared panel styling
tagged <- function(p, lab = NULL, title.hjust = 0.5) {
  p <- p + theme(plot.tag = element_text(size = TAG.SIZE, face = "bold",
                                         hjust = 0, vjust = 1),
                 plot.tag.position = c(0, 1),
                 plot.title = element_text(size = BASE.SIZE, face = "bold",
                                           hjust = title.hjust),
                 plot.margin = margin(6, 6, 4, 6))
  if (!is.null(lab)) p <- p + labs(tag = lab)
  p
}

pdf.device <- function() {
  if (!isTRUE(capabilities("cairo"))) {
    warning("cairo not available - fonts may not embed correctly in the PDF.",
            call. = FALSE)
    grDevices::pdf
  } else grDevices::cairo_pdf
}

save.figure <- function(plot, name, width, height, panels = NULL) {
  dev <- pdf.device()
  ggsave(file.path(FIG.DIR, paste0(name, ".pdf")), plot,
         width = width, height = height, units = "mm", device = dev)
  ggsave(file.path(FIG.DIR, paste0(name, ".png")), plot,
         width = width, height = height, units = "mm", dpi = 400, bg = "white")
  if (!is.null(panels)) {
    for (nm in names(panels)) {
      ggsave(file.path(FIG.DIR, "panels", paste0(name, nm, ".pdf")),
             panels[[nm]]$plot, width = panels[[nm]]$size[1],
             height = panels[[nm]]$size[2], units = "mm", device = dev)
    }
  }
  message("Written: ", normalizePath(file.path(FIG.DIR, paste0(name, ".pdf")),
                                     mustWork = FALSE))
}

# Collapse the deferred warning list into unique messages 
finish.up <- function(name) {
  w <- warnings()
  if (length(w)) {
    tab <- sort(table(names(w)), decreasing = TRUE)
    message("\n--- warning summary (", length(w), " total) ---")
    for (i in seq_along(tab)) message("  ", tab[[i]], "x  ", names(tab)[i])
    message("--- end warning summary ---\n")
  }
  writeLines(capture.output(sessionInfo()),
             file.path(ANALYSIS.DIR, paste0("sessionInfo-", name, ".txt")))
}


## ---- shared metadata formatting ---------------------------------------------
# Stage from DPI, used when the sample sheet has no Stage column.

stage.from.dpi <- function(dpi, recovered.from = 35) {
  out <- rep(NA_character_, length(dpi))
  out[dpi == 0]                   <- "Baseline"
  out[dpi %in% c(3, 4)]           <- "Early"
  out[dpi == 7]                   <- "Middle"
  out[dpi > 7 & dpi < recovered.from]  <- "Late"
  out[dpi >= recovered.from]      <- "Recovered"
  out
}

finish.meta <- function(meta) {
  meta <- meta %>%
    mutate(Outcome = str_to_title(str_trim(as.character(Outcome))),
           Group   = case_when(Animal.ID %in% fast.animals ~ "Fast",
                               Animal.ID %in% slow.animals ~ "Slow",
                               TRUE                        ~ "Non-responder"),
           Stage   = factor(Stage,   levels = names(cols.stage)),
           Outcome = factor(Outcome, levels = names(shapes.outcome)),
           Group   = factor(Group,   levels = names(shapes.group)))
  if (any(is.na(meta$Outcome)))
    stop("Unrecognised Outcome value(s); expected Survived/Succumbed.",
         call. = FALSE)
  meta <- as.data.frame(meta)
  rownames(meta) <- meta$Sample.ID
  meta
}

as.feature.matrix <- function(df, sample.ids, what = "counts") {
  id.col <- names(df)[1]
  feats  <- as.character(df[[id.col]])
  num    <- df[, setdiff(names(df), id.col), drop = FALSE]
  num[]  <- suppressWarnings(lapply(num, function(v) as.numeric(as.character(v))))
  mat    <- as.matrix(num)
  rownames(mat) <- feats

  # orientation A: features in rows, samples in columns
  nA <- length(intersect(colnames(mat), sample.ids))
  # orientation B: samples in rows
  nB <- length(intersect(feats, sample.ids))

  if (nB > nA) {
    message("  detected samples-in-rows layout; transposing")
    mat <- t(mat)
  }
  if (!length(intersect(colnames(mat), sample.ids)))
    stop("No ", what, " columns matched the sample sheet.\n",
         "  first ", what, " IDs : ",
         paste(head(colnames(mat), 5), collapse = ", "), "\n",
         "  first sheet IDs: ", paste(head(sample.ids, 5), collapse = ", "),
         call. = FALSE)
  nas <- sum(is.na(mat))
  if (nas)
    message("  note: ", nas, " value(s) are NA after numeric coercion (",
            round(100 * nas / length(mat), 2), "%)")
  mat
}


# Per-animal table ('animals' sheet): one row per NHP with its Outcome.
load.animals <- function() {
  an <- read.any(c("animals", "NHP"), what = "animals")
  c.id  <- pick.col(an, c("NHP", "NHP.ID", "animal id", "animal"), "animal ID")
  c.out <- pick.col(an, c("Outcome", "survival"), "outcome")
  data.frame(Animal.ID = trimws(as.character(an[[c.id]])),
             Outcome   = trimws(as.character(an[[c.out]])),
             stringsAsFactors = FALSE)
}

# Outcome per sample: from the sample sheet if it has the column, otherwise
# joined from the 'animals' sheet by animal ID.
outcome.lookup <- function(ss, animal.ids) {
  c.out <- pick.col(ss, c("Outcome", "survival"), "outcome", required = FALSE)
  if (!is.na(c.out)) return(as.character(ss[[c.out]]))
  message("  no Outcome column; joining from the 'animals' sheet")
  an  <- load.animals()
  out <- an$Outcome[match(animal.ids, an$Animal.ID)]
  if (anyNA(out))
    stop("Animal(s) missing from the 'animals' sheet: ",
         paste(unique(animal.ids[is.na(out)]), collapse = ", "), call. = FALSE)
  out
}


## ---- dataset loaders --------------------------------------------------------
# Both return list(meta = data.frame, mat = matrix) aligned on Sample.ID, with
# a standard schema: Sample.ID, Animal.ID, DPI, Stage, Outcome, Group.

load.nanostring <- function() {
  message("Loading NanoString data")
  ss <- read.any(c("nanostring.samplesheet", "nanostring samplesheet"),
                 what = "nanostring.samplesheet")
  ct <- read.any(c("nanostring.thresholded", "nanostring thresholded"),
                 what = "nanostring.thresholded")

  c.sample <- pick.col(ss, c("Sample.ID", "sample", "sampleid", "sample name",
                             "ID"), "sample ID")
  c.animal <- pick.col(ss, c("NHP.ID", "animal id", "animal", "nhp", "subject"),
                       "animal ID")
  c.dpi    <- pick.col(ss, c("Study.day", "dpi", "day", "timepoint"), "DPI")
  c.stage  <- pick.col(ss, c("Stage", "disease stage"), "stage", required = FALSE)

  meta <- data.frame(
    Sample.ID = as.character(ss[[c.sample]]),
    Animal.ID = trimws(as.character(ss[[c.animal]])),
    DPI       = suppressWarnings(as.numeric(as.character(ss[[c.dpi]]))),
    stringsAsFactors = FALSE)
  meta$Outcome <- outcome.lookup(ss, meta$Animal.ID)

  meta$Stage <- if (!is.na(c.stage)) as.character(ss[[c.stage]]) else {
    message("  no Stage column; deriving from DPI (Late runs to 28)")
    stage.from.dpi(meta$DPI, recovered.from = 35)
  }

  dropped <- sum(.norm(meta$Sample.ID) %in% .norm(qc.fail))
  if (dropped) {
    message("  excluding ", dropped, " QC-failed sample(s)")
    meta <- meta[!(.norm(meta$Sample.ID) %in% .norm(qc.fail)), , drop = FALSE]
  }
  meta <- finish.meta(meta)

  # drop NanoString positive/negative spike-in probes if still present
  probe.col <- names(ct)[1]
  spikes <- grepl("POS_|NEG_", as.character(ct[[probe.col]]))
  if (any(spikes)) {
    message("  dropping ", sum(spikes), " spike-in control probe(s)")
    ct <- ct[!spikes, , drop = FALSE]
  }
  mat <- as.feature.matrix(ct, meta$Sample.ID, "count")

  keep <- intersect(colnames(mat), rownames(meta))
  meta <- meta[keep, , drop = FALSE]
  mat  <- mat[, keep, drop = FALSE]
  message("Loaded ", nrow(mat), " probes x ", ncol(mat), " samples.")
  list(meta = meta, mat = mat)
}

load.legendplex <- function() {
  message("Loading LEGENDplex data")
  ss <- read.any(c("legendplex.samplesheet", "legendplex samplesheet"),
                 what = "legendplex.samplesheet")
  cc <- read.any(c("legendplex.concentrations", "legendplex concentrations"),
                 what = "legendplex.concentrations")

  c.sample <- pick.col(ss, c("sample", "Sample.ID", "sampleid", "sample name",
                             "ID"), "sample ID")
  c.animal <- pick.col(ss, c("animal id", "NHP.ID", "animal", "nhp", "subject"),
                       "animal ID")
  c.dpi    <- pick.col(ss, c("dpi", "Study.day", "day", "timepoint"), "DPI")
  c.stage  <- pick.col(ss, c("Stage", "stage"), "stage", required = FALSE)

  meta <- data.frame(
    Sample.ID = as.character(ss[[c.sample]]),
    Animal.ID = trimws(as.character(ss[[c.animal]])),
    DPI       = suppressWarnings(as.numeric(as.character(ss[[c.dpi]]))),
    stringsAsFactors = FALSE)
  meta$Outcome <- outcome.lookup(ss, meta$Animal.ID)
  meta$Stage <- if (!is.na(c.stage)) as.character(ss[[c.stage]]) else {
    message("  no Stage column; deriving from DPI (Recovered starts at 28)")
    stage.from.dpi(meta$DPI, recovered.from = 28)
  }
  meta <- finish.meta(meta)

  # average technical replicates if the sheet still carries them
  rep.col <- pick.col(cc, c("replicate", "rep", "tech rep"), "replicate",
                      required = FALSE)
  samp.col <- pick.col(cc, c("sample", "Sample.ID", "sampleid", "ID"), "sample",
                       required = FALSE)
  if (!is.na(rep.col) && !is.na(samp.col)) {
    message("  averaging technical replicates")
    cc <- cc %>%
      group_by(.data[[samp.col]]) %>%
      summarise(across(where(is.numeric), \(v) mean(v, na.rm = TRUE)),
                .groups = "drop") %>%
      select(-any_of(c(rep.col, "X"))) %>%
      as.data.frame()
  }

  mat <- as.feature.matrix(cc, meta$Sample.ID, "concentration")

  drop.now <- intersect(rownames(mat), dup.analytes)
  if (length(drop.now)) {
    message("  dropping ", length(drop.now), " duplicated analyte(s)")
    mat <- mat[!(rownames(mat) %in% drop.now), , drop = FALSE]
  }
  # strip a trailing panel suffix only if one is actually present
  if (any(grepl("_", rownames(mat)))) {
    stripped <- sub("_[^_]+$", "", rownames(mat))
    if (!anyDuplicated(stripped)) {
      message("  stripping panel suffixes from analyte names")
      rownames(mat) <- stripped
    } else {
      warning("Panel suffixes left in place: stripping them would create ",
              "duplicate analyte names.", call. = FALSE)
    }
  }

  keep <- intersect(colnames(mat), rownames(meta))
  meta <- meta[keep, , drop = FALSE]
  mat  <- mat[, keep, drop = FALSE]
  message("Loaded ", nrow(mat), " analytes x ", ncol(mat), " samples.")
  list(meta = meta, mat = mat)
}

# Report the animals present, so a mis-joined sheet is obvious immediately
report.cohort <- function(meta) {
  a <- meta[!duplicated(meta$Animal.ID), ]
  message("Animals by group: ",
          paste(names(table(a$Group)), table(a$Group), sep = "=",
                collapse = ", "))
}

# Everything defined above is shared. Each figure section below starts from
# this state: objects a section creates (meta, cmat, labels.stage, FIG.WIDTH,
# ...) are removed before the next one runs, so nothing leaks between figures.
.shared <- c(ls(all.names = TRUE), ".shared")
clear.figure <- function()
  rm(list = setdiff(ls(globalenv(), all.names = TRUE), .shared), envir = globalenv())
.shared <- c(.shared, "clear.figure")


# =============================================================================
# NanoString - Figure 4
#
# Reads data.xlsx (sheets: animals, nanostring.samplesheet,
# nanostring.thresholded, virology.blood) and the IPA results in ipa/.
#
# Panels:
#   A  PCA of the seven viremic animals, colored by DPI
#   B  Viremia onset and response pace by outcome
#   C  Fast responders: each stage vs Baseline
#   D  Slow responders: each stage vs Baseline
#   E  Overlap of the two initial responses
#   F  IPA radar
# =============================================================================


FIG.WIDTH  <- 180   # mm
FIG.HEIGHT <- 230   # mm


labels.stage <- c(Baseline  = "Baseline",
                  Early     = "3-4 DPI",
                  Middle    = "7 DPI",
                  Late      = "10-28 DPI",
                  Recovered = "Recovered")

fig4 <- list()


## ---- data --------------------------------------------------------------------
ns   <- load.nanostring()
meta <- ns$meta
cmat <- ns$mat
report.cohort(meta)


## ---- panel A: PCA of viremic animals ----------------------------------------
keep.ids <- c(fast.animals, slow.animals)
meta.sub <- meta[meta$Animal.ID %in% keep.ids, , drop = FALSE]
if (!nrow(meta.sub))
  stop("No samples from the viremic animals; check Animal.ID formatting.",
       call. = FALSE)

pca <- cmat[, meta.sub$Sample.ID, drop = FALSE] %>%
  format.counts() %>% cpm.log() %>% run.pca(meta.sub)

fig4$A <- pca$pca %>%
  arrange(Group == "Fast") %>%            # fast responders drawn on top
  pca.panel(pca$pcs, shape.var = "Outcome", shape.values = shapes.outcome)


## ---- panel B: viremia kinetics (two plots, one panel) -----------------------
# Viremia onset is the first DPI with detectable genomes in blood; pace is the
# responder group and outcome comes from the 'animals' sheet.
vb <- read.any("virology.blood", what = "virology.blood")
c.nhp <- pick.col(vb, c("NHP", "NHP.ID", "animal"), "animal ID")
c.gq  <- pick.col(vb, c("GEq/mL", "GEq.mL", "GEq"), "GEq/mL")
vb <- vb[!is.na(vb[[c.gq]]) & vb[[c.gq]] > 0, , drop = FALSE]
vk <- load.animals()
vk$GEq <- vapply(vk$Animal.ID, function(a) {
  dpi <- vb$DPI[vb[[c.nhp]] == a]
  if (length(dpi)) min(dpi) else NA_real_
}, numeric(1))
vk$Pace <- case_when(vk$Animal.ID %in% fast.animals ~ "Fast",
                     vk$Animal.ID %in% slow.animals ~ "Slow",
                     TRUE                           ~ "Non-responder")

d <- vk[vk$Pace != "Non-responder", , drop = FALSE]
d$Outcome <- factor(str_to_title(d$Outcome), levels = c("Survived", "Succumbed"))
d$Pace    <- factor(d$Pace, levels = c("Fast", "Slow"))
if (anyNA(d$GEq) || anyNA(d$Outcome) || anyNA(d$Pace))
  warning("viremia onset could not be determined for some animal(s); they are ",
          "not plotted in panel B.", call. = FALSE)

pal <- c(Survived = "#2C6E8E", Succumbed = "#B4472B") 

# spacing
fan <- function(g, width = 0.16) {
  ave(seq_along(g), g, FUN = function(i)
    if (length(i) > 1) seq(-width, width, length.out = length(i)) else 0)
}
# group counts drive the labels, so they cannot drift from the data
n.lab <- function(f) {
  lv  <- levels(f)
  cnt <- table(factor(f, levels = lv))
  paste0(lv, "\n(n=", as.integer(cnt), ")")
}

base.theme <- theme_classic(base_size = BASE.SIZE - 1, base_family = FONT.FAMILY) +
  theme(axis.line  = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        axis.text  = element_text(size = BASE.SIZE - 2),
        axis.title = element_text(size = BASE.SIZE - 1),
        legend.title = element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = BASE.SIZE - 2))

dE <- d
dE$xnum <- as.numeric(dE$Outcome) + fan(interaction(dE$Outcome, dE$GEq))
medE <- aggregate(GEq ~ Outcome, dE, median)
medE$xnum <- as.numeric(medE$Outcome)

# Annotations are MWU p val for viremia and Fisher's exact for pace
fig4$B <- ggplot(dE, aes(xnum, GEq, fill = Outcome)) +
  geom_segment(data = medE,
               aes(x = xnum - 0.30, xend = xnum + 0.30, y = GEq, yend = GEq,
                   colour = Outcome),
               inherit.aes = FALSE, linewidth = 0.7, show.legend = FALSE) +
  geom_point(shape = 21, size = 2.6, colour = "white", stroke = 0.4) +
  annotate("text", x = 1.5, y = 10.55, label = "italic(p)== 0.41",
           parse = TRUE, size = 2.7, family = FONT.FAMILY) +
  scale_x_continuous(breaks = c(1, 2), labels = n.lab(dE$Outcome),
                     limits = c(0.45, 2.55)) +
  scale_y_continuous(breaks = c(7, 10), limits = c(6.2, 10.8)) +
  scale_fill_manual(values = pal) + scale_colour_manual(values = pal) +
  labs(x = NULL, y = "Viremia onset (DPI)") +
  base.theme +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-6, 0, 0, 0))

dF <- d
dF$xnum <- as.numeric(dF$Outcome) + fan(interaction(dF$Outcome, dF$Pace))
dF$ynum <- ifelse(dF$Pace == "Fast", 1, 0)

fig4$B2 <- ggplot(dF, aes(xnum, ynum, fill = Outcome)) +
  geom_point(shape = 21, size = 2.6, colour = "white", stroke = 0.4) +
  annotate("text", x = 1.5, y = 1.42, label = "italic(p)== 0.49",
           parse = TRUE, size = 2.7, family = FONT.FAMILY) +
  scale_x_continuous(breaks = c(1, 2), labels = n.lab(dF$Outcome),
                     limits = c(0.45, 2.55)) +
  scale_y_continuous(breaks = c(0, 1), labels = c("Slow", "Fast"),
                     limits = c(-0.6, 1.6)) +
  scale_fill_manual(values = pal) +
  labs(x = NULL, y = "Response pace") +
  base.theme + theme(legend.position = "none")


## ---- panels C/D: stage vs baseline, by responder group ----------------------
de.by.group <- function(group.name, ipa.file) {
  m <- meta[meta$Group == group.name & !is.na(meta$Stage), , drop = FALSE]
  m$Stage <- droplevels(m$Stage)
  need <- c("Baseline", "Early", "Middle", "Late")
  if (!all(need %in% levels(m$Stage)))
    stop("Group '", group.name, "' is missing stage(s): ",
         paste(setdiff(need, levels(m$Stage)), collapse = ", "), call. = FALSE)

  model <- model.matrix(~ 0 + Stage, data = m)
  colnames(model) <- levels(m$Stage)
  contrs <- limma::makeContrasts(
    contrasts = c("Early - Baseline", "Middle - Baseline", "Late - Baseline"),
    levels = colnames(model))

  fit <- cmat[, rownames(m), drop = FALSE] %>%
    format.counts() %>% run.de.counts(model, contrs, block = m$Animal.ID)

  res <- colnames(contrs) %>%
    lapply(format.results, ebayes.output = fit) %>%
    do.call(rbind, .) %>%
    mutate(Stage = factor(Stage, levels = names(cols.stage)))

  res %>%
    select(Gene, Stage, lfc, padj) %>%
    reshape2::melt(id.vars = c("Gene", "Stage")) %>%
    mutate(Group = paste(Stage, variable, sep = ".")) %>%
    reshape2::dcast(Gene ~ Group, value.var = "value") %>%
    write.csv(ipa.file, row.names = FALSE)
  res
}

res.fast <- de.by.group("Fast", file.path(ANALYSIS.DIR, "ipa-input-fast.csv"))
res.slow <- de.by.group("Slow", file.path(ANALYSIS.DIR, "ipa-input-slow.csv"))

fig4$C <- volcano.plot(res.fast, labels.stage = labels.stage) +
  ggtitle("Fast responders")
fig4$D <- volcano.plot(res.slow, labels.stage = labels.stage) +
  ggtitle("Slow responders")


## ---- panel E: overlap of the two initial responses --------------------------
sig.genes <- function(res, stage)
  res %>% filter(Stage == stage, Significant) %>% pull(Gene) %>% unique()

groups1 <- sig.genes(res.fast, "Early")    # fast responders, Early
groups2 <- sig.genes(res.slow, "Middle")   # slow responders, Middle


venn2 <- function(a, b, set.names = c("Fast", "Slow"),
                  fills = c("steelblue", "firebrick")) {
  circle <- function(x0, r = 1, n = 200) {
    t <- seq(0, 2 * pi, length.out = n)
    data.frame(x = x0 + r * cos(t), y = r * sin(t))
  }
  polys <- rbind(transform(circle(-0.45), Set = set.names[1]),
                 transform(circle( 0.45), Set = set.names[2]))
  polys$Set <- factor(polys$Set, levels = set.names)
  counts <- data.frame(x = c(-1.0, 0, 1.0), y = 0,
                       lab = as.character(c(length(setdiff(a, b)),
                                            length(intersect(a, b)),
                                            length(setdiff(b, a)))))
  titles <- data.frame(x = c(-0.75, 0.75), y = 1.22, lab = set.names)
  ggplot(polys, aes(x, y)) +
    geom_polygon(aes(fill = Set, group = Set), alpha = 0.5,
                 colour = "grey25", linewidth = 0.4) +
    geom_text(data = counts, aes(label = lab), size = 3.6, fontface = "bold") +
    geom_text(data = titles, aes(label = lab), size = 3.3, fontface = "bold") +
    scale_fill_manual(values = setNames(fills, set.names), guide = "none") +
    coord_equal(clip = "off", xlim = c(-1.52, 1.52), ylim = c(-1.08, 1.36)) +
    theme_void(base_size = BASE.SIZE) + theme(legend.position = "none")
}
fig4$E <- venn2(groups1, groups2)


## ---- panel F: IPA pathway radar ---------------------------------------------
for (pkg in c("fmsb", "ggplotify", "reshape2"))
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Panel F needs the '", pkg, "' package: install.packages('", pkg, "')",
         call. = FALSE)

# ipa-shortlist.csv: pathways to plot (Canonical Pathways, Nickname)
# ipa_comps.tsv:     IPA comparison z-scores (Canonical Pathways, Early, Middle)
ipa.shortlist <- read.any("ipa.shortlist", files = file.path(IPA.DIR, "ipa-shortlist.csv"),
                          what = "ipa.shortlist")
ipa.comps <- read.any("ipa.comps", files = file.path(IPA.DIR, "ipa_comps.tsv"),
                      what = "ipa.comps")

ipa.mat <- ipa.comps %>%
  right_join(ipa.shortlist, by = "Canonical.Pathways") %>%
  mutate(Nickname = factor(Nickname, levels = ipa.shortlist$Nickname)) %>%
  arrange(Nickname) %>%
  select(-Canonical.Pathways) %>%
  replace_na(list(Middle = 0, Early = 0)) %>%
  column_to_rownames("Nickname") %>%
  as.matrix() %>% t()

lims <- rbind(Max = rep(4, ncol(ipa.mat)), Min = rep(0, ncol(ipa.mat)))
colnames(lims) <- colnames(ipa.mat)
ipa.mat <- rbind(lims, ipa.mat)[c("Max", "Min", "Early", "Middle"), , drop = FALSE]
ipa.df  <- as.data.frame(ipa.mat)
ipa.df[] <- lapply(ipa.df, as.numeric)


# labels are wrapped onto two lines
wrap.label <- function(x, width = 13) {
  vapply(x, function(s) {
    if (grepl("\n", s, fixed = TRUE)) return(s)
    paste(strwrap(s, width = width), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}
varnames <- wrap.label(str_replace_all(colnames(ipa.df), "\\\\n", "\n"))

fig4$F <- ggplotify::as.ggplot(~{

  op <- par(pty = "s", mar = c(2.5, 5.5, 2.5, 5.5), xpd = NA,
            family = FONT.FAMILY)
  on.exit(par(op), add = TRUE)
  fmsb::radarchart(ipa.df, seg = 4, pcol = cols.stage[c("Early", "Middle")],
                   plty = 1, plwd = 1.8, cglty = 1, cglcol = "lightgrey",
                   vlabels = rep("", length(varnames)),
                   pty = 16, pdensity = 0, calcex = 0.55, palcex = 0.55)

  theta <- (seq(90, 450, length.out = length(varnames) + 1) * pi / 180)[seq_along(varnames)]
  r.lab <- 1.08
  for (i in seq_along(varnames))
    text(r.lab * cos(theta[i]), r.lab * sin(theta[i]), varnames[i], cex = 0.58,
         adj = c((1 - round(cos(theta[i]), 2)) / 2, (1 - round(sin(theta[i]), 2)) / 2))
  # The two IPA columns are the fast responders' Early stage and the slow
  # responders' Middle stage, so the legend names the group, not the stage.
  radar.labels <- c("Fast", "Slow")
  usr <- par("usr")
  lx  <- usr[2] - 0.10 * diff(usr[1:2])
  ly  <- usr[4] - c(0, 0.13 * diff(usr[3:4]))   # raise 0.13 to spread further
  for (i in 1:2)
    legend(x = lx, y = ly[i], legend = radar.labels[i],
           col = cols.stage[c("Early", "Middle")][i],
           lty = 1, lwd = 2, bty = "n", cex = 0.7, seg.len = 1.4, xpd = NA)
}) + theme(aspect.ratio = 1) + coord_cartesian(clip = "off")


## ---- assemble ----------------------------------------------------------------

layout <- "
AAAAAABBBBCCCC
DDDDDDDDDDDDDD
EEEEEEEEEEEEEE
FFFFFFGGGGGGGG
"

fig4$full <-
  tagged(fig4$A, "A") + tagged(fig4$B, "B") + tagged(fig4$B2) +
  tagged(fig4$C, "C") + tagged(fig4$D, "D") +
  tagged(fig4$E, "E") + tagged(fig4$F, "F") +
  plot_layout(design = layout, heights = c(1.15, 1, 1, 1.45))

save.figure(fig4$full, "Figure4", FIG.WIDTH, FIG.HEIGHT,
            panels = list(A  = list(plot = fig4$A,  size = c(90, 80)),
                          B  = list(plot = fig4$B,  size = c(55, 75)),
                          B2 = list(plot = fig4$B2, size = c(48, 75)),
                          C  = list(plot = fig4$C,  size = c(180, 70)),
                          D  = list(plot = fig4$D,  size = c(180, 70)),
                          E  = list(plot = fig4$E,  size = c(75, 70)),
                          F  = list(plot = fig4$F,  size = c(100, 95))))
finish.up("fig4")

clear.figure()


# =============================================================================
# LEGENDplex - Figure 5
#
# Reads data.xlsx (sheets legendplex.samplesheet, legendplex.concentrations).
#
# Panels:
#   A  PCA of log2 concentrations, coloured by DPI
#   B  Top PCA loadings
#   C  Succumbed vs survived, 7-21 DPI
#   D  CRP trajectory from baseline
#   E  Galectin-9 trajectory from baseline
#   F  Baseline vs recovered
# =============================================================================


FIG.WIDTH  <- 200   # mm (180 is the usual two-column max; see note)
FIG.HEIGHT <- 185   # mm

GENE.D <- "CRP"            # analytes plotted as trajectories
GENE.E <- "Galectin9"

fig5 <- list()


## ---- data --------------------------------------------------------------------
lp   <- load.legendplex()
meta <- lp$meta
cmat <- lp$mat
report.cohort(meta)
lmat <- log2.safe(cmat)


## ---- panels A/B: PCA ---------------------------------------------------------
pca <- run.pca(lmat, meta)

fig5$A <- pca.panel(pca$pca, pca$pcs,
                    shape.var = "Outcome", shape.values = shapes.outcome)

fig5$B <- pca$contribs %>%
  mutate(Magnitude = sqrt(PC1^2 + PC2^2)) %>%
  slice_max(n = 30, order_by = Magnitude, with_ties = FALSE) %>%
  ggplot() +
  geom_vline(xintercept = 0, linetype = 3, col = "lightgrey") +
  geom_hline(yintercept = 0, linetype = 3, col = "lightgrey") +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), col = "darkgrey",
               alpha = 0.9, arrow = arrow(length = unit(0.05, "in"))) +
  ggrepel::geom_text_repel(aes(PC1, PC2, label = Gene), size = 2,
                           fontface = "bold", max.overlaps = Inf) +
  labs(x = pca$pcs[1], y = pca$pcs[2]) +
  theme(legend.position = "none")


## ---- a one-contrast DE helper -----------------------------------------------
de.contrast <- function(m, design.var, contrast.str) {
  m <- as.data.frame(m)
  m[[design.var]] <- droplevels(factor(m[[design.var]]))
  model <- model.matrix(reformulate(paste0("0 + ", design.var)), data = m)
  colnames(model) <- levels(m[[design.var]])
  contrs <- limma::makeContrasts(contrasts = contrast.str,
                                 levels = colnames(model))
  fit <- lmat[, m$Sample.ID, drop = FALSE] %>%
    run.de.conc(model, contrs, block = m$Animal.ID)
  colnames(contrs) %>%
    lapply(format.results, ebayes.output = fit) %>%
    do.call(rbind, .)
}


## ---- panel C: succumbed vs survived, 7-21 DPI --------------------------------
m.mid <- meta %>% filter(DPI >= 7, DPI <= 21)
res.outcome <- de.contrast(m.mid, "Outcome", "Succumbed - Survived")
write.csv(res.outcome, file.path(ANALYSIS.DIR, "legendplex-outcome-7to21dpi.csv"),
          row.names = FALSE)

# Contrast is described 
fig5$C <- volcano.plot(res.outcome, facet = FALSE, n.up = 10,
                       xlim = c(-4, 8), xbreaks = c(-4, 0, 4, 8),
                       ylim = c(NA, 20))


## ---- panel F: baseline vs recovered ------------------------------------------

m.ends <- meta %>% filter(Stage %in% c("Baseline", "Recovered"))
res.stage <- de.contrast(m.ends, "Stage", "Baseline - Recovered")
write.csv(res.stage, file.path(ANALYSIS.DIR, "legendplex-baseline-vs-recovered.csv"),
          row.names = FALSE)

fig5$F <- volcano.plot(res.stage, facet = FALSE, n.up = 10, n.down = 10,
                       xlim = c(-4, 8), xbreaks = c(-4, 0, 4, 8),
                       ylim = c(NA, 100))


## ---- panels D/E: per-analyte trajectories ------------------------------------
long.expr <- as.data.frame(t(cmat)) %>%
  rownames_to_column("Sample.ID") %>%
  pivot_longer(-Sample.ID, names_to = "Gene", values_to = "Expression") %>%
  left_join(meta[, c("Sample.ID", "Animal.ID", "DPI", "Stage", "Outcome", "Group")],
            by = "Sample.ID")

df.lfc <- long.expr %>%
  group_by(Animal.ID, Gene) %>%
  mutate(baseline = {
    b <- Expression[DPI == 0]
    if (length(b) == 0) NA_real_ else b[1]   # guard animals with no day-0 draw
  }) %>%
  ungroup() %>%
  mutate(logFC.from.day0 = log2((Expression + 1) / (baseline + 1)))

no.base <- df.lfc %>% filter(is.na(baseline)) %>% pull(Animal.ID) %>% unique()
if (length(no.base))
  message("Note: no day-0 sample for ", paste(no.base, collapse = ", "),
          "; omitted from panels D and E.")

line.plot <- function(gene, show.legend = FALSE) {
  df.sub <- df.lfc %>% filter(Gene == gene, !is.na(logFC.from.day0))
  if (!nrow(df.sub))
    stop("Analyte '", gene, "' not found. Available: ",
         paste(sort(unique(df.lfc$Gene)), collapse = ", "), call. = FALSE)
  p <- ggplot(df.sub, aes(DPI, logFC.from.day0, group = Animal.ID)) +
    geom_hline(yintercept = 0, linetype = 3, col = "lightgrey") +
    geom_line(aes(colour = Outcome), linewidth = 0.5, alpha = 0.8) +
    geom_point(aes(shape = Outcome, fill = Outcome),
               size = 2, colour = "grey20", stroke = 0.4) +
    scale_shape_manual(values = shapes.outcome) +
    scale_fill_manual(values = cols.outcome) +
    scale_colour_manual(values = cols.outcome) +
    scale_x_continuous(breaks = breaks.dates) +
    labs(title = gene, x = "DPI", y = "FC from baseline (log2)",
         shape = NULL, fill = NULL, colour = NULL)
  if (show.legend)
    p + theme(legend.position = "right",
              legend.text = element_text(size = BASE.SIZE - 1),
              legend.key.height = unit(4, "mm"))
  else p + theme(legend.position = "none")
}

fig5$D <- line.plot(GENE.D, show.legend = TRUE)
fig5$E <- line.plot(GENE.E)


## ---- assemble ----------------------------------------------------------------
# Letters are matched to plots ALPHABETICALLY, in the order added.
layout <- "
AAAABBBB
CCCCDDDD
EEEEFFFF
"

fig5$full <-
  tagged(fig5$A, "A", title.hjust = 0) + tagged(fig5$B, "B", title.hjust = 0) +
  tagged(fig5$C, "C", title.hjust = 0) + tagged(fig5$D, "D", title.hjust = 0) +
  tagged(fig5$E, "E", title.hjust = 0) + tagged(fig5$F, "F", title.hjust = 0) +
  plot_layout(design = layout, heights = c(1.1, 1, 1))

save.figure(fig5$full, "Figure5", FIG.WIDTH, FIG.HEIGHT,
            panels = list(A = list(plot = fig5$A, size = c(105, 78)),
                          B = list(plot = fig5$B, size = c(105, 78)),
                          C = list(plot = fig5$C, size = c(100, 72)),
                          D = list(plot = fig5$D, size = c(100, 72)),
                          E = list(plot = fig5$E, size = c(100, 72)),
                          F = list(plot = fig5$F, size = c(100, 72))))
finish.up("fig5")

clear.figure()


# =============================================================================
# NanoString - Supplemental Figure S1
#
# Panels:
#   A  PCA of ALL animals, non-responders included (Figure 4A shows only the
#      seven viremic animals). Shape encodes responder group.
#   B  Succumbed vs survived at each group's initial response timepoint
#      (Fast: Early / 3-4 DPI; Slow: Middle / 7 DPI)
# =============================================================================


FIG.WIDTH  <- 180   # mm
FIG.HEIGHT <- 95    # mm

figS1 <- list()


## ---- data --------------------------------------------------------------------
ns   <- load.nanostring()
meta <- ns$meta
cmat <- ns$mat
report.cohort(meta)


## ---- panel A: PCA of all animals --------------------------------------------
pca <- cmat %>% format.counts() %>% cpm.log() %>% run.pca(meta)

figS1$A <- pca$pca %>%
  arrange(Group != "Non-responder") %>%   # responders drawn on top
  pca.panel(pca$pcs, shape.var = "Group", shape.values = shapes.group)


## ---- panel B: outcome at each group's initial response -----------------------

m.init <- meta %>%
  filter((Group == "Fast" & Stage == "Early") |
         (Group == "Slow" & Stage == "Middle")) %>%
  as.data.frame()
rownames(m.init) <- m.init$Sample.ID
m.init$Outcome <- droplevels(m.init$Outcome)

message("Panel B: ", nrow(m.init), " samples (",
        sum(m.init$Group == "Fast"), " Fast/Early, ",
        sum(m.init$Group == "Slow"), " Slow/Middle); outcomes: ",
        paste(names(table(m.init$Outcome)), table(m.init$Outcome),
              sep = "=", collapse = ", "))

if (nlevels(m.init$Outcome) < 2) {
  warning("Panel B needs both outcomes at the initial timepoint; only ",
          paste(levels(m.init$Outcome), collapse = ", "), " present.",
          call. = FALSE)
  figS1$B <- ggplot() + theme_void()
} else {
  model <- model.matrix(~ 0 + Outcome, data = m.init)
  colnames(model) <- levels(m.init$Outcome)
  contrs <- limma::makeContrasts(contrasts = "Succumbed - Survived",
                                 levels = colnames(model))

  fit <- cmat[, rownames(m.init), drop = FALSE] %>%
    format.counts() %>% run.de.counts(model, contrs, block = m.init$Animal.ID)

  res.init <- colnames(contrs) %>%
    lapply(format.results, ebayes.output = fit) %>%
    do.call(rbind, .)
  write.csv(res.init, file.path(ANALYSIS.DIR, "nanostring-outcome-initial.csv"),
            row.names = FALSE)

  # positive fold change = higher in animals that succumbed
  figS1$B <- volcano.plot(res.init, facet = FALSE, n.up = 8, n.down = 8) +
    ggtitle("Succumbed vs survived, initial response")
}


## ---- assemble ----------------------------------------------------------------
layout <- "
AAAAAAABBBBBBB
"

figS1$full <- tagged(figS1$A, "A") + tagged(figS1$B, "B") +
  plot_layout(design = layout)

save.figure(figS1$full, "FigureS1", FIG.WIDTH, FIG.HEIGHT,
            panels = list(A = list(plot = figS1$A, size = c(105, 85)),
                          B = list(plot = figS1$B, size = c(90, 85))))
finish.up("figS1")

clear.figure()


# =============================================================================
# LEGENDplex - Supplemental Figure S2
#
# Each disease stage versus Baseline, using ALL animals 
# =============================================================================


FIG.WIDTH  <- 180   # mm
FIG.HEIGHT <- 75    # mm

CONTRASTS    <- c("Early - Baseline", "Middle - Baseline", "Late - Baseline")
N.LABEL.UP   <- 10  # analytes labelled per facet
N.LABEL.DOWN <- 0   # set to e.g. 5 to also label downregulated analytes

# Facet labels
labels.stage <- c(Baseline  = "Baseline",
                  Early     = "3-4 DPI",
                  Middle    = "7 DPI",
                  Late      = "10-21 DPI",
                  Recovered = "28-35 DPI")


## ---- data --------------------------------------------------------------------
lp   <- load.legendplex()
meta <- lp$meta
cmat <- lp$mat
report.cohort(meta)
lmat <- log2.safe(cmat)


## ---- differential expression: every stage vs Baseline ------------------------
m <- meta %>% filter(!is.na(Stage)) %>% as.data.frame()
rownames(m) <- m$Sample.ID
dropped <- nrow(meta) - nrow(m)
if (dropped)
  message("Note: ", dropped, " sample(s) at unstaged DPI excluded from the model.")

m$Stage <- droplevels(m$Stage)
message("Samples per stage: ",
        paste(names(table(m$Stage)), table(m$Stage), sep = "=", collapse = ", "))

model <- model.matrix(~ 0 + Stage, data = m)
colnames(model) <- levels(m$Stage)

# Keep only contrasts whose terms survived droplevels
want <- CONTRASTS[vapply(CONTRASTS, function(cc)
  all(trimws(strsplit(cc, "-", fixed = TRUE)[[1]]) %in% colnames(model)),
  logical(1))]
if (!length(want))
  stop("None of the requested contrasts are available. Stages present: ",
       paste(colnames(model), collapse = ", "), call. = FALSE)
if (length(want) < length(CONTRASTS))
  message("Note: skipping unavailable contrast(s): ",
          paste(setdiff(CONTRASTS, want), collapse = "; "))

contrs <- limma::makeContrasts(contrasts = want, levels = colnames(model))
fit <- lmat[, rownames(m), drop = FALSE] %>%
  run.de.conc(model, contrs, block = m$Animal.ID)

res.all <- colnames(contrs) %>%
  lapply(format.results, ebayes.output = fit) %>%
  do.call(rbind, .) %>%
  mutate(Stage = factor(Stage, levels = names(cols.stage)))

write.csv(res.all, file.path(ANALYSIS.DIR, "legendplex-allanimals-stage.csv"),
          row.names = FALSE)


## ---- volcano -----------------------------------------------------------------
figS2 <- volcano.plot(res.all, facet = TRUE,
                      n.up = N.LABEL.UP, n.down = N.LABEL.DOWN,
                      xlim = c(-4, 8), xbreaks = c(-4, 0, 4, 8),
                      ylim = c(NA, 32), labels.stage = labels.stage) +
  theme(plot.margin = margin(6, 6, 4, 6))

save.figure(figS2, "FigureS2", FIG.WIDTH, FIG.HEIGHT)
finish.up("figS2")

clear.figure()

message("\nAll figures written to: ", normalizePath(file.path(getwd(), FIG.DIR)))
