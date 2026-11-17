
#!/usr/bin/env Rscript

# ============================================================
# AI PRIMETTA V3.3 — STEP 41
# ESPRESSIONE DIFFERENZIALE SMALL RNA
#
# --audit  verifica file, campioni e pacchetti
# --run    esegue edgeR sull'intera matrice filtrabile
#
# NON modifica GUI, Docker, modelli o dati sorgente.
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1L ||
    !args[1] %in% c("--audit", "--run")) {
  stop("Utilizzo: Rscript script.R --audit|--run")
}

mode <- args[1]
root <- normalizePath(getwd())

count_file <- file.path(
  root, "outputs/prjna419326",
  "PRJNA419326_mirnaseq_srna_count_matrix_filtered.tsv"
)

design_file <- file.path(
  root, "outputs/prjna419326",
  "PRJNA419326_experimental_design_curated.tsv"
)

runs_dir <- file.path(
  root, "target_prediction_v33_vscode/runs"
)

targets <- c(
  "sRNA_00000037",
  "sRNA_00000435"
)

sample_group <- function(genotype, condition) {
  paste0(
    genotype, "_", condition,
    "_sRNA_rep", 1:3
  )
}

samples <- c(
  sample_group("B73", "WW"),
  sample_group("B73", "DS"),
  sample_group("H21", "WW"),
  sample_group("H21", "DS")
)

groups <- factor(
  rep(
    c("B73_WW", "B73_DS",
      "H21_WW", "H21_DS"),
    each = 3
  ),
  levels = c(
    "B73_WW", "B73_DS",
    "H21_WW", "H21_DS"
  )
)

design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)
rownames(design) <- samples

cat(strrep("=", 68), "\n")
cat("AI PRIMETTA V3.3 | STEP 41 - edgeR\n")
cat(strrep("=", 68), "\n")

# ============================================================
# 1. CONTROLLO FILE ORIGINALI
# ============================================================

stopifnot(
  file.exists(count_file),
  file.exists(design_file),
  dir.exists(runs_dir)
)

header <- strsplit(
  readLines(count_file, n = 1L),
  "\t",
  fixed = TRUE
)[[1]]

required <- c(
  "feature_id",
  "sequence",
  "length",
  "total_count",
  "samples_detected",
  samples
)

if (!all(required %in% header)) {
  stop("Colonne mancanti nella matrice count")
}

if (anyDuplicated(header)) {
  stop("Intestazioni duplicate nella matrice")
}

# ============================================================
# 2. VERIFICA DISEGNO SPERIMENTALE
# ============================================================

meta <- read.delim(
  design_file,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

necessary <- c(
  "run_accession",
  "sample_alias",
  "sequencing_type"
)

if (!all(necessary %in% names(meta))) {
  stop("Metadati incompleti")
}

small <- meta[
  meta$sequencing_type ==
    "small_rna_mirna_seq",
  ,
  drop = FALSE
]

if (
  nrow(small) != 12L ||
  anyDuplicated(small$run_accession) ||
  anyDuplicated(small$sample_alias) ||
  !setequal(small$sample_alias, samples)
) {
  stop("Corrispondenza metadati-campioni non valida")
}

if (qr(design)$rank != ncol(design)) {
  stop("Matrice del modello non a rango pieno")
}

cat("[PASS] 12 campioni e disegno 2x2 verificati\n")
cat("[PASS] Contrasti B73, H21 e interazione definiti\n")

# ============================================================
# 3. CONTROLLO PACCHETTI R
# ============================================================

required_packages <- c("edgeR", "statmod")

missing <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing)) {
  cat("[PACCHETTI MANCANTI]",
      paste(missing, collapse = ", "), "\n")
  stop("Installare i pacchetti prima del training statistico")
}

cat(
  "[PASS] edgeR versione:",
  as.character(packageVersion("edgeR")), "\n"
)

cat(
  "[PASS] statmod versione:",
  as.character(packageVersion("statmod")), "\n"
)

# ============================================================
# 4. MODALITA' AUDIT
# ============================================================

if (mode == "--audit") {
  cat("[OK] Audit completato in sola lettura\n")
  cat("[OK] Nessuna analisi differenziale ancora eseguita\n")
  quit(save = "no", status = 0L)
}

# ============================================================
# 5. LETTURA MIRATA DELLA MATRICE
# ============================================================

cat("[LETTURA] Matrice dei conteggi originali\n")

# Carichiamo soltanto feature_id e le 12 colonne count.
# Le altre colonne non sono necessarie per il modello.

classes <- ifelse(
  header == "feature_id",
  "character",
  ifelse(header %in% samples, "integer", "NULL")
)

data <- read.delim(
  count_file,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  colClasses = classes
)

ids <- data$feature_id

if (
  anyNA(ids) ||
  anyDuplicated(ids) ||
  !all(targets %in% ids)
) {
  stop("Identificativi feature non validi")
}

counts <- as.matrix(
  data[, samples, drop = FALSE]
)

rownames(counts) <- ids

rm(data, ids)
invisible(gc())

if (
  anyNA(counts) ||
  any(counts < 0) ||
  any(colSums(counts) == 0)
) {
  stop("Conteggi non validi")
}

cat(
  "[MATRICE] Feature totali:",
  nrow(counts), "\n"
)

# ============================================================
# 6. OGGETTO edgeR E PREFILTRAGGIO
# ============================================================

y <- edgeR::DGEList(
  counts = counts
)

rm(counts)
invisible(gc())

keep <- edgeR::filterByExpr(
  y,
  design = design
)

if (!all(targets %in% rownames(y)[keep])) {
  stop(
    "Una delle due small RNA non supera ",
    "il filtraggio di espressione"
  )
}

cat(
  "[FILTRO] Feature conservate:",
  sum(keep), "\n"
)

y <- y[
  keep, ,
  keep.lib.sizes = FALSE
]

# ============================================================
# 7. NORMALIZZAZIONE TMMwsp
# ============================================================

y <- edgeR::calcNormFactors(
  y,
  method = "TMMwsp"
)

cat("[NORMALIZZAZIONE] TMMwsp completata\n")

print(
  data.frame(
    campione = samples,
    gruppo = as.character(groups),
    lib_size = y$samples$lib.size,
    norm_factor = y$samples$norm.factors
  ),
  row.names = FALSE
)

# ============================================================
# 8. DISPERSIONE E MODELLO QUASI-LIKELIHOOD
# ============================================================

cat("[MODELLO] Stima dispersioni\n")

y <- edgeR::estimateDisp(
  y,
  design = design
)

cat("[MODELLO] Regressione quasi-likelihood robusta\n")

fit <- edgeR::glmQLFit(
  y,
  design = design,
  robust = TRUE
)

# ============================================================
# 9. DEFINIZIONE DEI TRE CONTRASTI
# ============================================================

# Ordine coefficienti:
# B73_WW, B73_DS, H21_WW, H21_DS

contrasts <- list(

  B73_DS_vs_WW = c(
    -1, 1, 0, 0
  ),

  H21_DS_vs_WW = c(
    0, 0, -1, 1
  ),

  INTERAZIONE_H21_vs_B73 = c(
    1, -1, -1, 1
  )
)

# ============================================================
# 10. CALCOLO DEI TEST STATISTICI
# ============================================================

results <- list()

for (contrast_name in names(contrasts)) {

  qlf <- edgeR::glmQLFTest(
    fit,
    contrast = contrasts[[contrast_name]]
  )

  tab <- qlf$table

  q_bh <- p.adjust(
    tab$PValue,
    method = "BH"
  )

  indices <- match(
    targets,
    rownames(tab)
  )

  for (i in seq_along(targets)) {

    index <- indices[i]

    results[[length(results) + 1L]] <- data.frame(

      srna_id = targets[i],

      contrasto = contrast_name,

      log2FC = tab$logFC[index],

      FC_stimato = 2^tab$logFC[index],

      p_value = tab$PValue[index],

      FDR_BH = q_bh[index],

      feature_testate = nrow(tab),

      stringsAsFactors = FALSE
    )
  }
}

result_table <- do.call(
  rbind,
  results
)

# ============================================================
# 11. CREAZIONE OUTPUT
# ============================================================

stamp <- format(
  Sys.time(),
  "%Y%m%dT%H%M%S",
  tz = "UTC"
)

outdir <- file.path(
  runs_dir,
  paste0(
    "STEP41_EDGER_DROUGHT_",
    stamp,
    "_",
    Sys.getpid()
  )
)

dir.create(
  outdir,
  recursive = FALSE,
  showWarnings = FALSE
)

if (!dir.exists(outdir)) {
  stop("Impossibile creare la cartella output")
}

out_tsv <- file.path(
  outdir,
  "CONTRASTI_2_SRNA.tsv"
)

write.table(
  result_table,
  out_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 12. REPORT DELLA NORMALIZZAZIONE
# ============================================================

norm_table <- data.frame(
  sample = samples,
  group = as.character(groups),
  library_size = y$samples$lib.size,
  normalization_factor = y$samples$norm.factors,
  effective_library_size =
    y$samples$lib.size *
    y$samples$norm.factors
)

write.table(
  norm_table,
  file.path(
    outdir,
    "NORMALIZZAZIONE_12_LIBRERIE.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 13. GRAFICO MDS DEI CAMPIONI
# ============================================================

png(
  file.path(
    outdir,
    "MDS_12_CAMPIONI.png"
  ),
  width = 1600,
  height = 1000,
  res = 150
)

palette_groups <- c(
  B73_WW = "#377eb8",
  B73_DS = "#e41a1c",
  H21_WW = "#4daf4a",
  H21_DS = "#984ea3"
)

limma::plotMDS(
  y,
  labels = samples,
  col = palette_groups[
    as.character(groups)
  ],
  cex = 0.65,
  main = "AI PRIMETTA V3.3 - MDS small RNA"
)

legend(
  "topright",
  legend = names(palette_groups),
  col = unname(palette_groups),
  pch = 16,
  cex = 0.85,
  bty = "n"
)

dev.off()

# ============================================================
# 14. RAPPORTO METODOLOGICO
# ============================================================

report_lines <- c(

  "AI PRIMETTA V3.3 - STEP41",

  paste(
    "Data UTC:",
    format(Sys.time(), tz = "UTC")
  ),

  paste(
    "R version:",
    R.version.string
  ),

  paste(
    "edgeR version:",
    as.character(packageVersion("edgeR"))
  ),

  paste(
    "Feature totali:",
    length(keep)
  ),

  paste(
    "Feature testate:",
    sum(keep)
  ),

  "Disegno: 4 gruppi, 3 repliche per gruppo",

  "Normalizzazione: TMMwsp",

  "Modello: binomiale negativo quasi-likelihood robusto",

  "Contrasti: B73 DS-WW; H21 DS-WW; interazione",

  "Correzione: Benjamini-Hochberg per singolo contrasto",

  "Nessuna replica esclusa",

  "Nessun gene bersaglio biologicamente validato",

  "Nessuna classificazione miRNA/siRNA confermata",

  "I risultati sono inferenze statistiche esplorative"
)

writeLines(
  report_lines,
  file.path(
    outdir,
    "RAPPORTO_STEP41.txt"
  )
)

# ============================================================
# 15. RISULTATI TERMINALE
# ============================================================

cat("\n")
cat(strrep("=", 68), "\n")
cat("[RISULTATI] CONTRASTI DIFFERENZIALI\n")
cat(strrep("=", 68), "\n")

print(
  result_table,
  row.names = FALSE,
  digits = 5
)

cat("\n[OUTPUT]", outdir, "\n")
cat("[OK] STEP41 completato\n")
cat("[OK] Nessuna replica esclusa\n")
cat("[OK] GUI, Docker e modelli ML/GNN invariati\n")
