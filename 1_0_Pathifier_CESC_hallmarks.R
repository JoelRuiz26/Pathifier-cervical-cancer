################################################################################
# Running PATHIFIER (Drier et al., 2013)
# Versión HALLMARK (MSigDB H)
################################################################################
suppressPackageStartupMessages({
  library(msigdbr); library(dplyr); library(pathifier)
})

setwd("~/Pathifier-cervical-cancer/")

# --- Load expression data (genes in rows; samples in columns) ---
exp.matrix <- readRDS("Data/counts.rds") %>% as.data.frame()

# --- Load sample metadata to derive phenotypes (Normal vs Tumor) ---
metadata <- readRDS("Data/metadata.rds")

# --- Retrieve HALLMARK gene sets from MSigDB ---
# category = "H"  -> Hallmark gene sets
hallmark_df <- msigdbr(
  species  = "Homo sapiens",
  category = "H"
) %>%
  dplyr::select(gs_name, gs_id, gene_symbol) %>%
  dplyr::distinct()

# --- Build 'gene_sets' matrix (GMT-like): row = pathway; col1 = name; col2 = id; col3.. = genes ---
gs_list   <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)    # pathway -> vector of genes
gs_id_map <- tapply(hallmark_df$gs_id, hallmark_df$gs_name, `[`, 1) # pathway -> stable ID (first)
max_len   <- max(lengths(gs_list))                                  # widest set dictates matrix width

gene_sets <- t(vapply(names(gs_list), function(pw) {
  genes <- unique(gs_list[[pw]])
  c(pw, gs_id_map[[pw]], genes, rep(NA_character_, max_len - length(genes)))
}, FUN.VALUE = character(2 + max_len)))
mode(gene_sets) <- "character"

# --- Build 'gs' list exactly as in the original script (each element = 1-column matrix of genes) ---
gs <- vector("list", length = nrow(gene_sets))
for (i in seq_len(nrow(gene_sets))) {
  a <- as.vector(gene_sets[i, 3:ncol(gene_sets)])    # gene columns only
  a <- a[!is.na(a) & a != ""]                        # drop NA/empty
  gs[[i]] <- matrix(a, ncol = 1)                     # 1-column matrix
}

# --- Pathway names and PATHWAYS container (same fields as original) ---
pathwaynames <- as.list(gene_sets[, 1])
PATHWAYS <- list(); PATHWAYS$gs <- gs; PATHWAYS$pathwaynames <- pathwaynames

# =========================
# Prepare data and parameters
# =========================

stopifnot("Gene" %in% colnames(exp.matrix))
rownames(exp.matrix) <- exp.matrix$Gene
exp.matrix <- exp.matrix[, setdiff(colnames(exp.matrix), "Gene"), drop = FALSE]
exp.matrix <- as.matrix(exp.matrix)

stopifnot(all(c("specimenID","sample_type") %in% colnames(metadata)))
sample_map <- setNames(metadata$sample_type, metadata$specimenID)
common_samples <- intersect(colnames(exp.matrix), names(sample_map))
exp.matrix <- exp.matrix[, common_samples, drop = FALSE]
normals <- as.vector(sample_map[common_samples] == "Solid Tissue Normal")

# --- Calculate min_std on normal samples (25th percentile of gene-wise SD) ---
N.exp.matrix <- exp.matrix[, as.logical(normals), drop = FALSE]
if (ncol(N.exp.matrix) == 0L) stop("No normal samples found in metadata.")
rsd <- apply(N.exp.matrix, 1, sd, na.rm = TRUE)
min_std <- as.numeric(quantile(rsd, 0.25, na.rm = TRUE))

# --- Calculate min_exp globally (10th percentile of all expression values) ---
min_exp <- as.numeric(quantile(as.vector(exp.matrix), 0.10, na.rm = TRUE))

# --- Filter low-value genes and floor at min_exp ---
over   <- apply(exp.matrix, 1, function(x) x > min_exp)
G.over <- apply(over, 2, mean)
G.over <- names(G.over)[G.over > 0.10]
exp.matrix <- exp.matrix[G.over, , drop = FALSE]
exp.matrix[exp.matrix < min_exp] <- min_exp

# --- Keep up to 5000 most variable genes (heuristic) ---
gene_var <- apply(exp.matrix, 1, var)
V <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(5000L, length(gene_var)))]
V <- V[!is.na(V)]
exp.matrix <- exp.matrix[V, , drop = FALSE]
genes <- rownames(exp.matrix)
allgenes <- as.vector(rownames(exp.matrix))

# --- Pack into DATASET ---
DATASET <- list(); DATASET$allgenes <- allgenes; DATASET$normals <- normals; DATASET$data <- exp.matrix

# =========================
# Light sanity checks
# =========================
cat("DATASET dimensions (genes x samples): ", nrow(DATASET$data), " x ", ncol(DATASET$data), "\n", sep = "")
cat("Normals / Tumors: ", sum(DATASET$normals), " / ", sum(!DATASET$normals), "\n", sep = "")

genes_in_sets <- unique(unlist(lapply(PATHWAYS$gs, function(m) as.vector(m[,1]))))
overlap_genes <- length(intersect(DATASET$allgenes, genes_in_sets))
cat("Unique genes in gene sets: ", length(genes_in_sets),
    " | Overlap with DATASET$allgenes: ", overlap_genes, "\n", sep = "")

set_sizes_post <- vapply(PATHWAYS$gs, function(m) sum(DATASET$allgenes %in% as.vector(m[,1])), integer(1))
cat("Pathways (total): ", length(PATHWAYS$gs),
    " | Empty (post-overlap): ", sum(set_sizes_post == 0),
    " | Median size (post-overlap): ", median(set_sizes_post), "\n", sep = "")

# =========================
# Run Pathifier
# =========================
PDS <- quantify_pathways_deregulation(
  DATASET$data,
  DATASET$allgenes,
  PATHWAYS$gs,
  PATHWAYS$pathwaynames,
  DATASET$normals, 
  maximize_stability = TRUE,
  attempts = 10,
  logfile = "logfile_hallmark.txt",
  min_std = min_std,
  min_exp = min_exp
)

# =========================
# Clean-up and save
# =========================
rm(gene_sets, exp.matrix, allgenes, DATASET, PATHWAYS, rsd, V, over, G.over, 
   N.exp.matrix, gs, genes, min_exp, min_std, pathwaynames)

save.image("/STORAGE/csbig/jruiz/1_1_PDS_HALLMARK.RData")
message("DONE")
