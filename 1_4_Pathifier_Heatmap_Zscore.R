###############################################################################
# Top-20 KEGG_MEDICUS REFERENCE and VARIANT heatmaps (tumor median PDS z-score)
# Author: Joel Ruiz (template) + small modifications
###############################################################################

options(stringsAsFactors = FALSE)
setwd("~/Pathifier-cervical-cancer/")

###############################################################################
# 1. Load required packages
###############################################################################
suppressPackageStartupMessages({
  if (!requireNamespace("gplots", quietly = TRUE))       install.packages("gplots")
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
  library(gplots)
  library(RColorBrewer)
})

###############################################################################
# 2. Parameters
###############################################################################
K    <- 20L   # number of top pathways to plot
CLIP <- 3.0   # clip z-scores to [-CLIP, +CLIP] for color scale
OUT_DIR <- "Output_plots"

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

###############################################################################
# 3. Load Pathifier results and build PDS matrix
###############################################################################
load("/STORAGE/csbig/jruiz/1_1_PDS.RData")   # object "PDS" and "normals" should exist
stopifnot(exists("PDS"))

# PDS$scores is a list: each element is a numeric vector (samples)
# We bind them into a matrix: rows = pathways, columns = samples.
PDSmatrix <- t(mapply(FUN = c, PDS$scores))

# Make sure rows have names (pathway IDs)
if (is.null(rownames(PDSmatrix))) {
  rownames(PDSmatrix) <- names(PDS$scores)
}

###############################################################################
# 4. Get "normals" phenotype (TRUE = normal, FALSE = tumor)
###############################################################################
# If "normals" was not saved in the RData, rebuild it from metadata
if (!exists("normals") || length(normals) != ncol(PDSmatrix)) {
  stopifnot(exists("metadata"))
  smap <- setNames(metadata$sample_type, metadata$specimenID)
  common <- intersect(colnames(PDSmatrix), names(smap))
  PDSmatrix <- PDSmatrix[, common, drop = FALSE]
  normals <- smap[common] == "Solid Tissue Normal"
}
stopifnot(any(normals), any(!normals))  # at least one normal and one tumor

###############################################################################
# 5. Compute PDS z-scores per pathway across ALL samples
###############################################################################
# Row-wise z-score: for each pathway, subtract mean and divide by sd across samples
zAll <- t(scale(t(PDSmatrix)))
zAll[!is.finite(zAll)] <- NA   # handle NaN/Inf if any

###############################################################################
# 6. Helper function to build a top-K heatmap for a given pathway type
###############################################################################
make_top_heatmap <- function(type_prefix,        # "REFERENCE" or "VARIANT"
                             K,
                             zAll,
                             normals,
                             out_file,
                             main_title) {
  
  # 6.1 Subset rows that belong to the given type
  #     Pathway names look like "KEGG_MEDICUS_REFERENCE_..."
  pattern <- paste0("^KEGG_MEDICUS_", type_prefix, "_")
  sel_rows <- grepl(pattern, rownames(zAll))
  z_sub <- zAll[sel_rows, , drop = FALSE]
  
  if (nrow(z_sub) == 0L) {
    stop("No pathways found for type: ", type_prefix)
  }
  
  # 6.2 Rank pathways by tumor median z-score (directional)
  tumor_z <- z_sub[, !normals, drop = FALSE]
  tumor_median <- apply(tumor_z, 1, median, na.rm = TRUE)
  
  ord  <- order(tumor_median, decreasing = TRUE, na.last = NA)
  keep <- ord[seq_len(min(K, length(ord)))]
  zTop <- z_sub[keep, , drop = FALSE]
  
  # 6.3 Reorder samples: normals first, tumors second
  o_cols <- c(which(normals), which(!normals))
  zTop   <- zTop[, o_cols, drop = FALSE]
  normals_reord <- normals[o_cols]
  
  # 6.4 Clip extreme z-scores for nicer colours
  zTop[zTop >  CLIP] <-  CLIP
  zTop[zTop < -CLIP] <- -CLIP
  
  # 6.5 Row and column clustering
  row.distance <- if (nrow(zTop) > 1) dist(zTop, method = "euclidean") else NULL
  row.cluster  <- if (!is.null(row.distance)) hclust(row.distance, method = "ward.D2") else NULL
  col.distance <- if (ncol(zTop) > 1) dist(t(zTop), method = "euclidean") else NULL
  col.cluster  <- if (!is.null(col.distance)) hclust(col.distance, method = "ward.D2") else NULL
  
  # 6.6 Colour palette and side bar for samples
  my_palette <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(1000))
  my_breaks  <- seq(-CLIP, CLIP, length.out = 1001)
  colLabels  <- ifelse(normals_reord, "#377EB8", "#E41A1C")   # blue = normal, red = tumor
  
  # 6.7 Clean row names: remove "KEGG_MEDICUS_REFERENCE_" or "KEGG_MEDICUS_VARIANT_"
  rn <- rownames(zTop)
  rn <- sub("^KEGG_MEDICUS_REFERENCE_", "", rn)
  rn <- sub("^KEGG_MEDICUS_VARIANT_",   "", rn)
  rownames(zTop) <- rn
  
  # 6.8 Draw heatmap
  png(out_file,
      width = 4200, height = 2200, units = "px",
      res = 300, pointsize = 10)
  
  heatmap.2(zTop,
            main = main_title,
            density.info = "none",
            trace = "none",
            margins = c(10, 40),   # extra space on left/right for row labels
            col = my_palette,
            breaks = my_breaks,
            Rowv = if (!is.null(row.cluster)) as.dendrogram(row.cluster) else FALSE,
            Colv = if (!is.null(col.cluster)) as.dendrogram(col.cluster) else FALSE,
            keysize = 1.0,
            ColSideColors = colLabels,
            na.color = "grey95",
            cexRow = 0.6,         # row label size (pathway names)
            cexCol = 0.6          # column label size (sample IDs)
  )
  
  # Legend for normals vs tumors
  par(lend = 1, xpd = NA)
  legend("topright",
         legend = c("Normals","Tumors"),
         col = c("#377EB8","#E41A1C"),
         lty = 1, lwd = 5, bty = "n", inset = 0.02)
  
  dev.off()
  message("Saved: ", out_file)
}

###############################################################################
# 7. Build heatmaps for REFERENCE and VARIANT (top-20 each)
###############################################################################

make_top_heatmap(
  type_prefix = "REFERENCE",
  K          = K,
  zAll       = zAll,
  normals    = normals,
  out_file   = file.path(OUT_DIR, "heatmap_top20_REFERENCE_tumorMedianZ.png"),
  main_title = "Top 20 KEGG_MEDICUS REFERENCE pathways (tumor median PDS z-score)"
)

make_top_heatmap(
  type_prefix = "VARIANT",
  K          = K,
  zAll       = zAll,
  normals    = normals,
  out_file   = file.path(OUT_DIR, "heatmap_top20_VARIANT_tumorMedianZ.png"),
  main_title = "Top 20 KEGG_MEDICUS VARIANT pathways (tumor median PDS z-score)"
)
