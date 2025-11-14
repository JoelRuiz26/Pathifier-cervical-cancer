###############################################################################
## Heatmap for Pathifier results in R
### Author: Joel Ruiz. Taken and modified from:
#           Angel Garcia-Campos https://github.com/AngelCampos
###############################################################################

###############################################################################
### Installing and/or loading required packages
###############################################################################
setwd("~/Pathifier-cervical-cancer/")
suppressPackageStartupMessages({
  if (!require("gplots"))         { install.packages("gplots", dependencies = TRUE);         library(gplots) }
  if (!require("RColorBrewer"))   { install.packages("RColorBrewer", dependencies = TRUE);   library(RColorBrewer) }
})

###############################################################################
## Load Pathifier results and turn into a matrix
###############################################################################
load("/STORAGE/csbig/jruiz/1_1_PDS.RData")
PDSmatrix <- t(mapply(FUN = c, PDS$scores))

###############################################################################
## Creating Custom Palette
###############################################################################

my_palette <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(n = 1000))

###############################################################################
## Clustering Methods
###############################################################################

row.distance = dist(PDSmatrix, method = "euclidean")
row.cluster  = hclust(row.distance, method = "ward.D2")

col.distance = dist(t(PDSmatrix), method = "euclidean")
col.cluster  = hclust(col.distance, method = "ward.D2")

###############################################################################
## Assign Column labels (Optional)
###############################################################################

colLabels <- as.character(normals)
colLabels[colLabels == "TRUE"]  <- "#377EB8"
colLabels[colLabels == "FALSE"] <- "#E41A1C"

###############################################################################
## Row names (pathway names)
###############################################################################
if (is.null(rownames(PDSmatrix))) rownames(PDSmatrix) <- names(PDS$scores)
rownames(PDSmatrix) <- sub("^KEGG_MEDICUS_", "", rownames(PDSmatrix))

# [CAMBIO] Subsets para KEGG REFERENCE y VARIANT (solo para plotting)
PDS_ref <- PDSmatrix[grepl("^REFERENCE_", rownames(PDSmatrix)), , drop = FALSE]  # [CAMBIO]
PDS_var <- PDSmatrix[grepl("^VARIANT_",   rownames(PDSmatrix)), , drop = FALSE]  # [CAMBIO]

###############################################################################
## 1) HEATMAP GLOBAL (todas las vías, sin nombres)
###############################################################################

png("Output_plots/1_3_heatmap_GLOBAL.png",
    width = 3400, height = 2500, units = "px",
    res = 300,
    pointsize = 12)

hm <- heatmap.2(PDSmatrix,
                main = "PDS-Heatmap (Global)",
                density.info = "none",
                trace = "none",
                margins = c(3, 3),          # [CAMBIO] márgenes pequeños (sin etiquetas)
                col = my_palette,
                Rowv = as.dendrogram(row.cluster),
                Colv = as.dendrogram(col.cluster),
                keysize = 1.0,
                ColSideColors = colLabels,
                na.color = "grey95",
                labRow = NA,                # [CAMBIO] sin nombres de filas
                labCol = NA)                # [CAMBIO] sin nombres de columnas

par(lend = 1, xpd = NA)
legend("topright",
       legend = c("Normals", "Tumors"),
       col = c("dodgerblue", "firebrick1"),
       lty = 1, lwd = 5, bty = "n", inset = 0.02)

dev.off()

###############################################################################
## 2) HEATMAP KEGG REFERENCE
###############################################################################

# [CAMBIO] Recalcular clustering de filas solo para REFERENCE (columnas son las mismas)
row.distance.ref <- dist(PDS_ref, method = "euclidean")      # [CAMBIO]
row.cluster.ref  <- hclust(row.distance.ref, method = "ward.D2")  # [CAMBIO]

png("Output_plots/1_3_heatmap_REFERENCE.png",
    width = 3400, height = 2500, units = "px",
    res = 300,
    pointsize = 12)

hm_ref <- heatmap.2(PDS_ref,
                    main = "PDS-Heatmap (KEGG REFERENCE)",
                    density.info = "none",
                    trace = "none",
                    margins = c(3, 3),      # [CAMBIO] sin etiquetas
                    col = my_palette,
                    Rowv = as.dendrogram(row.cluster.ref),
                    Colv = as.dendrogram(col.cluster),
                    keysize = 1.0,
                    ColSideColors = colLabels,
                    na.color = "grey95",
                    labRow = NA,            # [CAMBIO]
                    labCol = NA)            # [CAMBIO]

par(lend = 1, xpd = NA)
legend("topright",
       legend = c("Normals", "Tumors"),
       col = c("dodgerblue", "firebrick1"),
       lty = 1, lwd = 5, bty = "n", inset = 0.02)

dev.off()

###############################################################################
## 3) HEATMAP KEGG VARIANT
###############################################################################

# [CAMBIO] Clustering para VARIANT
row.distance.var <- dist(PDS_var, method = "euclidean")      # [CAMBIO]
row.cluster.var  <- hclust(row.distance.var, method = "ward.D2")  # [CAMBIO]

png("Output_plots/1_3_heatmap_VARIANT.png",
    width = 3400, height = 2500, units = "px",
    res = 300,
    pointsize = 12)

hm_var <- heatmap.2(PDS_var,
                    main = "PDS-Heatmap (KEGG VARIANT)",
                    density.info = "none",
                    trace = "none",
                    margins = c(3, 3),      # [CAMBIO]
                    col = my_palette,
                    Rowv = as.dendrogram(row.cluster.var),
                    Colv = as.dendrogram(col.cluster),
                    keysize = 1.0,
                    ColSideColors = colLabels,
                    na.color = "grey95",
                    labRow = NA,            # [CAMBIO]
                    labCol = NA)            # [CAMBIO]

par(lend = 1, xpd = NA)
legend("topright",
       legend = c("Normals", "Tumors"),
       col = c("dodgerblue", "firebrick1"),
       lty = 1, lwd = 5, bty = "n", inset = 0.02)

dev.off()
