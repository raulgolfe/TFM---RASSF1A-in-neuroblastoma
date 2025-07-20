# ------------------------------------------------------------------------------
# 0. INSTALACIÓN Y CARGA DE PAQUETES
# ------------------------------------------------------------------------------

if (!requireNamespace("CellChat", quietly = TRUE)) {
  install.packages("remotes")
  remotes::install_github("jinworks/CellChat")
}

if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("biomaRt")
}
library(biomaRt)
library(CellChat)
library(patchwork)
library(Seurat)
library(dplyr)
library(scran)
library(SingleCellExperiment)
library(Matrix)
library(ggridges)
library(ggplot2)

# ------------------------------------------------------------------------------
# 1. CARGA DE DATOS Y PRIMERAS VISUALIZACIONES
# ------------------------------------------------------------------------------

setwd("PATH_ARCHIVO_RDS")
nbatlas <- readRDS("seuratObj_NBAtlas_share_v20240130.rds")

FeaturePlot(nbatlas, features = "RASSF1", raster=FALSE)
VlnPlot(nbatlas, features = "RASSF1", pt.size = 0, raster = FALSE)
DotPlot(nbatlas, features = "RASSF1")

# ------------------------------------------------------------------------------
# 2. ANÁLISIS DE EXPRESIÓN DE RASSF1 EN CÉLULAS TUMORALES Y ANÁLISIS DIFERENCIAL
# ------------------------------------------------------------------------------

# Subset solo células tumorales
nb_tumor <- subset(nbatlas, idents = "Neuroendocrine")

# Obtener expresión de RASSF1 en tumorales
expr_vals_tumor <- FetchData(nb_tumor, vars = "RASSF1")$RASSF1

# Clasificar por presencia/ausencia de expresión
high_cells <- WhichCells(nb_tumor, expression = RASSF1 > 0)
low_cells <- WhichCells(nb_tumor, expression = RASSF1 == 0)
rassf1_group <- rep(NA, length(expr_vals_tumor))
names(rassf1_group) <- names(expr_vals_tumor)
rassf1_group[high_cells] <- "high"
rassf1_group[low_cells] <- "low"

# Añadir metadata
nb_tumor$RASSF1_group <- rassf1_group
nb_tumor <- subset(nb_tumor, cells = c(high_cells, low_cells))  # elimina las que son NA
Idents(nb_tumor) <- nb_tumor$RASSF1_group

# Análisis diferencial
diff_genes_tumor <- FindMarkers(nb_tumor, ident.1 = "high", ident.2 = "low")
head(diff_genes_tumor)


# ------------------------------------------------------------------------------
# 3. ANÁLISIS SEGÚN RIESGO (DENSIDAD)
# ------------------------------------------------------------------------------

df_high <- data.frame(
  Risk = nb_tumor$Risk[high_cells],
  RASSF1_expr = FetchData(nb_tumor, vars = "RASSF1")[high_cells, 1]
)
df_high_filtered <- subset(df_high, Risk != "unknown")
ggplot(df_high_filtered, aes(x = RASSF1_expr, color = Risk, fill = Risk)) +
  geom_density(alpha = 0.3) +
  theme_minimal() +
  labs(title = "Densidad de expresión de RASSF1 por riesgo",
       x = "Expresión de RASSF1", y= "Densidad") +
  theme(plot.title = element_text(hjust = 0.5))

# ------------------------------------------------------------------------------  
# 4. ANÁLISIS DE CICLO CELULAR EN CÉLULAS TUMORALES  
# ------------------------------------------------------------------------------

# Cargar listas de genes de ciclo celular proporcionadas por Seurat
s.genes <- Seurat::cc.genes$s.genes
g2m.genes <- Seurat::cc.genes$g2m.genes

# Ejecutar CellCycleScoring
nb_tumor <- CellCycleScoring(
  object = nb_tumor,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = FALSE
)

# Consultar resultados de fase celular
table(nb_tumor$Phase)
head(nb_tumor@meta.data[, c("S.Score", "G2M.Score", "Phase")])


# ------------------------------------------------------------------------------
# 5. CELLCCHAT: PREPARACIÓN DE METADATA Y CREACIÓN DE OBJETOS
# ------------------------------------------------------------------------------

library(CellChat)
library(Seurat)
library(Matrix)
data(CellChatDB.human)

# Copiar metadata y crear columna agrupadora
meta_cellchat <- nbatlas@meta.data
meta_cellchat$celltype_rassf1 <- meta_cellchat$Cell_type
tumor_cells <- WhichCells(nbatlas, idents = "Neuroendocrine")
tumor_high <- intersect(tumor_cells, WhichCells(nbatlas, expression = RASSF1 > 0))
tumor_low  <- intersect(tumor_cells, WhichCells(nbatlas, expression = RASSF1 == 0))
meta_cellchat[tumor_high, "celltype_rassf1"] <- "Tumor_high"
meta_cellchat[tumor_low,  "celltype_rassf1"] <- "Tumor_low"

cellchat_all <- createCellChat(
  object = nbatlas@assays$RNA@counts,
  meta = meta_cellchat,
  group.by = "celltype_rassf1"
)
cellchat_all@DB <- CellChatDB.human

cellchat_all <- subsetData(cellchat_all)
cellchat_all <- identifyOverExpressedGenes(cellchat_all)
cellchat_all <- identifyOverExpressedInteractions(cellchat_all)
cellchat_all <- computeCommunProb(cellchat_all)
cellchat_all <- filterCommunication(cellchat_all, min.cells = 10)
cellchat_all <- computeCommunProbPathway(cellchat_all)
cellchat_all <- aggregateNet(cellchat_all)


# ------------------------------------------------------------------------------
# 6. CREACIÓN DE OBJETOS CELLCCHAT POR SUBGRUPOS (redundante pero no modificar)
# ------------------------------------------------------------------------------

# Para el grupo Tumor_high vs resto
meta_high <- nbatlas@meta.data
cells_to_keep_high <- setdiff(colnames(nbatlas), tumor_low)
nb_high_vs_rest <- subset(nbatlas, cells = cells_to_keep_high)
meta_high <- nb_high_vs_rest@meta.data
meta_high$celltype_high <- meta_high$Cell_type
tumor_high_subset <- intersect(rownames(meta_high), tumor_high)
meta_high[tumor_high_subset, "celltype_high"] <- "Tumor_high"
cellchat_high <- createCellChat(
  object = nb_high_vs_rest@assays$RNA@counts,
  meta = meta_high,
  group.by = "celltype_high"
)
cellchat_high@DB <- CellChatDB.human
cellchat_high <- subsetData(cellchat_high)
cellchat_high <- identifyOverExpressedGenes(cellchat_high)
cellchat_high <- identifyOverExpressedInteractions(cellchat_high)
cellchat_high <- computeCommunProb(cellchat_high)
cellchat_high <- filterCommunication(cellchat_high, min.cells = 10)
cellchat_high <- computeCommunProbPathway(cellchat_high)
cellchat_high <- aggregateNet(cellchat_high)

# Para el grupo Tumor_low vs resto
meta_low <- nbatlas@meta.data
cells_to_keep_low <- setdiff(colnames(nbatlas), tumor_high)
nb_low_vs_rest <- subset(nbatlas, cells = cells_to_keep_low)
meta_low <- nb_low_vs_rest@meta.data
meta_low$celltype_low <- meta_low$Cell_type
tumor_low_subset <- intersect(rownames(meta_low), tumor_low)
meta_low[tumor_low_subset, "celltype_low"] <- "Tumor_low"
cellchat_low <- createCellChat(
  object = nb_low_vs_rest@assays$RNA@counts,
  meta = meta_low,
  group.by = "celltype_low"
)
cellchat_low@DB <- CellChatDB.human
cellchat_low <- subsetData(cellchat_low)
cellchat_low <- identifyOverExpressedGenes(cellchat_low)
cellchat_low <- identifyOverExpressedInteractions(cellchat_low)
cellchat_low <- computeCommunProb(cellchat_low)
cellchat_low <- filterCommunication(cellchat_low, min.cells = 10)
cellchat_low <- computeCommunProbPathway(cellchat_low)
cellchat_low <- aggregateNet(cellchat_low)

# ------------------------------------------------------------------------------
# 7. ANÁLISIS Y VISUALIZACIÓN DE PATHWAYS Y ROLES
# ------------------------------------------------------------------------------

load("cellchat_high.RData")
load("cellchat_low.RData")

pathways.show <- unique(cellchat_high@netP$pathways)
cellchat_analisis_high <- netAnalysis_computeCentrality(cellchat_high, slot.name = "netP") 
cellchat_analisis_low <- netAnalysis_computeCentrality(cellchat_low, slot.name = "netP") 
netAnalysis_signalingRole_network(cellchat_analisis_high, signaling = pathways.show, width = 15, height = 10)
netAnalysis_signalingRole_network(cellchat_analisis_low, signaling = pathways.show, width = 15, height = 10)

ht1 <- netAnalysis_signalingRole_heatmap(cellchat_analisis_high, pattern = "outgoing",width = 12.5, height = 15)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_analisis_high, pattern = "incoming",width = 12.5, height = 15)
ht1 + ht2
ht1 <- netAnalysis_signalingRole_heatmap(cellchat_analisis_low, pattern = "outgoing",width = 12.5, height = 15)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_analisis_low, pattern = "incoming",width = 12.5, height = 15)
ht1 + ht2

groupSize <- as.numeric(table(cellchat_high@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_high@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_high@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

groupSize <- as.numeric(table(cellchat_low@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_low@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_low@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


# ------------------------------------------------------------------------------
# 8. SUMATORIO DE PROBABILIDADES GLOBALES E INTENSIDAD OUTGOING/INCOMING
# ------------------------------------------------------------------------------

# Para cellchat_analisis_high
prob_high <- cellchat_analisis_high@netP$prob
prob_sum_high <- apply(prob_high, c(1,2), sum)
outgoing_global_high <- rowSums(prob_sum_high)
incoming_global_high <- colSums(prob_sum_high)

# Para cellchat_analisis_low
prob_low <- cellchat_analisis_low@netP$prob
prob_sum_low <- apply(prob_low, c(1,2), sum)
outgoing_global_low <- rowSums(prob_sum_low)
incoming_global_low <- colSums(prob_sum_low)

outgoing_global_high
incoming_global_high
outgoing_global_low
incoming_global_low

df <- data.frame(
  Group = rep(c("Tumor_high", "Tumor_low"), each = 2),
  Interaction = rep(c("Outgoing", "Incoming"), times = 2),
  Value = c(1.685259e-05, 4.750311e-06, 3.208252e-06, 6.226863e-07)
)
ggplot(df, aes(x = Interaction, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(name = "Intensidad", limits = c(0, max(df$Value)*1.1)) +
  labs(title = "Comparación interacción Outgoing e Incoming",
       x = "Dirección de la interacción",
       fill = "Grupo") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
