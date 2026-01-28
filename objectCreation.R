library(Seurat)
library(dplyr)
library(plotly)
library(openxlsx)
library(Matrix)
library(patchwork)

setwd('path/to/sample/folder')

# =============================================================================
# LOAD DATA AND CREATE INIITAL OBJECT
# =============================================================================

# Load original data
barcodes <- read.delim("matrix/barcodes.tsv", header = FALSE)
features <- read.delim("matrix/features.tsv", header = FALSE)
matrix <- readMM("matrix/matrix.mtx")

colnames(matrix) <- barcodes[,1]
rownames(matrix) <- features[,2]

# Create object with genes present in at least 2 cells and
# cells that have at least 30 genes
so.pv <- CreateSeuratObject(matrix, min.cells = 2, min.features = 30)
so.pv

so.pv <- NormalizeData(so.pv)
so.pv <- FindVariableFeatures(so.pv)
so.pv <- ScaleData(so.pv)
so.pv <- RunPCA(so.pv, verbose = FALSE)
so.pv <- RunUMAP(so.pv, dims=1:10, n.components = 2, n.neighbors = 20, 
                 min.dist = 0.5, verbose = FALSE, return.model = TRUE)
so.pv <- FindNeighbors(so.pv)
so.pv <- FindClusters(so.pv, resolution = 0.3, algorithm = 4, 
                      random.seed = 42)
DimPlot(so.pv, reduction = 'pca')
DimPlot(so.pv, reduction = 'umap')

# =============================================================================
# EXPORT DATA FOR SCRUBLET
# =============================================================================

# Create output directory
output_dir <- "matrix/scrublet_input"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Extract filtered counts matrix (genes x cells)
counts_matrix <- so.pv@assays$RNA$counts

# Export matrix in Matrix Market format
writeMM(counts_matrix, file = file.path(output_dir, "matrix.mtx"))

# Export barcodes (cell names)
barcodes <- colnames(counts_matrix)
write.table(barcodes, 
            file = file.path(output_dir, "barcodes.tsv"),
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

# Export features (gene names)
features <- rownames(counts_matrix)
write.table(features, 
            file = file.path(output_dir, "features.tsv"),
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

# =============================================================================
# CREATE FINAL SCRUBLET OBJECT
# =============================================================================

## Run Scrublet in Python ##

# Apply Scrublet filtering in R

scrublet_results <- read.csv(scrublet_csv_path, stringsAsFactors = FALSE)

# Convert predicted_doublet to logical if it's character
if (is.character(scrublet_results$predicted_doublet)) {
scrublet_results$predicted_doublet <- as.logical(scrublet_results$predicted_doublet)
}

# Get singlet barcodes from Scrublet
singlet_barcodes <- scrublet_results$barcode[scrublet_results$predicted_doublet == FALSE]

# Find indices of singlets in the current matrix
barcode_matches <- match(singlet_barcodes, barcodes)

# Remove NAs (barcodes in Scrublet but not in current matrix)
valid_matches <- barcode_matches[!is.na(barcode_matches)]

# Filter original matrix to only singlets
matrix_scrublet <- matrix[, valid_matches]
barcodes_scrublet <- barcodes[valid_matches, drop = FALSE]
features_scrublet <- features

colnames(matrix_scrublet) <- barcodes_scrublet
rownames(matrix_scrublet) <- features_scrublet

so.scrublet <- CreateSeuratObject(matrix_scrublet)
so.scrublet

so.scrublet <- NormalizeData(so.scrublet)
so.scrublet <- FindVariableFeatures(so.scrublet)
so.scrublet <- ScaleData(so.scrublet)
so.scrublet <- RunPCA(so.scrublet, verbose = FALSE)
so.scrublet <- RunUMAP(so.scrublet, dims=1:10, n.components = 2, n.neighbors = 20, 
                       min.dist = 0.5, verbose = FALSE, return.model = TRUE)
so.scrublet <- FindNeighbors(so.scrublet)
so.scrublet <- FindClusters(so.scrublet, resolution = 0.3, algorithm = 4, 
                            random.seed = 42)
DimPlot(so.scrublet, reduction = 'pca')
DimPlot(so.scrublet, reduction = 'umap')

# =============================================================================
# CLEAN DATA MANUALLY
# =============================================================================

# Visualize gene and read distribution
VlnPlot(so.pv, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2,
        group.by = 'orig.ident')

# Filter out cells with top 5% of genes/reads
qfeature <- quantile(so.pv$nFeature_RNA, probs = c(0.0, 0.95))
qcount <- quantile(so.pv$nCount_RNA, probs = c(0.0, 0.95))
so.pv.clean <- subset(so.pv, subset = nFeature_RNA > qfeature[1] & 
                        nFeature_RNA < qfeature[2] & nCount_RNA > qcount[1] & 
                        nCount_RNA < qcount[2])

# Visualize filtered distribution
hist(so.pv.clean$nFeature_RNA, breaks = 100, main = "nFeature Distribution")
hist(so.pv.clean$nCount_RNA, breaks = 100, main = "nCount Distribution")

summary(so.pv.clean$nCount_RNA)
summary(so.pv.clean$nFeature_RNA)

# =============================================================================
# CREATE FINAL OBJECT
# =============================================================================

so.final <- NormalizeData(so.pv.clean)
so.final <- FindVariableFeatures(so.final)
so.final <- ScaleData(so.final)
so.final <- RunPCA(so.final, verbose = FALSE)
so.final <- RunUMAP(so.final, dims=1:10, n.components = 2, n.neighbors = 20, 
                    min.dist = 0.5, verbose = FALSE, return.model = TRUE)
so.final <- FindNeighbors(so.final)
so.final <- FindClusters(so.final, resolution = 0.3, algorithm = 4, 
                         random.seed = 42)
DimPlot(so.final, reduction = 'pca')
DimPlot(so.final, reduction = 'umap')

# =============================================================================
# PLOT COMPARISON
# =============================================================================

wrap_plots(
  DimPlot(so.pv, 
          reduction = "umap") +
    NoAxes() +
    ggtitle("No Filter"),
  
  DimPlot(so.scrublet, 
          reduction = "umap") +
    NoAxes() +
    ggtitle("Scrublet"),
  
  DimPlot(so.final, 
          reduction = "umap") +
    NoAxes() +
    ggtitle("Manual"),
  
  FeaturePlot(so.pv, 
              features = 'nCount_RNA') +
    NoAxes() +
    ggtitle("No Filter Count"),
  
  FeaturePlot(so.scrublet, 
              features = 'nCount_RNA') +
    NoAxes() +
    ggtitle("Scrublet Count"),
  
  FeaturePlot(so.final, 
              features = 'nCount_RNA') +
    NoAxes() +
    ggtitle("Manual Count"),
  
  ncol = 3
) + plot_layout(guides = "collect")

VlnPlot(so.final, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, 
        pt.size = 0, group.by = 'orig.ident')

# =============================================================================
# SAVE OBJECT
# =============================================================================

saveRDS(so.final, 
        "./S.O.final.rds")
