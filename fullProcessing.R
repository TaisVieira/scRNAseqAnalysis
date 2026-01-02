library(Seurat)
library(dplyr)
library(plotly)
library(clustree)
library(openxlsx)

setwd('./R analysis/Scripts')
source('./util_funcs.R')
source('./filterEmptyDrops.R')
source('./filterRRNA.R')
source('./filterScrublet.R')
source('./filterGenes.R')
source('./analyzeSlingshot.R')
setwd('./P_vivax/')

########################## FILTERING AND OBJECT CREATION ########################## 

# Run EmptyDrops filtering
emptydrops_results <- run_emptydrops_filtering(
  matrix_path = "matrix/ERR5087439_matrix.mtx",
  barcodes_path = "matrix/ERR5087439_barcodes.tsv",
  features_path = "matrix/ERR5087439_features.tsv",
  lower = 40,
  output_dir = "matrix/matrix_filtered_emptydrops"
)

# Run Scrublet in Python
# Apply Scrublet filtering in R
scrublet_results <- apply_scrublet_filtering(
  matrix_filtered = emptydrops_results$matrix_filtered,
  barcodes_filtered = emptydrops_results$barcodes_filtered,
  features = emptydrops_results$features,
  scrublet_csv_path = "scrublet_doublet_scores.csv",
  output_dir = "matrix/matrix_filtered_scrublet"
)

# Apply quality filters and create Seurat object
quality_results <- apply_quality_filters(
  matrix_final = scrublet_results$matrix_final,
  barcodes_final = scrublet_results$barcodes_final,
  features = scrublet_results$features,
  min_genes = 30,
  min_cells = 2,
  min_umis_per_gene = 2,
  project_name = 'gam_50k',
  output_dir = "matrix/matrix_final"
)

# Print complete summary
print_complete_filtering_summary(
  emptydrops_stats = emptydrops_results$stats,
  rrna_stats = rrna_results$stats,
  scrublet_stats = scrublet_results$stats,
  quality_stats = quality_results$stats
)

# Normalization and dimension reduction
so.spz <- quality_results$seurat_obj
so.spz <- NormalizeData(so.spz)
so.spz <- ScaleData(so.spz)
so.spz <- FindVariableFeatures(so.spz, nfeatures = 1160)
so.spz <- RunPCA(so.spz, npcs=35, verbose = FALSE)
so.spz <- RunUMAP(so.spz, dims=1:10, n.components = 2, n.neighbors = 20, min.dist = 0.5, 
                  verbose = FALSE, return.model = TRUE)
so.spz <- FindNeighbors(so.spz)
so.spz <- FindClusters(so.spz, resolution = 0.5, algorithm = 4, random.seed = 42)
