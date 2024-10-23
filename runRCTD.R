# R 4.4.1

library(spacexr)
library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# set up single cell reference
ref <- readRDS("single_cell.Rds")

# extract information to pass to the RCTD Reference function
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$new_annotation)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)

# set up query with the RCTD function SpatialRNA
sp <- readRDS("spatialRNA.seurat.bin50.rds")

### read spatial data
# sp <- readRDS(sp_path)
# counts <- slide.seq[["Spatial"]]$counts
counts <- sp@assays$RNA@counts
coords <- sp@meta.data[c('x','y')]
nUMI <- sp@meta.data$nCount_Spatial
names(nUMI) <- rownames(sp@meta.data)
puck <- SpatialRNA(coords, counts, nUMI)

### run RCTD in doublet mode
myRCTD <- create.RCTD(puck, reference, max_cores = 14)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')

#Save the data
saveRDS(myRCTD,"RCTD_output.rds")
save.image("RCTD_workspace.RData")

# Add RCTD results to Seurat object
sp <- AddMetaData(sp, metadata = myRCTD@results$results_df)

# Retrieve the normalize doublet weights in data frame form
norm_weights <- get_doublet_weights(myRCTD)


# Add RCTD results as a new assay
sp[["rctd_full"]] <- CreateAssayObject(data = t(as.matrix(norm_weights)))
if (length(sp@assays$rctd_full@key) == 0) {
  sp@assays$rctd_full@key <- "rctd_full_"
}

# Plotting the tumor cell enrichment weights in the samples
DefaultAssay(sp) <- "rctd_full"
cell_types <- c("Tumor Cells")

FeaturePlot(sp, features = cell_types,reduction="spatial", 
            pt.size=1, cols=c("lightgrey","#ff0000" ))

#Save the image in pdf form
ggsave(file = paste0(f, '.pdf'), device = "pdf", path = "Tumor_weight_featurePlot/",
       width = 11.69, height = 8.27, units = "in")
