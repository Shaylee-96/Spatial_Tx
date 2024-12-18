library(spacexr)
library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# set up reference
ref <- readRDS("singleCell.Rds")

# extract information to pass to the RCTD Reference function
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$new_annotation)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)

#Visium data
sp <- readRDS("sample1.rds")

### read Visium spatial data
counts <- sp@assays$Spatial@counts
#coords <- sp@images$slice1@coordinates
coords <-  GetTissueCoordinates(sp)
#coords <- coords[,c('imagerow','imagecol')]
colnames(coords) <- c("x", "y")
coords[is.na(colnames(coords))] <- NULL
nUMI <- colSums(counts) #set total counts per pixel as uMIs
puck <- SpatialRNA(coords, counts, nUMI)
### run RCTD
myRCTD <- create.RCTD(puck, reference, max_cores = 14)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
# saveRDS(myRCTD,file=paste0('./RCTD/',sample,'.RCTDout.rds'))
saveRDS(myRCTD,"RCTD_output/sample1_RCTD_output.rds")
save.image("RCTD_output/sample1_RCTD_workspace.RData")

sp <- AddMetaData(sp, metadata = myRCTD@results$weights)
  
# Normalize weights
results <- myRCTD@results
norm_weights <- normalize_weights(results$weights)

sp[["rctd_full"]] <- CreateAssayObject(data = t(as.matrix(norm_weights)))
if (length(sp@assays$rctd_full@key) == 0) {
  sp@assays$rctd_full@key <- "rctd_full_"
}

# Plotting
DefaultAssay(sp) <- "rctd_full"
cell_types <- c("rctdfull_Tumorcells")
  
FeaturePlot(sp, features = cell_types,reduction="spatial", 
            pt.size=3
              #cols=c("lightgrey","#ff0000")
) + scale_color_gradientn(colours = c("lightgrey","#ff0000"),  limits = c(0, 1)) +
  scale_y_reverse()# To flip the y-axis in order to aligh with H&E images
  
ggsave(file = fname, device = "pdf", path = "../Tumor_weight_featurePlot/",
        width = 11.69, height = 8.27, units = "in")

