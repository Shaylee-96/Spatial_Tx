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