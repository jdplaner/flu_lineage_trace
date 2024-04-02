# annotate by canonical markers
# 20240401
# Joe Planer


# working directory and other addresses with potentially identifiable information have been redacted
# as they would not be useful in reproducing analyses

# dependencies
library(Seurat)
library(ggplot2)
library(tidyverse)

# set wd
# setwd([redacted])

# read in raw data and metadata to be added
scrna<-readRDS("scrna/20240106_merged_data_clustered.rds")

# plot clustering and metadata
DimPlot(scrna, raster =TRUE, raster.dpi = c(2048,2048), label = TRUE, pt.size = 3, group.by= "seurat_clusters") + NoLegend()
#ggsave("annotation/global/20240106_global_dimplot.pdf")

# qc
VlnPlot(scrna, "scrublet_score", pt.size = 0) + NoLegend()
#ggsave("annotation/global/20240106_scrublet_score_global_seuratclusters.pdf")

VlnPlot(scrna, "hybrid_score", pt.size = 0)+ NoLegend()
#ggsave("annotation/global/20240106_hybrid_score_global_seuratclusters.pdf")

#pdf("annotation/global/20240106_scrublet_call_global_clusters.pdf")
barplot(t(table(scrna$seurat_clusters, scrna$scrublet_call)/rowSums(table(scrna$seurat_clusters, scrna$scrublet_call))))
#dev.off()

#pdf("annotation/global/20240106_hybrid_call_global_seurtclusters.pdf")
barplot(t(table(scrna$seurat_clusters, scrna$hybrid_call)/rowSums(table(scrna$seurat_clusters, scrna$hybrid_call))))
#dev.off()

### compartment assignment
# major compartment markers
FeaturePlot(scrna, "Epcam", order = T, raster = TRUE, raster.dpi = c(1024,1024), pt.size = 1)
#ggsave("annotation/global/20240106_epcam.pdf")

FeaturePlot(scrna, "Pecam1", order = T, raster = TRUE, raster.dpi = c(1024,1024), pt.size = 1)
#ggsave("annotation/global/20240106_pecam1.pdf")

FeaturePlot(scrna, "Col1a1", order = T, raster = TRUE, raster.dpi = c(1024,1024), pt.size = 1)
#ggsave("annotation/global/20240106_col1a1.pdf")

FeaturePlot(scrna, "Ptprc", order = T, raster = TRUE, raster.dpi = c(1024,1024), pt.size = 1)
#ggsave("annotation/global/20240106_ptprc.pdf")

DotPlot(scrna, c("Epcam","Pecam1","Col1a1","Ptprc"), assay = "RNA",cols = c("grey80","darkblue")) + scale_y_discrete(limits=rev)
#ggsave("annotation/global/20240106_lineage_markers_dotplot.pdf", height = 8, width = 6)

# annotation into csv file outside of R

# read in lineage assignments
lineage_assignments<-read.csv("annotation/global/20240106_global_lineage_assignments.csv")

# add lineage metadata
scrna$lineage<-scrna$seurat_clusters %>% plyr::mapvalues(from = lineage_assignments$cluster, to = lineage_assignments$lineage)

# lineage colors
lin_cols<-c("#AF2733","#2952AC","#599649","#FDC151","grey80")

DimPlot(scrna, group.by = "lineage", cols = lin_cols)

# generate lineage subsetted objects
#subset(scrna, lineage == "endothelium") %>% saveRDS(., "scrna/lineage/20240107_endothelium.RDS")
#subset(scrna, lineage == "epithelium") %>% saveRDS(., "scrna/lineage/20240107_epithelium.RDS")
#subset(scrna, lineage == "immune") %>% saveRDS(., "scrna/lineage/20240107_immune.RDS")
#subset(scrna, lineage == "mesenchyme") %>% saveRDS(., "scrna/lineage/20240107_mesenchyme.RDS")
#subset(scrna, lineage == "other") %>% saveRDS(., "scrna/lineage/20240107_other_lineage.RDS")

