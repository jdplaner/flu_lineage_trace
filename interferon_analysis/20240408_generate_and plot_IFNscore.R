# generate and add IFN module score to metadata, plot
# 20240408
# Joe Planer

# working directory and other addresses with potentially identifiable information have been redacted
# as they would not be useful in reproducing analyses

# dependencies, wd, seed, date, colors ------------------------------------

# dependencies
library(Seurat)
library(tidyverse)
library(viridis)

# set.seed and wd
set.seed(42)
# setwd([redacted])

# colors
days<-c("darkblue","darkred","firebrick1","maroon","purple3","royalblue", "forestgreen","goldenrod")

# read in raw data and marker list
scrna<-readRDS("scrna/20240108_scrna_main_annotated.RDS")
scrna$lineage<-factor(scrna$lineage, levels = c("Endothelium","Epithelium","Mesenchyme","Myeloid","Lymphoid"))

# find celltype markers
# this removes iMONs since their transcriptional signature is heavily biased toward ISGs and using their marker genes in the blacklist would basically remove ISGs from the analysis
Idents(scrna)<-"celltype"
all_celltype_markers<-FindAllMarkers(scrna, assay = "RNA", min.pct = 0.50, only.pos = T) 
top100_markers<-all_celltype_markers %>% mutate(diff = pct.1 - pct.2) %>% filter( cluster != "iMON") %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 100)

# read in list of ISGs, blacklist cell type marker genes
IFN_stimulated_genes<-#read.csv([redacted])[,2] # see ISG_list.csv in github directory
gene_intersection<-intersect(row.names(scrna), IFN_stimulated_genes)
blacklist<-which(gene_intersection %in% top100_markers$gene)
non_marker_genes<-gene_intersection[-blacklist]

# add module score and plot
scrna<-AddModuleScore(scrna,
                      features = list(non_marker_genes),
                      nbin = 10,
                      ctrl = 100,
                      name = "IFN_stimulated")

scrna@meta.data %>%
  ggplot(., aes(x = lineage, y = IFN_stimulated1, fill = sac_day)) + 
  geom_boxplot(outlier.shape = NA, 
               linewidth = 0.75, 
               alpha = 0.75) + 
  scale_fill_manual(values=days) + 
  theme_classic()
#ggsave("figures/Figure6/20240212_IFN_score_boxplots.pdf", height = 6, width = 9)

ggplot() +
  rasterize(geom_point(data = as.data.frame(scrna[['umap']]@cell.embeddings), aes(x = umap_1, y = umap_2, colour = scrna$IFN_stimulated1), size = 0.25), dpi = 300) + 
  scale_color_viridis(option = "A") + 
  coord_equal() +
  theme_void()

FeaturePlot(scrna, "IFN_stimulated1", raster = T, raster.dpi = c(1024,1024), pt.size = 1.25, max.cutoff = "q99", order = T) + scale_color_viridis(option  = "B", begin = 1, end = 0) +
  coord_equal() + theme_void()
#ggsave("figures/Figure6/20240212_IFN_score_umap.pdf", height = 4, width = 6)
