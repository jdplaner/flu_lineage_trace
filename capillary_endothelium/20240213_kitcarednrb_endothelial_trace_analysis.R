# kit/car4/ednrb endothelial trace analysis
# 20240213
# Joe Planer

# dependencies, wd, date --------------------------------------------------

# dependencies
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggbeeswarm)
library(raster)

# wd, seed, date
setwd("~/Morrisey Lab Dropbox/Morrisey Lab/Joe Planer/Collaborations/TKN/TKN_Ki67/scRNAseq/20240105//")
date<-paste(format(Sys.Date(), format = "%Y%m%d"))

# import objects and preprocess -------------------------------------------
# import scrna objects
scrna<-read_rds("scrna/KitCar4Ednrb_objects/20231114_kitmcm_dataset_annotated_main.RDS")
endothelium<-read_rds("scrna/KitCar4Ednrb_objects/20231114_endothelium_annotated_main.RDS")
capillary_endothelium<-subset(endothelium, subtype %in% c(unique(endothelium$subtype) %>% sort() %>% magrittr::extract(1:5)))
DefaultAssay(capillary_endothelium)<-"SCT"

# import color assignments
colors<-read_csv("metadata/20240108_color_assignments.csv")
endothelial_celltype_colors<-c("red4","tomato","lightpink","#EEC76C","darkolivegreen3","#869FDC") # "#7B4A73"

# recluster capillary endothelium, find markers
recluster_lineage<-function(lin, res){
  lin %>% 
    RunPCA(npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
    FindClusters(resolution = res, verbose = FALSE)
}

capillary_endothelium<-recluster_lineage(capillary_endothelium, 0.4)
capillary_endothelium<-PrepSCTFindMarkers(capillary_endothelium)

endo_markers<-FindAllMarkers(capillary_endothelium, min.pct = 0.5)
top_endo_markers<-endo_markers %>% mutate(pctdiff = pct.1 - pct.2) %>% group_by(cluster) %>% slice_max(., order_by = pctdiff, n= 25)

# annotate capillary endothelium
endo_lookup<-as.data.frame(cbind(as.numeric(c(0,seq(1:9))),c("Cap1_a","Cap1_b","iCap_a","Cap1_c","other","Cap2","iCap_b","other","other","other")))
colnames(endo_lookup)<-c("cluster","subtype")
capillary_endothelium$subtype<-capillary_endothelium$seurat_clusters %>% plyr::mapvalues(from = endo_lookup$cluster, to = endo_lookup$subtype) 

# id cell types of interest, subset to remove 'other' cells
#endothelium_levels<-c("Cap1_a","Cap1_b","iCap_1","iCap_2","Cap2","Arterial_endothelium","Venous_endothelium","Lymphatic_endothelium")
capillary_endothelium_levels<-c("Cap1_a","Cap1_b","Cap1_c","iCap_a","iCap_b","Cap2")
capillary_endothelium_identified<-subset(capillary_endothelium, subset = subtype %in% capillary_endothelium_levels)
capillary_endothelium_identified$subtype<-factor(capillary_endothelium_identified$subtype, levels = capillary_endothelium_levels)

# add umap loadings
# add umap loadings to metadata
capillary_endothelium_identified<-AddMetaData(capillary_endothelium_identified, metadata = as.numeric(capillary_endothelium_identified@reductions$umap@cell.embeddings[,1]), col.name = "umap_1")
capillary_endothelium_identified<-AddMetaData(capillary_endothelium_identified, metadata = as.numeric(capillary_endothelium_identified@reductions$umap@cell.embeddings[,2]), col.name = "umap_2")

# plot traced cells
plot_trace<-function(lin,cre_driver, levels, maximum){
  tmp_traced<-capillary_endothelium@meta.data %>% 
    filter(trace_call == "Traced", lineage == lin) %>% 
    dplyr::count(orig.ident,cre,subtype) %>% pivot_wider(names_from = subtype, values_from = n, values_fill = 0) %>% filter(cre == cre_driver)
  tmp_total<-capillary_endothelium@meta.data %>% 
    filter(trace_call %in% c("Traced","Untraced"), lineage == lin) %>% 
    dplyr::count(orig.ident,cre,subtype) %>% pivot_wider(names_from = subtype, values_from = n, values_fill = 0) %>% filter(cre == cre_driver)
  
  tmp_prop<-tmp_traced[,levels]/tmp_total[,levels] * 100
  tmp_prop<- pivot_longer(tmp_prop, cols = levels)
  tmp_prop$name<-factor(tmp_prop$name, levels = levels)
  
  ggplot(data = tmp_prop) +
    geom_bar(data= tmp_prop %>% group_by(name) %>% summarize(mean_se(value)), aes(x = name, y = y), stat= "identity",size = 0.5, color = "black", fill = "white") + 
    geom_errorbar(data= tmp_prop %>% group_by(name) %>% summarize(mean_se(value)), aes(x=name, ymin=ymin, ymax=ymax), colour="black", size=0.5, width = 0.25) +
    geom_beeswarm(aes(x = name, y = value), size = 3.5)+
    stat_summary(aes(x = name, y = value), fun = mean, geom = "crossbar", width = 0.25, size = 0.25)+
    xlab(label = NULL)+
    ylab(paste0("% ",cre_driver," traced cells"))+
    scale_y_continuous(limits = c(0,maximum))+
    theme_classic()+
    theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5))
}

plot_trace("endothelium","car4",capillary_endothelium_levels, 60)
#ggsave("figures/Figure5/20240827_car4_endothelium_trace_percent.pdf", height = 4, width = 8)

plot_trace("endothelium","ednrb",capillary_endothelium_levels,40)
#ggsave("figures/Figure5/20240827_ednrb_endothelium_trace_percent.pdf", height = 4, width = 8)

plot_trace("endothelium","kit",capillary_endothelium_levels,70)
#ggsave("figures/Figure5/20240827_kit_endothelium_trace_percent.pdf", height = 4, width = 8)

ggplot() +
  rasterize(geom_point(data = capillary_endothelium_identified@meta.data, aes(x = umap_1, y = umap_2), colour = "grey80", size=1), dpi = 300) + 
  rasterize(geom_point(data = capillary_endothelium_identified@meta.data %>% filter(cre == "car4", trace_call == "Traced"), aes(x = umap_1, y = umap_2), colour = "darkgreen", size = 2.5), dpi = 300) + 
  coord_equal() +
  theme_void()
ggsave("figures/Figure5/20240213_car4_endothelial_trace_dimplot.pdf", height = 4, width = 6)

ggplot() +
  rasterize(geom_point(data = capillary_endothelium_identified@meta.data, aes(x = umap_1, y = umap_2), colour = "grey80", size=1), dpi = 300) + 
  rasterize(geom_point(data = capillary_endothelium_identified@meta.data %>% filter(cre == "ednrb", trace_call == "Traced"), aes(x = umap_1, y = umap_2), colour = "darkgreen", size = 2.5), dpi = 300) + 
  coord_equal() +
  theme_void()
ggsave("figures/Figure5/20240213_ednrb_endothelial_trace_dimplot.pdf", height = 4, width = 6)

ggplot() +
  rasterize(geom_point(data = capillary_endothelium_identified@meta.data, aes(x = umap_1, y = umap_2), colour = "grey80", size=1), dpi = 300) + 
  rasterize(geom_point(data = capillary_endothelium_identified@meta.data %>% filter(cre == "kit", trace_call == "Traced"), aes(x = umap_1, y = umap_2), colour = "darkorange", size = 2.5), dpi = 300) + 
  coord_equal() +
  theme_void()
ggsave("figures/Figure5/20240213_kit_endothelial_trace_dimplot.pdf", height = 4, width = 6)



# plot canonical and iCAP markers 
VlnPlot(capillary_endothelium_identified, "Gpihbp1", group.by= "subtype", assay = "RNA", pt.size = 0, cols = endothelial_celltype_colors)+ NoLegend() +
  theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1))
#ggsave("figures/kitcarednrb/endothelial/20231213_kitcarednrb_endothelial_Gpihbp1.pdf", height = 4, width = 6) 

VlnPlot(capillary_endothelium_identified, "Ednrb", group.by= "subtype", assay = "RNA", pt.size = 0, cols = endothelial_celltype_colors)+ NoLegend()+
  theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1))
#ggsave("figures/kitcarednrb/endothelial/20231213_kitcarednrb_endothelial_Ednrb.pdf", height = 4, width = 6) 

VlnPlot(capillary_endothelium_identified, "Sparcl1", group.by= "subtype", assay = "RNA", pt.size = 0, cols = endothelial_celltype_colors)+ NoLegend()+
  theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1))
#ggsave("figures/kitcarednrb/endothelial/20231213_kitcarednrb_endothelial_Sparcl1.pdf", height = 4, width = 6) 

VlnPlot(capillary_endothelium_identified, "Ntrk2", group.by= "subtype", assay = "RNA", pt.size = 0, cols = endothelial_celltype_colors)+ NoLegend()+
  theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1))
#ggsave("figures/kitcarednrb/endothelial/20231213_kitcarednrb_endothelial_Ntrk2.pdf", height = 4, width = 6) 

DimPlot(capillary_endothelium_identified, group.by = "subtype", cols = endothelial_celltype_colors, pt.size = 4, raster= T, raster.dpi = c(1024,1024)) +
  NoLegend() + theme_void()
#ggsave("figures/kitcarednrb/endothelial/20231213_kitcarednrb_endothelial_reclustered_dimplot.pdf", height = 4, width = 6) 

