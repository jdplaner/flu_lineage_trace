# recluster capillary endothelium and make descriptive plots
# 20240408
# Joe Planer

# working directory and other addresses with potentially identifiable information have been redacted
# as they would not be useful in reproducing analyses

# dependencies
library(Seurat)
library(tidyverse)

# date, color palettes
date<-paste(format(Sys.Date(), format = "%Y%m%d"))
days<-c("darkblue","darkred","firebrick1","lightcoral","purple3","royalblue","forestgreen", "goldenrod")
endothelial_cols<-c("#EC958B","#EEC76C","#869FDC")

# set wd and seed
# setwd([redacted])
set.seed(42)

# read in scrna and colors, add color metadata to myeloid data
scrna<-read_rds("scrna/20240108_scrna_main_annotated.RDS")
scrna<-subset(scrna, celltype %in% c("CAP1","CAP2")) # subset capillary endothelium
color_lookup<-read_csv("metadata/20240108_color_assignments.csv")
experimental_metadata<-read_csv("metadata/20240106_minimal_metadata.csv")

# functions ---------------------------------------------------------------
recluster_lineage<-function(lin, res,n){
  # lin: cellular lineage specified in the 'lineage' metadata field
  # res: clustering resolution
  # n: number of PCs to include
  
  subset(scrna, lineage == lin)  %>% 
    RunPCA(npcs = n, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:n, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:n, verbose = FALSE) %>%
    FindClusters(resolution = res, verbose = FALSE)
}

SplitMetaPlot <- function(object,feature=NULL,ncol=3,col=NULL){
  if(is.null(col)){
    stop("Please provide color palette")
  }
  
  groups <- object@meta.data %>% pull(!!feature) %>% unique() %>% sort()
  if(length(groups) > 35) {
    stop("Too Many Groups")
  }
  plots <- groups %>% map(~DimPlot(object = object, reduction = "umap",label=F,pt.size = 5, cells.highlight = object@meta.data %>% filter(!!sym(feature)==!!.x) %>% rownames(), raster = TRUE, raster.dpi = c(1024,1024))  + 
                            coord_equal() + 
                            theme_void() + 
                            ggtitle(.x) + 
                            NoLegend()
  ) 
  
  if(!is.null(col)){
    plots <- map2(plots,col[1:length(groups)], ~.x + scale_color_manual(values=c('lightgrey',.y)))
  }
  patchwork::wrap_plots(plots,ncol=ncol)
}

# recluster lineage -----------------------------------------------
endothelium<-recluster_lineage("Endothelium", 0.6,30)
DimPlot(endothelium, label = T, group.by = "seurat_clusters")

endothelium<-PrepSCTFindMarkers(endothelium)
Idents(endothelium)<-"seurat_clusters"
epithelium_markers<-FindAllMarkers(endothelium, only.pos = T, min.pct = 0.5) %>% mutate(diff = pct.1 - pct.2)

top25_markers<-epithelium_markers %>% group_by(cluster) %>% slice_max(order_by = diff, n = 25)
VlnPlot(endothelium, c("nFeature_RNA","nCount_RNA","hybrid_score"), group.by = "seurat_clusters", pt.size = 0)

# annotation --------------------------------------------------------------
# re-annotate based on capillary re-clustering
endothelium_annotation<-read_csv("annotation/sub-lineage/capillary_endothelium/20240130_endothelial_subtype_assignments.csv")
endothelium$celltype<-endothelium$seurat_clusters %>% plyr::mapvalues(from = endothelium_annotation$cluster, to = endothelium_annotation$celltype)
endothelium$subtype<-endothelium$seurat_clusters %>% plyr::mapvalues(from = endothelium_annotation$cluster, to = endothelium_annotation$subtype)
endothelium<-subset(endothelium, celltype %in% c("CAP1","iCAP","CAP2")) # removes doublets and poor quality cells
endothelium$celltype<-factor(endothelium$celltype, levels = c("CAP1","iCAP","CAP2")) 
DimPlot(endothelium, group.by = "subtype", label = T)

# generate even data set
Idents(endothelium)<-"sac_day"
endothelium_even<-subset(endothelium, downsample = 2000)
Idents(endothelium)<-"celltype"

# descriptive plots ----------------------------------------
# dpi metaplots
SplitMetaPlot(endothelium_even, "sac_day",col = days)
#ggsave("figures/Figure5/20240131_endothelium_even_metaplot_downsample3k.pdf")

# dimplot
DimPlot(endothelium, group.by = "celltype", raster = T, raster.dpi = c(1024,1024), pt.size = 2.5, cols = endothelial_cols) + theme_void()
#ggsave("figures/Figure5/20240131_endothelial_celltype_dimplot.pdf", height = 4, width = 6)

# stacked bargraph
endothelium@meta.data %>% group_by(celltype, sac_day) %>% dplyr::summarize(n=n()) %>%
  ggplot(aes(x = sac_day, y = n, fill = celltype)) + 
  geom_bar(colour = "black", linewidth = 0.25, position = "fill", stat='identity') + 
  scale_fill_manual(values = endothelial_cols) +
  labs(x = "Day Post Infection", y = "% of total", fill= "Cell type") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))
#ggsave("figures/Figure5/20240131_endothelium_stacked_barplot_celltype.pdf", height = 4, width = 8)


# trace fraction --------------------------------------------
# trace umap plot
endothelium_tracecallid<-subset(endothelium, trace_call %in% c("Traced","Untraced"))
DimPlot(endothelium_tracecallid, group.by= "trace_call", cols = c("red","blue"),pt.size=3, raster = T, raster.dpi = c(1024,1024)) + coord_equal()+ theme_void() 
#ggsave("figures/Figure5/20240131_endothelial_proliferation_umap.pdf", width = 4, height = 4)

# trace barplot
endothelium@meta.data %>% filter(trace_call %in% c("Traced","Untraced")) %>% group_by(celltype,trace_call) %>% summarize(n = n()) %>%
  ggplot(aes(x = celltype, y = n, fill = trace_call)) + 
  geom_bar(colour = "black", size = 0.5, position = "fill", stat='identity') + 
  scale_fill_manual(values = c("red","blue")) +
  labs(x = "Day Post Infection", y = "% of total", fill= "Cell type") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 0.5, vjust = 0.5, face = "plain"))
#ggsave("figures/Figure5/20240131_endothelial_proliferation_barplot.pdf", width = 6, height = 4)


# plot markers --------------------------------------------
marker_genes<-c("Gpihbp1","Kit","Sparcl1","Ntrk2","H2-Aa","Ifngr1","Car4","Ednrb")

DotPlot(endothelium, group.by ="celltype", marker_genes, assay = "RNA", cols = c("lightgrey","darkblue"), scale = F) + scale_x_discrete(limits=rev) + coord_flip()
#ggsave("figures/Figure5/20240131_capillary_endothelium_celltype_markergene_dotplot.pdf", height = 4, width = 6)

FeaturePlot(endothelium, c("Sparcl1","Ntrk2","H2-Aa","Ifngr1"), cols = c("lightgrey","darkblue"), raster = T, raster.dpi = c(512,512), pt.size = 1)
#ggsave("figures/Figure5/20240131_capillary_endothelium_iCAP_markergene_umap.pdf", height = 4, width = 6)

# save files --------------------------------------------------------------
saveRDS(endothelium, "scrna/lineage/20240208_capillary_endothelium_reclustered.RDS")

