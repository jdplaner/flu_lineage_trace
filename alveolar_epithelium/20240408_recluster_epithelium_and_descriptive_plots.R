# recluster alveolar epithelium and make descriptive plots
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
alveolar_cols<-c("#FFD252","#EE8537","#341547","#B9282C")

# set wd and seed
# setwd([redacted])
set.seed(42)

# read in scrna and colors, add color metadata to myeloid data
scrna<-read_rds("scrna/20240108_scrna_main_annotated.RDS")
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

# recluster myeloid lineage -----------------------------------------------
epithelium<-recluster_lineage("Epithelium", 1.2,20)
DimPlot(epithelium, label = T, group.by = "seurat_clusters")
#table(epithelium$seurat_clusters, epithelium$sac_day)

epithelium<-PrepSCTFindMarkers(epithelium)
Idents(epithelium)<-"seurat_clusters"
epithelium_markers<-FindAllMarkers(epithelium, only.pos = T, min.pct = 0.5) %>% mutate(diff = pct.1 - pct.2)

# generate even data set
Idents(epithelium)<-"sac_day"
epithelium_even<-subset(epithelium, downsample = 1000)
Idents(epithelium)<-"celltype"

# annotation --------------------------------------------------------------
# re-annotate based on myeloid re-clustering
epithelium_annotation<-read_csv("annotation/sub-lineage/alveolar_epithelium/20240124_epithelial_subtype_assignments.csv")
epithelium$celltype<-epithelium$seurat_clusters %>% plyr::mapvalues(from = epithelium_annotation$cluster, to = epithelium_annotation$celltype)
epithelium$subtype<-epithelium$seurat_clusters %>% plyr::mapvalues(from = epithelium_annotation$cluster, to = epithelium_annotation$subtype)
epithelium<-subset(epithelium, celltype %in% c("AT1","AT2","AT1_AT2","Alveolar_transitional","Krt5","Ciliated","Secretory"))
epithelium$celltype<-factor(epithelium$celltype, levels = c("Ciliated","Secretory","AT2","AT1_AT2","AT1","Alveolar_transitional","Krt5")) # removes doublets and poor quality cells
#DimPlot(epithelium, group.by = "subtype", label = T)

# descriptive plots for epithelium ----------------------------------------
# dpi metaplots
SplitMetaPlot(epithelium_even, "sac_day",col = days)
#ggsave("figures/FigureS4/20240124_epithelium_even_metaplot_downsample1k.pdf")

# dimplot
DimPlot(epithelium, group.by = "celltype", raster = T, raster.dpi = c(1024,1024), pt.size = 2.5, cols = color_lookup %>% filter(lineage=="Epithelium") %>% pull(celltype_color)) + theme_void()
#ggsave("figures/FigureS4/20240124_epithelial_celltype_dimplot.pdf", height = 4, width = 6)

# stacked bargraph
epithelium@meta.data %>% group_by(celltype, sac_day) %>% dplyr::summarize(n=n()) %>%
  ggplot(aes(x = sac_day, y = n, fill = celltype)) + 
  geom_bar(colour = "black", linewidth = 0.25, position = "fill", stat='identity') + 
  scale_fill_manual(values = unique(color_lookup %>% filter(lineage=="Epithelium") %>% pull(celltype_color))) +
  labs(x = "Day Post Infection", y = "% of total", fill= "Cell type") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))
#ggsave("figures/FigureS4/20240124_epithelium_stacked_barplot_celltype.pdf", height = 4, width = 8)


# alveolar epithelial fraction --------------------------------------------
# subset and plot alveolar epithelial compartment
alveolar_epithelium<-subset(epithelium, celltype %in% c("AT1","AT2","Alveolar_transitional","AT1_AT2"))
alveolar_epithelium$celltype<-factor(alveolar_epithelium$celltype, levels = c("AT1","AT2","Alveolar_transitional","AT1_AT2"))
DimPlot(alveolar_epithelium, group.by = "celltype", cols = alveolar_cols, raster =T, raster.dpi = c(1024,1024), pt.size =3)+ theme_void()
#ggsave("figures/Figure4/20240124_alveolar_epithelial_celltypes_dimplot.pdf", width = 4, height = 2)

# trace umap plot
alveolar_epithelium_tracecallid<-subset(alveolar_epithelium, trace_call %in% c("Traced","Untraced"))
DimPlot(alveolar_epithelium_tracecallid, group.by= "trace_call", cols = c("red","blue"),pt.size=3, raster = T, raster.dpi = c(1024,1024)) + coord_equal()+ theme_void() 
#ggsave("figures/Figure4/20240124_epithelial_proliferation_umap.pdf", width = 4, height = 4)

# trace barplot
alveolar_epithelium@meta.data %>% filter(trace_call %in% c("Traced","Untraced")) %>% group_by(celltype,trace_call) %>% summarize(n = n()) %>%
  ggplot(aes(x = celltype, y = n, fill = trace_call)) + 
  geom_bar(colour = "black", linewidth = 0.5, position = "fill", stat='identity') + 
  scale_fill_manual(values = c("red","blue")) +
  labs(x = "Day Post Infection", y = "% of total", fill= "Cell type") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 0.5, vjust = 0.5, face = "plain"))
#ggsave("figures/Figure4/20240124_epithelial_proliferation_barplot.pdf", width = 6, height = 4)

# alveolar epithelium stacked bargraph
alveolar_epithelium@meta.data %>% group_by(celltype, sac_day) %>% dplyr::summarize(n=n()) %>%
  ggplot(aes(x = sac_day, y = n, fill = celltype)) + 
  geom_bar(colour = "black", linewidth = 0.25, position = "fill", stat='identity') + 
  scale_fill_manual(values = alveolar_cols) +
  labs(x = "Day Post Infection", y = "% of total", fill= "Cell type") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))
#ggsave("figures/Figure4/20240124_alveolar_epithelium_stacked_barplot_celltype.pdf", height = 4, width = 8)

# alveolar epithelial intermediates stacked bargraph
AlvTrans_by_library<-as.data.frame(table(alveolar_epithelium$orig.ident, alveolar_epithelium$celltype)/rowSums(table(alveolar_epithelium$orig.ident, alveolar_epithelium$celltype))) %>%
  filter(Var2 %in% c("Alveolar_transitional","AT1_AT2"))
AlvTrans_by_library<-cbind(AlvTrans_by_library, AlvTrans_by_library$Var1 %>% plyr::mapvalues(from = experimental_metadata$orig.ident, to = experimental_metadata$sac_day))
colnames(AlvTrans_by_library)<-c("orig.ident","celltype","freq","sac_day")
AlvTrans_by_library$sac_day<-factor(AlvTrans_by_library$sac_day, levels = c("0","6","11","19","25","42","90","366"))

AlvTrans_by_library %>% group_by(celltype,sac_day) %>% summarise(mean = mean(freq), sd = sd(freq)) %>%
  ggplot(aes(x = sac_day, y = mean, fill = celltype)) +
  geom_bar(stat = "identity", position = position_dodge(), color ="black")+
  geom_errorbar(aes(ymin=mean+sd, ymax=mean-sd), width = 0.4, position= position_dodge(0.9)) +
  scale_fill_manual(values = c("#341547","#B9282C")) +
  xlab(label = "Day post influenza")+
  ylab("% of alveolar epithelial cells")+
  scale_y_continuous(n.breaks = 5, limits = c(0,0.3))+
  theme_classic()+
  theme(axis.text.x=element_text(angle = 0, hjust = 0.5, size = 10),
        axis.text.y=element_text(angle = 0, hjust = 1, size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))+
  NoLegend()
#ggsave("figures/Figure4/20240124_alveolar_epithelium_intermediates_stacked_barplot_celltype.pdf", height = 4, width = 8)
# save files --------------------------------------------------------------
#saveRDS(epithelium, "scrna/lineage/20240124_epithelium_reclustered.RDS")
#saveRDS(alveolar_epithelium, "scrna/lineage/20240124_alveolar_epithelium_reclustered.RDS")

