# capillary endothelial trajectory and subtype analysis
# 20240408
# Joe Planer

# working directory and other addresses with potentially identifiable information have been redacted
# as they would not be useful in reproducing analyses

# dependencies, wd, seed, date, colors ------------------------------------

# dependencies
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggridges)
library(ggrastr)
library(slingshot)
library(SingleCellExperiment)
library(stringi)
library(tradeSeq)
library(viridis)
library(ComplexHeatmap)

# wd, seed, date
setwd("~/Morrisey Lab Dropbox/Morrisey Lab/Joe Planer/Collaborations/TKN/TKN_Ki67/scRNAseq/20240105/")
set.seed(42)
date<-paste(format(Sys.Date(), format = "%Y%m%d"))


# read in color schemes and scrna object ----------------------------------
celltype_color_scheme<-as_tibble(cbind(c("CAP1","iCAP","CAP2"),c("#EC958B","#EEC76C","#869FDC")))
colnames(celltype_color_scheme)<-c("celltype","color")
subtype_color_scheme<-as_tibble(cbind(c("CAP1_a","CAP1_b","CAP1_c","CAP1_d","iCAP_a","iCAP_c","iCAP_b","CAP2"),
                                      c("red4","tomato","lightpink","#7B4A73","#EEC76C","navajowhite3","darkolivegreen3","#869FDC")))
colnames(subtype_color_scheme)<-c("subtype","color")

# read in scrna
# setwd([redacted])
DefaultAssay(scrna)<-"RNA"
scrna<-ScaleData(scrna)

scrna<-AddMetaData(scrna,
                   metadata = scrna$celltype %>% plyr::mapvalues(from = celltype_color_scheme$celltype, to = celltype_color_scheme$color),
                   col.name = "celltype_color")
scrna<-AddMetaData(scrna,
                   metadata = scrna$subtype %>% plyr::mapvalues(from = subtype_color_scheme$subtype, to = subtype_color_scheme$color),
                   col.name = "subtype_color")

scrna$subtype<-factor(scrna$subtype, levels = subtype_color_scheme$subtype)

# functions ---------------------------------------------------------------
pst_analysis<-function(scrna, start, end, subtype_list, threshold=0.1, pval= 0.05, fc=0.5){
  # subset specified cellular subtypes
  sub<-subset(scrna, subset = subtype %in% subtype_list)
  
  # generate pseudotime curve
  lin <- getLineages(data = sub@reductions$umap@cell.embeddings,
                     clusterLabels = sub$celltype, 
                     start.clus = start, 
                     end.clus = end)
  crv<-SlingshotDataSet(getCurves(lin))
  
  # add pseudotime embeddings to subset data
  sub<-AddMetaData(sub,metadata = crv@curves$Lineage1$lambda/max(crv@curves$Lineage1$lambda),col.name = "pst")
  
  # fit GAM on abundance-filtered data
  #counts<-as.matrix(sub@assays$RNA@data)
  #counts<-counts[rowSums(counts > 1) > ncol(counts)*threshold, crv@reducedDim %>% rownames()]
  #sce <- fitGAM(counts = as.matrix(counts), sds = crv, parallel =F)
  
  # association test
  #at<-associationTest(sce, nPoints = 8)
  
  # return list
  return(list("crv" = crv,
              "scrna" = sub#,
              #"sce" = sce,
              #"at" = at
              ))
}

# plot subtype makeup ~ experimental day ----------------------------------

scrna@meta.data %>% group_by(subtype, sac_day) %>% dplyr::summarize(n=n()) %>%
  ggplot(aes(x = sac_day, y = n, fill = subtype)) + 
  geom_bar(colour = "black", linewidth = 0.25, position = "fill", stat='identity') + 
  scale_fill_manual(values = subtype_color_scheme$color) +
  labs(x = "Day Post Infection", y = "% of total", fill= "Cell type") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))
#ggsave("figures/FigureS5/20240209_stacked_barplot_capillary_endothelial_subtypes.pdf", height = 4, width = 8)

# generate pseudotime trajectories ----------------------------------------

# note that this uses celltype rather than subtypes as the grouping variable
CAP1_CAP2_trajectory<-pst_analysis(scrna, 
                                    start = "CAP1", end = "CAP2", 
                                    subtype_list = c("CAP1_a","CAP1_b","iCAP_a","iCAP_b","CAP2"), 
                                    threshold = 0.1)

# plot trajectories -------------------------------------------------------
ggplot() +
  rasterize(geom_point(data = as.data.frame(scrna@reductions$umap@cell.embeddings),aes(x = umap_1, y = umap_2), colour = "grey90"), dpi = 300) + 
  rasterize(geom_point(data = as.data.frame(CAP1_CAP2_trajectory$scrna@reductions$umap@cell.embeddings), aes(x = umap_1, y = umap_2), colour = CAP1_CAP2_trajectory$scrna$subtype_color), dpi = 300) + 
  geom_path(data = as.data.frame(slingCurves(CAP1_CAP2_trajectory$crv, as.df = T)), aes(x = umap_1, y = umap_2), linewidth = 1) +
  coord_equal() +
  theme_classic()
#ggsave(paste0("compartment_specific/capillary_endothelium/trajectory_analysis/",date,"_","CAP1_CAP2_trajectory_line_subtype.pdf"))

# add pseudotime embeddings to metadata and plot ---------------------------
scrna_CAP1_CAP2_trajectory<-subset(scrna, subtype %in% c("CAP1_a","CAP1_b","iCAP_a","iCAP_b","CAP2"))
scrna_CAP1_CAP2_trajectory<-AddMetaData(scrna_CAP1_CAP2_trajectory, 
                   metadata = CAP1_CAP2_trajectory$crv@curves$Lineage1$lambda/max(CAP1_CAP2_trajectory$crv@curves$Lineage1$lambda), 
                   col.name = "pst")

# plot subtype distribution along pseudotime ~ experimental day (downsampled by experimental day normalize curve heights)
Idents(scrna_CAP1_CAP2_trajectory)<-"sac_day"
scrna_even<-subset(scrna_CAP1_CAP2_trajectory, downsample = 850)
ggplot(scrna_even@meta.data, aes(x = pst, y = sac_day, fill = subtype, height = ..count..)) + 
  geom_density_ridges2(alpha = 0.75, stat = "density", scale =1.5, panel_scaling = T)+
  scale_fill_manual(values = subtype_color_scheme %>% filter(subtype %in% unique(scrna_even$subtype)) %>% pull(color)) +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(-0.05,1.05))+
  labs(x="Pseudotime",y = "Experimental Day",fill = "Endothelial subtype")+
  theme_classic(base_size = 12)
#ggsave(paste0("compartment_specific/capillary_endothelium//trajectory_analysis/",date,"_","capillary_endothelium_subtype_ridgeplot.pdf"))
