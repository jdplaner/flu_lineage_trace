# descriptive plots for global level analyses
# 20240402
# Joe Planer

# working directory and other addresses with potentially identifiable information have been redacted
# as they would not be useful in reproducing analyses

# dependencies
library(Seurat)
library(tidyverse)
library(matrixStats)

# set.seed and wd
set.seed(42)
# setwd([redacted])

# colors
days<-c("darkblue","darkred","firebrick1","maroon","purple3","royalblue","forestgreen","goldenrod")

# read in raw data and color assignments
scrna<-readRDS("scrna/20240108_scrna_main_annotated.RDS")
color_assignments<-read.csv("metadata/20240108_color_assignments.csv")

# evenly downsample data set
Idents(scrna)<-"sac_day"
scrna_even<-subset(scrna, downsample= 8000)
Idents(scrna)<-"celltype"

# metaplot function
SplitMetaPlot <- function(object,feature=NULL,ncol=3,col=NULL){
  if(is.null(col)){
    stop("Please provide color palette")
  }
  
  
  groups <- object@meta.data %>% pull(!!feature) %>% unique() %>% sort()
  if(length(groups) > 35) {
    stop("Too Many Groups")
  }
  plots <- groups %>% map(~DimPlot(object = object, reduction = "umap",label=F,pt.size = 4, cells.highlight = object@meta.data %>% filter(!!sym(feature)==!!.x) %>% rownames(), raster = TRUE, raster.dpi = c(1024,1024))  + 
                            coord_equal() + 
                            theme_void() + 
                            ggtitle(.x) + 
                            NoLegend()
  ) 
  
  if(!is.null(col)){
    plots <- map2(plots,col[1:length(groups)], ~.x + scale_color_manual(values=c('grey',.y)))
  }
  patchwork::wrap_plots(plots,ncol=ncol)
}

######################################################################################################################################################

### DimPlots and Metaplot
DimPlot(scrna, group.by = "lineage", pt.size = 1, raster = T, raster.dpi = c(1024,1024), cols = unique(color_assignments$lineage_color))
#ggsave("figures/FigureS2/20240108_global_dimplot_groupby_lineage.pdf")

DimPlot(scrna, group.by = "celltype", pt.size = 1, raster = T, raster.dpi = c(1024,1024), cols = color_assignments$celltype_color)
#ggsave("figures/Figure2/20240108_global_dimplot_groupby_celltype.pdf")

SplitMetaPlot(scrna_even, "sac_day", col = days)
#ggsave("figures/FigureS2/20240108_global_metaplot_by_sacday.pdf")
rm(scrna_even) %>% gc()

######################################################################################################################################################
### Trace plots
# generate trace
scrna_trace<-subset(scrna, trace_call %in% c("Traced","Untraced"))

# ki67 trace
DimPlot(scrna_trace, group.by = "trace_call",cols = c("red","blue"), raster = T, raster.dpi = c(2048,2048), pt.size = 2)
#ggsave("figures/FigureS2/20240108_ki67_tracecall.pdf")

rm(scrna_trace) %>% gc()

######################################################################################################################################################
### Lineage marker plots
# lineage markers
FeaturePlot(scrna, "rna_Ptprc", raster = T, raster.dpi = c(2048,2048), order = T, pt.size = 1.5, cols = c("grey","darkblue"))
#ggsave("figures/FigureS2/20240108_global_featureplot_Ptprc.pdf")

FeaturePlot(scrna, "rna_Col1a1", raster = T, raster.dpi = c(2048,2048), order = T, pt.size = 1.5, cols = c("grey","darkblue"))
#ggsave("figures/FigureS2/20240108_global_featureplot_Col1a1.pdf")

FeaturePlot(scrna, "rna_Epcam", raster = T, raster.dpi = c(2048,2048),order = T, pt.size = 1.5, cols = c("grey","darkblue"))
#ggsave("figures/FigureS2/20240108_global_featureplot_Epcam.pdf")

FeaturePlot(scrna, "rna_Pecam1", raster = T, raster.dpi = c(2048,2048), order = T, pt.size = 1.5, cols = c("grey","darkblue"))
#ggsave("figures/FigureS2/20240108_global_featureplot_Pecam1.pdf")

######################################################################################################################################################
# plot lineage composition by day
scrna$lineage<-factor(scrna$lineage, levels = c("Endothelium","Epithelium","Mesenchyme","Myeloid","Lymphoid"))

scrna@meta.data %>% group_by(lineage, sac_day) %>% summarize(n = n()) %>%
  ggplot(aes(x = sac_day, y = n, fill = lineage)) + 
  geom_bar(colour = "black", size = 0.25, position = "fill", stat='identity') + 
  scale_fill_manual(values = unique(color_assignments$lineage_color)) +
  labs(x = "Day Post Infection", y = "% of total", fill= "Cell type") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))
#ggsave("figures/FigureS2/20240108_lineage_breakdown_by_day.pdf", height = 4, width = 6)

######################################################################################################################################################
### plot celltype markers
# manually curated markers
markers_of_interest<-c("Bmx","Slc6a2","Mmrn1","Gpihbp1","Kit","Ednrb","Car4",
                       "Foxj1","Scgb3a2","Sftpc","Lamp3","Hopx","Cldn4","Krt8","Krt5",
                       "Aspn","Cnn1","Dcn","Pdgfra","Pdgfrb","Msln",
                       "Siglecf","Itgax","C1qb","Cx3cr1","Ly6c2","Adgre4","Flt3","Xcr1","Cd209a","Ccr7","S100a8","Cxcr2",
                       "Cd79a","Ms4a1","Jchain","Cd3e","Nkg7")

# generate dotplot
DefaultAssay(scrna)<-"SCT"
DotPlot(scrna, features = markers_of_interest, group.by = "celltype", cols = c("grey80","darkblue")) + 
  scale_y_discrete(limits=rev) + 
  theme(axis.text.x=element_text(angle = 90, hjust = 1))
#ggsave("figures/Figure2/20240108_celltype_markergenes_global.pdf", height = 8, width = 12)

######################################################################################################################################################
### plot proliferation by day in the short arm at the cell type level

# reset levels
scrna$celltype<-factor(scrna$celltype, levels = c("Arterial_endothelium", "Venous_endothelium","Lymphatic_endothelium","CAP2","CAP1",
                                                            "Ciliated","Secretory","Krt5","AT1_AT2","Alveolar_transitional","AT2","AT1",
                                                            "Adventitial_fibroblast","VSMC","Peribronchial_fibroblast","AF2","AF1","Mesothelium",
                                                            "Neutrophil","cDC1","cDC2","maDC","pMON","cMON","iMON","iMAC","aMAC",
                                                            "B_lymphocyte","Plasma_cell","T_lymphocyte","NK_cell"))

# generate counts data
timecourse_counts<-scrna@meta.data %>% 
  filter(., trace_call %in% c("Traced","Untraced") & sac_day %in% c("0","6","11","19","25")) %>% 
  select(celltype, orig.ident, sac_day, trace_call)

# generate summary data
summary_data<- timecourse_counts %>% 
  group_by(celltype, sac_day, orig.ident) %>% 
  select(orig.ident, celltype, sac_day, trace_call) %>% 
  dplyr::count(trace_call) %>% 
  pivot_wider(., names_from = trace_call, values_from = n, values_fill = 0) %>% 
  mutate(ratio = 100*Traced/(Untraced+Traced)) %>% 
  select(celltype, sac_day, ratio) %>%
  pivot_wider(., names_from = celltype, values_from = ratio)

# change sac_day from factor to numeric
summary_data$sac_day<-as.numeric(as.character(summary_data$sac_day))

# plot capillary endothelium with loess line
summary_data %>% 
  ggplot(., aes(x = as.numeric(sac_day))) +
  geom_smooth(aes(y = CAP1), color = "#EC958B", se = FALSE) +
  geom_smooth(aes(y = CAP2), color = "#869FDC", se = FALSE) +
  #  geom_smooth(aes(y = Mesenchyme), color = "#2952AC", se = FALSE) +
  #  geom_jitter(aes(y = Cap1), color = "#EC958B", width = 0.35, shape =1, size = 2, stroke = 1.5) +
  #  geom_jitter(aes(y = Cap1_2), color = "#7B4A73", width = 0.35, shape =1, size = 2, stroke = 1.5) +
  #  geom_jitter(aes(y = Cap2), color = "#869FDC", width = 0.35, shape =1, size = 2, stroke = 1.5) +
  #  geom_jitter(aes(y = Mesenchyme), color = "#2952AC", width = 0.35, shape =1, size = 2, stroke = 1.5) +
  xlab("Day")+
  ylab("% tdTomato traced cells")+
  scale_y_continuous(n.breaks = 4, limits = c(-1,30))+
  scale_x_continuous(n.breaks = 6, limits = c(-1,27)) +
  theme_classic()
#ggsave("figures/Figure2/20240108_absolute_proliferation_capillary_endothelium_no_datapoints.pdf", height = 3, width = 4)

timecourse_counts %>% filter(celltype %in% c("CAP1","CAP2")) %>%
  group_by(sac_day, celltype) %>% summarise(n=n()) %>%
  ggplot(., aes(x = as.numeric(sac_day), y = n, fill = celltype)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  scale_fill_manual(values = c("#869FDC","#EC958B")) + scale_y_reverse()+
  #  geom_jitter(aes(y = aMAC), color = "#E76154", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  #  geom_jitter(aes(y = iMAC), color = "#ED8946", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  #  geom_jitter(aes(y = iMON), color = "#FFD06E", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  #  geom_jitter(aes(y = pMON), color = "grey50", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  xlab("Day")+
  ylab("# tdTomato traced cells")+
  #  scale_y_continuous(n.breaks = 5, limits = c(-1,30))+
  #  scale_x_continuous(n.breaks = 6, limits = c(-1,27)) +
  theme_classic()
#ggsave("figures/Figure2/20240108_capillary_endothelium_stacked_barplot.pdf", height = 3, width = 4)

# plot alveolar epithelium with loess line
summary_data %>% 
  ggplot(., aes(x = as.numeric(sac_day))) +
  geom_smooth(aes(y = AT1), color = "#FFD252", se = FALSE) +
  geom_smooth(aes(y = Alveolar_transitional), color = "#341547", se = FALSE) +
  geom_smooth(aes(y = AT2), color = "#EE8537", se = FALSE) +
  geom_smooth(aes(y = AT1_AT2), color = "#B9282C", se = FALSE) +
  #  geom_jitter(aes(y = AT1), color = "#FFD252", width = 0.35, shape =1, size = 2, stroke = 1.5) +
  #  geom_jitter(aes(y = Alveolar_transitional), color = "#341547", width = 0.35, shape =1, size = 2, stroke = 1.5) +
  #  geom_jitter(aes(y = AT2), color = "#EE8537", width = 0.35, shape =1, size = 2, stroke = 1.5) +
  #  geom_jitter(aes(y = AT1_AT2), color = "#B9282C", width = 0.35, shape =1, size = 2, stroke = 1.5) +
  xlab("Day")+
  ylab("% tdTomato traced cells")+
  scale_y_continuous(n.breaks = 5, limits = c(-1,90))+
  scale_x_continuous(n.breaks = 6, limits = c(-1,27)) +
  theme_classic()
#ggsave("figures/Figure2/20240108_absolute_proliferation_alveolar_epithelium_no_datapoints.pdf", height = 3, width = 4)

timecourse_counts %>% filter(celltype %in% c("AT1","Alveolar_transitional","AT1_AT2","AT2")) %>%
  group_by(sac_day, celltype) %>% summarise(n=n()) %>%
  ggplot(., aes(x = as.numeric(sac_day), y = n, fill = celltype)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  scale_fill_manual(values = c("#B9282C","#341547","#EE8537","#FFD252")) + scale_y_reverse()+
  #  geom_jitter(aes(y = aMAC), color = "#E76154", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  #  geom_jitter(aes(y = iMAC), color = "#ED8946", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  #  geom_jitter(aes(y = iMON), color = "#FFD06E", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  #  geom_jitter(aes(y = pMON), color = "grey50", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  xlab("Day")+
  ylab("# tdTomato traced cells")+
  #  scale_y_continuous(n.breaks = 5, limits = c(-1,30))+
  #  scale_x_continuous(n.breaks = 6, limits = c(-1,27)) +
  theme_classic()
#ggsave("figures/Figure2/20240108_alveolar_epithelium_stacked_barplot.pdf", height = 3, width = 4)

# plot myeloid immune with loess line
summary_data %>% 
  ggplot(., aes(x = as.numeric(sac_day))) +
  geom_smooth(aes(y = aMAC), color = "#E76154", se = FALSE) +
  geom_smooth(aes(y = iMAC), color = "#ED8946", se = FALSE) +
  geom_smooth(aes(y = iMON), color = "#FFD06E", se = FALSE) +
  geom_smooth(aes(y = pMON), color = "#B75347", se = FALSE) +
  geom_smooth(aes(y = cMON), color = "#75884b", se = FALSE) +
  #  geom_jitter(aes(y = aMAC), color = "#E76154", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  #  geom_jitter(aes(y = iMAC), color = "#ED8946", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  #  geom_jitter(aes(y = iMON), color = "#FFD06E", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  #  geom_jitter(aes(y = pMON), color = "grey50", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  xlab("Day")+
  ylab("% tdTomato traced cells")+
  scale_y_continuous(n.breaks = 5, limits = c(-1,102))+
  scale_x_continuous(n.breaks = 6, limits = c(-1,27)) +
  theme_classic()
#ggsave("figures/Figure2/20240108_absolute_proliferation_myeloid_immune_no_datapoints.pdf", height = 3, width = 4)

timecourse_counts %>% filter(celltype %in% c("aMAC","iMAC","iMON","cMON","pMON")) %>%
  group_by(sac_day, celltype) %>% summarise(n=n()) %>%
  ggplot(., aes(x = as.numeric(sac_day), y = n, fill = celltype)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  scale_fill_manual(values = c("#B75347","#75884b","#FFD06E","#ED8946","#E76154")) + scale_y_reverse()+
  #  geom_jitter(aes(y = aMAC), color = "#E76154", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  #  geom_jitter(aes(y = iMAC), color = "#ED8946", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  #  geom_jitter(aes(y = iMON), color = "#FFD06E", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  #  geom_jitter(aes(y = pMON), color = "grey50", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  xlab("Day")+
  ylab("# tdTomato traced cells")+
  #  scale_y_continuous(n.breaks = 5, limits = c(-1,30))+
  #  scale_x_continuous(n.breaks = 6, limits = c(-1,27)) +
  theme_classic()
#ggsave("figures/Figure2/20240108_myeloid_immune_stacked_barplot.pdf", height = 3, width = 4)

# plot mesenchyme with loess line
summary_data %>% 
  ggplot(., aes(x = as.numeric(sac_day))) +
  geom_smooth(aes(y = AF1), color = "#EEC76C", se = FALSE) +
  geom_smooth(aes(y = AF2), color = "#5C66A8", se = FALSE) +
  geom_smooth(aes(y = Peribronchial_fibroblast), color = "#6E9968", se = FALSE) +
  geom_smooth(aes(y = Adventitial_fibroblast), color = "#97C583", se = FALSE) +
  geom_smooth(aes(y = VSMC), color = "#808FE0", se = FALSE) +
  #  geom_jitter(aes(y = aMAC), color = "#E76154", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  #  geom_jitter(aes(y = iMAC), color = "#ED8946", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  #  geom_jitter(aes(y = iMON), color = "#FFD06E", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  #  geom_jitter(aes(y = pMON), color = "grey50", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  xlab("Day")+
  ylab("% tdTomato traced cells")+
  scale_y_continuous(n.breaks = 5, limits = c(-1,30))+
  scale_x_continuous(n.breaks = 6, limits = c(-1,27)) +
  theme_classic()
#ggsave("figures/Figure2/20240108_absolute_proliferation_mesenchyme_no_datapoints.pdf", height = 3, width = 4)

timecourse_counts %>% filter(celltype %in% c("AF1","AF2","Peribronchial_fibroblast","Adventitial_fibroblast","VSMC")) %>%
  group_by(sac_day, celltype) %>% summarise(n=n()) %>%
  ggplot(., aes(x = as.numeric(sac_day), y = n, fill = celltype)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  scale_fill_manual(values = c("#97C583","#808FE0","#6E9968","#5C66A8","#EEC76C")) + scale_y_reverse()+
  #  geom_jitter(aes(y = aMAC), color = "#E76154", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  #  geom_jitter(aes(y = iMAC), color = "#ED8946", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  #  geom_jitter(aes(y = iMON), color = "#FFD06E", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  #  geom_jitter(aes(y = pMON), color = "grey50", width = 0.35, shape =16, size = 1, stroke = 1.5) +
  xlab("Day")+
  ylab("# tdTomato traced cells")+
  #  scale_y_continuous(n.breaks = 5, limits = c(-1,30))+
  #  scale_x_continuous(n.breaks = 6, limits = c(-1,27)) +
  theme_classic()
#ggsave("figures/Figure2/20240108_mesenchyme_stacked_barplot.pdf", height = 3, width = 4)


######################################################################################################################################################
### plot proliferation by day in the short arm at the lineage level (represented here as a % of maximum for that lineage)

# generate counts data and consolidate Myeloid and Lymphoid immune compartments into 'Immune'
timecourse_counts<-scrna@meta.data %>% 
  filter(., trace_call %in% c("Traced","Untraced") & sac_day %in% c("0","6","11","19","25")) %>% 
  select(lineage, orig.ident, sac_day, trace_call)

timecourse_counts$lineage<-str_replace_all(timecourse_counts$lineage, "Myeloid", "Immune")
timecourse_counts$lineage<-str_replace_all(timecourse_counts$lineage, "Lymphoid", "Immune")

# generate summary data
summary_data<- timecourse_counts %>% 
  group_by(lineage, sac_day, orig.ident) %>% 
  select(orig.ident, lineage, sac_day, trace_call) %>% 
  dplyr::count(trace_call) %>% 
  pivot_wider(., names_from = trace_call, values_from = n) %>% 
  mutate(ratio = 100*Traced/(Untraced+Traced)) %>% 
  select(lineage, sac_day, ratio) %>%
  pivot_wider(., names_from = lineage, values_from = ratio)

# change sac_day from factor to numeric
summary_data$sac_day<-as.numeric(as.character(summary_data$sac_day))

# calculate percent max
summary_data_percentmax<-summary_data %>% group_by(sac_day) %>% summarise(Immune = mean(Immune),
                                                                          Epithelium = mean(Epithelium),
                                                                          Endothelium = mean(Endothelium),
                                                                          Mesenchyme = mean(Mesenchyme)) %>% column_to_rownames("sac_day") %>% as.matrix()

summary_data_percentmax<-t(t(summary_data_percentmax)/colMaxs(summary_data_percentmax)) *100

# plot percent maximum proliferation data with loess line
summary_data_percentmax %>% as.tibble() %>% mutate(.,sac_day = c("0","6","11","19","25"))%>% 
  ggplot(., aes(x = as.numeric(sac_day))) +
  geom_smooth(aes(y = Endothelium), color = "#AF2733") +
  geom_smooth(aes(y = Epithelium), color = "#FDC151") +
  geom_smooth(aes(y = Immune), color = "#599649") +
  geom_smooth(aes(y = Mesenchyme), color = "#2952AC") +
  xlab("Day")+
  ylab("% tdTomato traced cells")+
  scale_y_continuous(n.breaks = 7, limits = c(-1,110))+
  scale_x_continuous(n.breaks = 4, limits = c(0,26)) +
  theme_classic()
#ggsave("figures/Figure2/20240108_percentmax_proliferation_by_single_cell.pdf", width = 6, height = 3)
