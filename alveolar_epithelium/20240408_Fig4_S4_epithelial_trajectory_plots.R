# alveolar epithelial trajectory and subtype analysis
# 20240404
# Joe Planer

# working directory and other addresses with potentially identifiable information have been redacted
# as they would not be useful in reproducing analyses

# dependencies, wd, seed, date, colors ------------------------------------

# dependencies
library(Seurat)
library(tidyverse)
library(dorothea)
library(decoupleR)
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
# setwd([redacted])
set.seed(42)
date<-paste(format(Sys.Date(), format = "%Y%m%d"))

# colors
days<-c("darkblue","darkred","firebrick1","lightcoral","purple3","royalblue", "forestgreen","goldenrod")
celltype_color_scheme<-as_tibble(cbind(c("AT1","AT2","Alveolar_transitional","AT1_AT2"),c("#FFD252","#EE8537","#341547","#B9282C")))
colnames(celltype_color_scheme)<-c("celltype","color")
subtype_color_scheme<-as_tibble(cbind(c("AT1_a","AT1_b","AT1_c","AT1_AT2","Alveolar_transitional","AT2_d","AT2_c","AT2_b","AT2_a"),
                            c("mediumpurple3","orchid","cornflowerblue","grey60","lightgreen","yellow","orange","lightpink","indianred3")))
colnames(subtype_color_scheme)<-c("subtype","color")

# read in scrna, scale RNA assay, and add color metadata ------------------
global<-read_rds("scrna/20240108_scrna_main_annotated.RDS")
scrna<-read_rds("scrna/lineage/20240124_alveolar_epithelium_reclustered.RDS")
DefaultAssay(scrna)<-"RNA"
scrna<-ScaleData(scrna)

scrna<-AddMetaData(scrna,
                     metadata = scrna$celltype %>% plyr::mapvalues(from = celltype_color_scheme$celltype, to = celltype_color_scheme$color),
                     col.name = "celltype_color")
scrna<-AddMetaData(scrna,
                     metadata = scrna$subtype %>% plyr::mapvalues(from = subtype_color_scheme$subtype, to = subtype_color_scheme$color),
                     col.name = "subtype_color")

scrna$subtype<-factor(scrna$subtype, levels = subtype_color_scheme$subtype)

# set up tf structures ----------------------------------------------------

# animal tfdb v3 list (PMCID: PMC6323978)
tf_all<-as.data.frame(read_tsv("../../../../../Sequencing_resources/Mus_musculus_TF_v3.txt")) %>% pull(Symbol)

net <- get_dorothea(organism='mouse', levels=c('A', 'B', 'C'))
net_ab<-net %>% filter(confidence %in% c('A','B') & mor >0 & target %in% row.names(scrna@assays$RNA@data))

# tf list
tfs<-unique(net_ab$source)
tf_list<-list()

for(i in tfs){
  tf_list[[i]]<-net_ab %>% 
    filter(source == i) %>% 
    pull(target) %>% 
    unique() %>% 
    intersect(., row.names(scrna@assays$RNA@data))
}

# functions ---------------------------------------------------------------
pst_analysis<-function(scrna, start, end, subtype_list, threshold=0.1, pval= 0.05, fc=0.5){
  # subset specified cellular subtypes
  sub<-subset(scrna, subset = subtype %in% subtype_list)
  
  # generate pseudotime curve
  lin <- getLineages(data = sub@reductions$umap@cell.embeddings,
                     clusterLabels = sub$subtype, 
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

# marker gene dot plots  ---------------------------
# find markers
scrna<-PrepSCTFindMarkers(scrna)

Idents(scrna)<-"celltype"
celltype_markers<-FindAllMarkers(scrna) %>% mutate(diff = pct.1 - pct.2, only.pos = T, min.pct = 0.25)

Idents(scrna)<-"subtype"
subtypes_markers<-FindAllMarkers(scrna) %>% mutate(diff = pct.1 - pct.2)

AT1s<-subset(scrna, celltype == "AT1")
Idents(AT1s)<-"subtype"
AT1_markers<-FindAllMarkers(AT1s) %>% mutate(diff = pct.1 - pct.2)

marker_genes_celltype<-celltype_markers %>% filter(cluster %in% c("AT1","AT2","Alveolar_transitional")) %>% group_by(cluster) %>% slice_max(order_by = diff, n = 5) %>% pull(gene)
scrna$celltype<-factor(scrna$celltype, levels = c("AT2","AT1_AT2","AT1","Alveolar_transitional"))

DotPlot(scrna, group.by ="celltype", unique(marker_genes_celltype), assay = "RNA", cols = c("lightgrey","darkblue"), scale = T) + coord_flip()
#ggsave("figures/Figure4/20240125_celltype_markergene_dotplot.pdf", height = 6, width = 6)


# plot subtype makeup ~ experimental day ----------------------------------

scrna@meta.data %>% group_by(subtype, sac_day) %>% dplyr::summarize(n=n()) %>%
  ggplot(aes(x = sac_day, y = n, fill = subtype)) + 
  geom_bar(colour = "black", linewidth = 0.25, position = "fill", stat='identity') + 
  scale_fill_manual(values = subtype_color_scheme$color) +
  labs(x = "Day Post Infection", y = "% of total", fill= "Cell type") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))
#ggsave("figures/FigureS4/20240129_stacked_barplot_alveolar_epithelial_subtypes.pdf", height = 4, width = 8)


# generate pseudotime trajectories ----------------------------------------
alveolar_transitional_trajectory<-pst_analysis(scrna, 
                                               start = "AT2_a", end = "AT1_a", 
                                               subtype_list = c("AT1_a","AT1_b","AT1_c","Alveolar_transitional","AT2_d","AT2_c","AT2_b","AT2_a"), 
                                               threshold = 0.1)

AT1_AT2_trajectory<-pst_analysis(scrna, 
                                 start = "AT2_a", end = "AT1_a", 
                                 subtype_list = c("AT1_a","AT1_b","AT1_c","AT1_AT2","AT2_d","AT2_c","AT2_b","AT2_a"), 
                                 threshold = 0.1)

# plot trajectories -------------------------------------------------------
ggplot() +
  rasterize(geom_point(data = as.data.frame(scrna@reductions$umap@cell.embeddings),aes(x = umap_1, y = umap_2), colour = "grey80", size=0.5), dpi = 300) + 
  rasterize(geom_point(data = as.data.frame(alveolar_transitional_trajectory$scrna@reductions$umap@cell.embeddings), aes(x = umap_1, y = umap_2), colour = alveolar_transitional_trajectory$scrna$celltype_color, size = 0.5), dpi = 300) + 
  geom_path(data = as.data.frame(slingCurves(alveolar_transitional_trajectory$crv, as.df = T)) %>% filter(Lineage == 1), aes(x = umap_1, y = umap_2), linewidth = 1) +
  coord_equal() +
  theme_void()
#ggsave(paste0("compartment_specific/alveolar_epithelium/trajectory_analysis/",date,"_","alveolar_transitional_trajectory_line_celltype.pdf"), height = 4, width = 6)

ggplot() +
  rasterize(geom_point(data = as.data.frame(scrna@reductions$umap@cell.embeddings),aes(x = umap_1, y = umap_2), colour = "grey80", size=0.5), dpi = 300) + 
  rasterize(geom_point(data = as.data.frame(AT1_AT2_trajectory$scrna@reductions$umap@cell.embeddings), aes(x = umap_1, y = umap_2), colour = AT1_AT2_trajectory$scrna$celltype_color, size = 0.5), dpi = 300) + 
  geom_path(data = as.data.frame(slingCurves(AT1_AT2_trajectory$crv, as.df = T)) %>% filter(Lineage == 1), aes(x = umap_1, y = umap_2), linewidth = 1) +
  coord_equal() +
  theme_void()
#ggsave(paste0("compartment_specific/alveolar_epithelium/trajectory_analysis/",date,"_","AT1_AT2_trajectory_line_celltype.pdf"), height = 4, width = 6)

ggplot() +
  rasterize(geom_point(data = as.data.frame(scrna@reductions$umap@cell.embeddings), aes(x = umap_1, y = umap_2), colour = scrna$subtype_color, size = 0.5), dpi = 300) + 
  geom_path(data = as.data.frame(slingCurves(alveolar_transitional_trajectory$crv, as.df = T)) %>% filter(Lineage == 1), aes(x = umap_1, y = umap_2), linewidth = 1) +
  coord_equal() +
  theme_void()
#ggsave(paste0("compartment_specific/alveolar_epithelium/trajectory_analysis/",date,"_","AT1_AT2_trajectory_line_allcells_by_subtype.pdf"), height = 4, width = 6)


# merge pseudotime embeddings for both lineages and plot distributions---------------------------
# calculate mean of pseudotime embeddings, add to metadata
pst_embeddings<-merge(as.data.frame(alveolar_transitional_trajectory$crv@curves$Lineage1$lambda) %>% rownames_to_column('cell'),
                      as.data.frame(AT1_AT2_trajectory$crv@curves$Lineage1$lambda) %>% rownames_to_column('cell'), 
                      by = 'cell', all.x = T, all.y = T) %>% column_to_rownames('cell') 

pst_embeddings<-pst_embeddings/as.numeric(apply(pst_embeddings, MARGIN = 2, function(x) max(x, na.rm=TRUE)))
pst_embeddings<-rowMeans(pst_embeddings, na.rm = T)

scrna<-AddMetaData(scrna, 
                   metadata = pst_embeddings[match(names(pst_embeddings), rownames(scrna@meta.data))], 
                   col.name = "mean_pst")

# plot subtype distribution along pseudotime ~ experimental day (downsampled by experimental day normalize curve heights)
Idents(scrna)<-"sac_day"
scrna_even<-subset(scrna, downsample = 350)
ggplot(scrna_even@meta.data, aes(x = mean_pst, y = sac_day, fill = subtype, height = ..count..)) + 
  geom_density_ridges2(alpha = 0.75, stat = "density", scale =1.5, panel_scaling = T)+
  scale_fill_manual(values = subtype_color_scheme %>% filter(subtype %in% unique(scrna_even$subtype)) %>% pull(color)) +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(-0.05,1.05))+
  labs(x="Pseudotime",y = "Experimental Day",fill = "Epithelial subtype")+
  theme_classic(base_size = 12)
#ggsave(paste0("compartment_specific/alveolar_epithelium/trajectory_analysis/",date,"_","alveolar_epithelial_subtype_ridgeplot.pdf"))


# plot variable genes and modules ~ pseudotime -------------------------------------
marker_gene_data<-as.data.frame(cbind(scrna@meta.data, t(scrna[['RNA']]$scale.data[c("Igfbp2","Gpm6a","Nnat","Ager","Agrn","Hspg2"),])))

data_main<-marker_gene_data %>% filter(subtype %in% c("AT1_a","AT1_b","AT1_c","Alveolar_transitional","AT1_AT2","AT2_a","AT2_b","AT2_c","AT2_d"))

c("Igfbp2","Gpm6a","Nnat","Ager","Agrn","Hspg2") %>% map(function(i){
  ggplot() + 
    rasterize(geom_point(data=data_main , aes_string(x = "mean_pst", y = i, color = "subtype"), shape = 16, ,alpha = 1, size = 1), dpi = 300) + 
    #stat_smooth(data=data %>% filter(celltype %in% c("AT1","Alveolar_transitional","AT2")), aes_string(x = "pst_normalized", y = i), colour = "#341547", span = 10)+
    stat_smooth(data=data_main, aes_string(x = "mean_pst", y = i), colour = "black", span = 5)+
    scale_color_manual(values = subtype_color_scheme$color[c(1,2,3,4,5,6,7,8,9)])+
    labs(x = "Pseudotime", y = "Expression", fill= "Alveolar epithelial subtype",) +
    ggtitle(i) +
    theme_classic(base_line_size = 0.5, base_size = 12) +
    NoLegend()
  #ggsave(paste0("compartment_specific/alveolar_epithelium/gene_vs_pst/",date,"_pseudotime_",i,".pdf"), height = 3, width = 4)
})

# plot gene expression programs for AT1-AT2 differentiation programs
collagen4_pathway_genes<-c("Col4a1","Col4a2","Col4a3","Col4a4","Col4a5")
laminin_pathway_genes<-c("Lama3","Lama5","Lamb2","Lamb3","Lamc1","Lamc2")
mhc2_genes<-c("H2-Aa","H2-Ab1","H2-Eb1","H2-Eb2","H2-Ob","H2-DMb1","H2-DMb2","H2-DMa","H2-Oa","Cd74")

gene_sets<-list("Collagen4" = collagen4_pathway_genes,
                "Laminin" = laminin_pathway_genes,
                "MHC-II" = mhc2_genes)

scrna<-AddModuleScore(scrna, features = list(gene_sets$Collagen4),name = "Collagen4_Module", assay = "RNA")
scrna<-AddModuleScore(scrna, features = list(gene_sets$Laminin),name = "Laminin_Module", assay = "RNA")
scrna<-AddModuleScore(scrna, features = list(gene_sets$`MHC-II`),name = "MHCII_Module", assay = "RNA")

data_main<-scrna@meta.data %>% filter(subtype %in% c("AT1_a","AT1_b","AT1_c","Alveolar_transitional","AT1_AT2","AT2_a","AT2_b","AT2_c","AT2_d"))

c("Collagen4_Module1","Laminin_Module1","`MHCII_Module1`") %>% map(function(i){
  ggplot() + 
    rasterize(geom_point(data=data_main , aes_string(x = "mean_pst", y = i, colour = "mean_pst"), alpha = 1, size = 0.5), dpi = 300) + 
    #stat_smooth(data=data %>% filter(celltype %in% c("AT1","Alveolar_transitional","AT2")), aes_string(x = "pst_normalized", y = i), colour = "#341547", span = 10)+
    stat_smooth(data=data_main %>% filter(celltype %in% c("AT1","AT1_AT2","AT2")), aes_string(x = "mean_pst", y = i), colour = "black", span = 0.01)+
    scale_colour_viridis(begin = 0, end = 1, option = "E")+
    labs(x = "Pseudotime", y = "Expression", fill= "Alveolar epithelial subtype",) +
    ggtitle(i) +
    theme_classic(base_line_size = 0.5, base_size = 12) +
    NoLegend()
  #ggsave(paste0("compartment_specific/alveolar_epithelium/gene_vs_pst/",date,"_pseudotime_",i,".pdf"), height = 3, width = 4)
})


data_main<-scrna@meta.data %>% mutate(Vegfa = scrna[['RNA']]$data['Vegfa',])
data_main %>%
ggplot() + 
  rasterize(geom_point(aes_string(x = "mean_pst", y = "Vegfa", colour = "mean_pst"), alpha = 1, size = 1), dpi = 300) + 
  #stat_smooth(data=data %>% filter(celltype %in% c("AT1","Alveolar_transitional","AT2")), aes_string(x = "pst_normalized", y = i), colour = "#341547", span = 10)+
  stat_smooth(aes_string(x = "mean_pst", y = "Vegfa"), colour = "black", span = 1)+
  scale_colour_viridis(begin = 0, end = 1, option = "E")+
  labs(x = "Pseudotime", y = "Expression", fill= "Alveolar epithelial subtype",) +
  ggtitle("Vegfa") +
  theme_classic(base_line_size = 0.5, base_size = 12) +
  NoLegend()
#ggsave(paste0("compartment_specific/alveolar_epithelium/gene_vs_pst/",date,"_pseudotime_Vegfa.pdf"), height = 3, width = 4)
