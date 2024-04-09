# monocyte-macrophage transcription factor expression and inflammatory marker plots
# 20240408
# Joe Planer

# working directory and other addresses with potentially identifiable information have been redacted
# as they would not be useful in reproducing analyses

# dependencies
library(Seurat)
library(tidyverse)
library(ggplot2)
library(slingshot)
library(SingleCellExperiment)
library(stringi)
library(tradeSeq)
library(viridis)
library(clusterProfiler)
library(enrichplot)
library(ComplexHeatmap)
library(dorothea)
library(decoupleR)
library(stringr)

# wd, seed, date
# setwd([redacted])
set.seed(42)
date<-paste(format(Sys.Date(), format = "%Y%m%d"))

# colors
days<-c("darkblue","darkred","firebrick1","lightcoral","purple3","royalblue", "forestgreen","goldenrod")
subtype_color_scheme<-as.data.frame(cbind(c("aMAC_a","aMAC_b","iMON_a","iMON_b","iMAC"),c("hotpink","royalblue3","green4","gold","orange")))
colnames(subtype_color_scheme)<-c("subtype","color")

# read in scrna, scale RNA data slot, add color to metadata
myeloid<-read_rds("scrna/lineage/20240109_myeloid_reclustered.RDS")
DefaultAssay(myeloid)<-"RNA"
myeloid<-ScaleData(myeloid)

color_assignments<-read.csv("metadata/20240108_myeloid_specific_color_assignments.csv")
myeloid<-AddMetaData(myeloid,
                     metadata = myeloid$celltype %>% plyr::mapvalues(from = color_assignments$celltype, to = color_assignments$celltype_color),
                     col.name = "celltype_color")
myeloid<-AddMetaData(myeloid,
                     metadata = myeloid$subtype %>% plyr::mapvalues(from = subtype_color_scheme$subtype, to = subtype_color_scheme$color),
                     col.name = "subtype_color")


# set up tf structures ----------------------------------------------------

# animal tfdb v3 list
tf_all<-as.data.frame(read_tsv("../../../../../Sequencing_resources/Mus_musculus_TF_v3.txt")) %>% pull(Symbol)

net <- get_dorothea(organism='mouse', levels=c('A', 'B', 'C'))
net_ab<-net %>% filter(confidence %in% c('A','B') & mor >0 & target %in% row.names(myeloid@assays$RNA@data))

# tf list
tfs<-unique(net_ab$source)
tf_list<-list()

for(i in tfs){
  tf_list[[i]]<-net_ab %>% 
    filter(source == i) %>% 
    pull(target) %>% 
    unique() %>% 
    intersect(., row.names(myeloid@assays$RNA@data))
}

# functions ---------------------------------------------------------------
pst_at<-function(scrna, start, end, subtype_list, threshold=0.1, pval= 0.05, fc=0.5){
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
  counts<-as.matrix(sub@assays$RNA@data)
  counts<-counts[rowSums(counts > 1) > ncol(counts)*threshold, crv@reducedDim %>% rownames()]
  sce <- fitGAM(counts = as.matrix(counts), sds = crv, parallel =F)
  
  # association test
  at<-associationTest(sce)
  
  # return list
  return(list("crv" = crv,
              "scrna" = sub,
              "sce" = sce,
              "at" = at
              ))
}
tf_score<-function(input, tf_db, tf_name){
  AddModuleScore(input, 
                 features = list(tf_db[[tf_name]]),
                 assay = "RNA",
                 nbin = length(tf_db[[tf_name]]),
                 ctrl = 10,
                 name = paste0("tf.",tf_name))
}

# iMON to aMAC differentiation analysis -----------------------------------
# iMON_aMAC trajectory analysis
iMON_aMAC<-pst_at(myeloid, start = "iMON_b", end = "aMAC_a", subtype_list = c("aMAC_a","aMAC_b","iMON_a","iMON_b"), threshold = 0.1)
var_tfs<-iMON_aMAC$at[which(rownames(iMON_aMAC$at) %in% tf_all),]

# generate tf expression heatmap ------------------------------------------
input_genes<-var_tfs %>% arrange(desc(waldStat)) %>% head(.,25) %>% rownames()
pst.ord <- order(iMON_aMAC$crv@curves$Lineage1$lambda, na.last = NA)
heatdata <- as.matrix(iMON_aMAC$scrna@assays$RNA@scale.data[input_genes[1:25], pst.ord])
color_match<-colnames(heatdata) %>% plyr::mapvalues(from = row.names(iMON_aMAC$scrna@meta.data), to = as.character(iMON_aMAC$scrna@meta.data$subtype))
sac_day_match<-colnames(heatdata) %>% plyr::mapvalues(from = row.names(iMON_aMAC$scrna@meta.data), to = as.character(iMON_aMAC$scrna@meta.data$sac_day))

ha <- HeatmapAnnotation(Subtype = color_match,
                        Day = sac_day_match,
                        col = list(Subtype = c("aMAC_a" = "hotpink", "aMAC_b" = "royalblue3", "iMON_a" = "green4","iMON_b" = "gold","iMAC" = "orange"),
                                   Day = c("0" = "darkblue", "6" = "darkred", "11" = "firebrick1",
                                           "19" = "maroon", "25" = "purple3", "42" = "royalblue", "90" = "forestgreen","366" = "goldenrod")))

#pdf(paste0("compartment_specific/myeloid/tf_expression/",date,"_iMON_aMAC_tf_pathway_heatmap_unscaled.pdf"), width = 12, height = 8)
  print(
    Heatmap(heatdata, 
            col = turbo(50),
            cluster_rows = T, 
            cluster_columns = F, 
            show_column_names = F, 
            row_names_gp = grid::gpar(fontsize = 8),
            top_annotation = ha)
  )
#dev.off()


# dot plot
Idents(iMON_aMAC$scrna)<-"subtype"
iMON_aMAC$markers<-FindAllMarkers(iMON_aMAC$scrna, only.pos = T, min.pct = 0.25) %>% mutate(diff = pct.1 - pct.2)

myeloid_marker_genes<-c("Ly6c2","Ifit3","Oasl1","C1qa","Apoe","Il1b","Spp1","Tnf","Siglecf","Krt79","Cidec")
DotPlot(iMON_aMAC$scrna, myeloid_marker_genes,
        assay = "RNA", 
        group.by = 'subtype',
        cols = c("grey","darkblue"), scale = T, scale.by = "size") +
  scale_x_discrete(limits=rev) + 
  scale_y_discrete(limits=rev) + 
  coord_flip()+
  theme(axis.text.x=element_text(angle = 90, hjust = 1))
#ggsave("figures/Figure3/20240129_iMON_aMAC_markergenes_subtype.pdf", height = 6, width = 6)


