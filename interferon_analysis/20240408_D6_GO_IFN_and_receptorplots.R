# identify IFN GO signature in endothelium cells
# 20240408
# Joe Planer


# working directory and other addresses with potentially identifiable information have been redacted
# as they would not be useful in reproducing analyses

# dependencies
library(Seurat)
library(tidyverse)
library(viridis)
library(clusterProfiler)
library(enrichplot)
#library(organism)
require(DOSE)

# date, seed, wd
date<-paste(format(Sys.Date(), format = "%Y%m%d"))
set.seed(42)
# setwd([redacted])

# read in scrna
scrna<-readRDS("scrna/20240108_scrna_main_annotated.RDS")
capillary_endothelium<-read_rds("scrna/lineage/20240208_capillary_endothelium_reclustered.RDS")
Cap1<-subset(capillary_endothelium, celltype == "CAP1")

# find markers
Cap1<-PrepSCTFindMarkers(Cap1)
Idents(Cap1)<-"sac_day"
Cap1_d6_markers<-FindMarkers(Cap1, ident.1 = "6", test.use ="MAST",only.pos = T)

# GO enrichment analysis
go_enrich <- enrichGO(gene = Cap1_d6_markers %>% arrange(desc(avg_log2FC)) %>% rownames(),
                      OrgDb = "org.Mm.eg.db", 
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

barplot(go_enrich, showCategory=5, x = "Count") + theme_classic()
#ggsave("figures/FigureS6/20240212_GO_enriched_terms_d6_Cap1.pdf", height = 4, width = 6)

DotPlot(scrna, c("Ifnar1","Ifnar2","Ifngr1","Ifngr2","Ifnlr1","Il10rb"), group.by ="celltype", assay = "RNA", cols = c("lightgrey","darkorange")) + scale_y_discrete(limits=rev)
#ggsave("figures/FigureS6/20240212_IFN_receptors.pdf", height =6, width = 6)

FeaturePlot(scrna, c("Ifnar1","Ifnar2"), cols = c("lightgrey","darkblue"), raster = T, raster.dpi = c(1024,1024), pt.size = 1, order =T) & theme_void()
#ggsave("figures/FigureS6/20240212_typr1_IFN_receptors_umap.pdf", height =4, width = 8)

FeaturePlot(scrna, c("Ifngr1","Ifngr2"), cols = c("lightgrey","darkred"), raster = T, raster.dpi = c(1024,1024), pt.size = 1, order =T) & theme_void()
#ggsave("figures/FigureS6/20240212_type2_IFN_receptors_umap.pdf", height =4, width = 8)

FeaturePlot(scrna, c("Ifnlr1","Il10rb"), cols = c("lightgrey","darkgreen"), raster = T, raster.dpi = c(1024,1024), pt.size = 1, order =T) & theme_void()
#ggsave("figures/FigureS6/20240212_type3_IFN_receptors_umap.pdf", height =4, width = 8)
