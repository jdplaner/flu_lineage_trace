# GO analysis iCAP cells
# 20240819
# Joe Planer

# dependencies, date, seed, wd --------------------------------------------

# dependencies
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(clusterProfiler)
library(viridis)
library(enrichplot)
require(DOSE)

# date, seed, wd
date<-paste(format(Sys.Date(), format = "%Y%m%d"))
set.seed(42)
setwd("~/Morrisey Lab Dropbox/Morrisey Lab/Joe Planer/Collaborations/TKN/TKN_Ki67/scRNAseq/20240105/")

# read in color scheme, data, scale, and downsample per identity class --------------------------
# colors
days<-c("darkblue","darkred","firebrick1","lightcoral","purple3","royalblue", "forestgreen","goldenrod")
celltype_color_scheme<-as_tibble(cbind(c("CAP1","iCAP","CAP2"),c("#EC958B","#EEC76C","#869FDC")))
colnames(celltype_color_scheme)<-c("celltype","color")
subtype_color_scheme<-as_tibble(cbind(c("CAP1_a","CAP1_b","CAP1_c","CAP1_d","iCAP_a","iCAP_c","iCAP_b","CAP2"),
                                      c("red4","tomato","lightpink","#7B4A73","#EEC76C","navajowhite3","darkolivegreen3","#869FDC")))
colnames(subtype_color_scheme)<-c("subtype","color")

# read in scrna, scale RNA assay
scrna<-read_rds("scrna/lineage/20240208_capillary_endothelium_reclustered.RDS")
DefaultAssay(scrna)<-"RNA"
scrna<-ScaleData(scrna)

scrna<-AddMetaData(scrna,
                   metadata = scrna$celltype %>% plyr::mapvalues(from = celltype_color_scheme$celltype, to = celltype_color_scheme$color),
                   col.name = "celltype_color")
scrna<-AddMetaData(scrna,
                   metadata = scrna$subtype %>% plyr::mapvalues(from = subtype_color_scheme$subtype, to = subtype_color_scheme$color),
                   col.name = "subtype_color")

scrna$subtype<-factor(scrna$subtype, levels = subtype_color_scheme$subtype)

# subset main cellular subtypes (note that CAP1_c/d are mostly present at 6 dpi and iCAP_c does not express Ntrk2 at high levels)
scrna_main<-subset(scrna, subtype %in% c("CAP1_a","CAP1_b","iCAP_a","iCAP_b","CAP2"))

# identify shared iCAP signature (compares iCAPs to their parent endothelial cell populations)
scrna_main_CAP1<-subset(scrna, subtype %in% c("CAP1_a","CAP1_b","iCAP_a"))
scrna_main_CAP1<-PrepSCTFindMarkers(scrna_main_CAP1)
Idents(scrna_main_CAP1)<-"celltype"
iCAP1_markers<-FindMarkers(scrna_main_CAP1, only.pos = T, ident.1 = "iCAP", ident.2 = "CAP1") %>% mutate(diff = pct.1 - pct.2)

scrna_main_CAP2<-subset(scrna, subtype %in% c("CAP2","iCAP_b"))
scrna_main_CAP2<-PrepSCTFindMarkers(scrna_main_CAP2)
Idents(scrna_main_CAP2)<-"celltype"
iCAP2_markers<-FindMarkers(scrna_main_CAP2, only.pos = T, ident.1 = "iCAP", ident.2 = "CAP2") %>% mutate(diff = pct.1 - pct.2)

iCAP_shared_signature<-intersect(rownames(iCAP1_markers),rownames(iCAP2_markers))

# create iCAP data frame
iCAP_shared_df<-as.data.frame(cbind(iCAP1_markers[iCAP_shared_signature, "avg_log2FC"],
                                    iCAP2_markers[iCAP_shared_signature, "avg_log2FC"]))
colnames(iCAP_shared_df)<-c("iCAP1_fc","iCAP2_fc")
rownames(iCAP_shared_df)<-iCAP_shared_signature
iCAP_shared_df$iCAP1_fc<-as.numeric(iCAP_shared_df$iCAP1_fc)
iCAP_shared_df$iCAP2_fc<-as.numeric(iCAP_shared_df$iCAP2_fc)

# filter data frame and arrange by mean fold change relative to parent cell population
iCAP_shared_df_filtered<-iCAP_shared_df %>% 
  mutate(avg_fc = rowMeans(.)) %>% 
  rownames_to_column("gene") %>%
  filter(str_detect(gene, "Rps|Rpl", negate = T)) %>% # removes ribosomal genes
  arrange(desc(avg_fc))

# GO enrichment analysis
go_enrich <- enrichGO(gene = iCAP_shared_df_filtered$gene,
                      OrgDb = "org.Mm.eg.db", 
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

barplot(go_enrich, showCategory=10, x = "Count") + theme_classic()
#ggsave("revisions_figures/iCAP_GO/20240818_GO_enriched_terms_iCAP.pdf", height = 4, width = 6)

# plot out genes in metabolic and antigen presentation pathways
VlnPlot(scrna_main, 
        go_enrich@result %>% filter(ID %in% c("GO:0006007","GO:0006735","GO:0061621")) %>% pull(geneID) %>% strsplit("/") %>% unlist() %>% unique(), 
        group.by = "subtype", 
        stack = T, flip =T, cols = rep("darkred",7 ))
#ggsave("revisions_figures/iCAP_GO/20240818_GO_metabolic_term_genes_iCAP.pdf", height = 4, width = 6)

VlnPlot(scrna_main, 
        "Cd74",
        group.by = "subtype",
        split.by = "sac_day",
        pt.size = 0,
        cols = days)
#ggsave("revisions_figures/iCAP_GO/20240818_Cd74_endothelium.pdf", height = 4, width = 6)

VlnPlot(scrna_main, 
        "Ciita",
        group.by = "subtype",
        split.by = "sac_day",
        pt.size = 0,
        cols = days)
#ggsave("revisions_figures/iCAP_GO/20240818_Ciita_endothelium.pdf", height = 4, width = 6)

