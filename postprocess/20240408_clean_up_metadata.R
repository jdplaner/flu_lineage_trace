# Clean up metadata for submission and website
# 20240408
# Joe Planer

# working directory and other addresses with potentially identifiable information have been redacted
# as they would not be useful in reproducing analyses

# dependencies
library(Seurat)
library(tidyverse)

# wd, seed, date
# setwd([redacted])
date<-paste(format(Sys.Date(), format = "%Y%m%d"))

# import scrna objects
scrna<-read_rds("scrna/20240108_scrna_main_annotated.RDS")
alveolar_epithelium<-read_rds("scrna/lineage/20240124_alveolar_epithelium_reclustered.RDS")
myeloid_immune<-read_rds("scrna/lineage/20240109_myeloid_reclustered.RDS")
capillary_endothelium<-read_rds("scrna/lineage/20240208_capillary_endothelium_reclustered.RDS")

# clean up metadata
clean_metadata<-function(input){
  input@meta.data$phase<-input@meta.data$Phase
  input@meta.data$sex<-input@meta.data$sex_f
  input@meta.data$tamoxifen_start_day<-input@meta.data$tam_day
  input@meta.data$sacrifice_day<-input@meta.data$sac_day
  input@meta.data$experimental_group<-input@meta.data$exp_group
  input@meta.data$sex<-gsub("1","f",input@meta.data$sex)
  input@meta.data$sex<-gsub("0","m",input@meta.data$sex)
  return(input)
}

# metadata fields to retain
final_metadata<-c("orig.ident","nCount_RNA","nFeature_RNA","percent.mito","phase","tamoxifen_start_day","sacrifice_day","sex","experimental_group","trace_call","lineage","celltype")
final_metadata_lineage<-c("orig.ident","nCount_RNA","nFeature_RNA","percent.mito","phase","tamoxifen_start_day","sacrifice_day","sex","experimental_group","trace_call","lineage","celltype","subtype")

# clean up metadata
scrna<-clean_metadata(scrna)
scrna@meta.data<-scrna@meta.data[,final_metadata]

alveolar_epithelium<-clean_metadata(alveolar_epithelium)
alveolar_epithelium@meta.data<-alveolar_epithelium@meta.data[,final_metadata_lineage]

myeloid_immune<-clean_metadata(myeloid_immune)
myeloid_immune@meta.data<-myeloid_immune@meta.data[,final_metadata_lineage]

capillary_endothelium<-clean_metadata(capillary_endothelium)
capillary_endothelium@meta.data<-capillary_endothelium@meta.data[,final_metadata_lineage]

# add 'subtype' metadata field to scrna to make compatible with website
scrna<-AddMetaData(scrna,metadata = scrna$celltype, col.name = "subtype")

# set levels for celltype and subtype
alveolar_epithelium$subtype<-factor(alveolar_epithelium$subtype, levels = c("AT1_a","AT1_b","AT1_c","AT1_AT2","Alveolar_transitional","AT2_d","AT2_c","AT2_b","AT2_a"))
myeloid_immune$subtype<-factor(myeloid_immune$subtype, levels = c("aMAC_a","aMAC_b","aMAC_proliferating","iMAC","iMON_a","iMON_b","cMON","pMON","cDC1","cDC2","cDC_proliferating", "maDC","pDC","Neutrophil_a","Neutrophil_b","Neutrophil_c"))
capillary_endothelium$subtype<-factor(capillary_endothelium$subtype, levels = c("CAP1_a","CAP1_b","CAP1_c","CAP1_d","iCAP_a","iCAP_b","iCAP_c","CAP2"))

# set idents
Idents(scrna)<-"celltype"
Idents(alveolar_epithelium)<-"celltype"
Idents(myeloid_immune)<-"celltype"
Idents(capillary_endothelium)<-"celltype"

# change umap orientation in alveolar epithelium 
# this allows cells to be rendered on the public facing website in the same orientation as in the paper, where they have been rotated 180 degrees for clarity (which does not change relative spatial positioning)
alveolar_epithelium[['old_umap']]<-alveolar_epithelium[['umap']]
new_umap<-alveolar_epithelium@reductions$umap@cell.embeddings
new_umap[,1]<-alveolar_epithelium@reductions$old_umap@cell.embeddings[,2]*-1
new_umap[,2]<-alveolar_epithelium@reductions$old_umap@cell.embeddings[,1]*-1
alveolar_epithelium[['umap']]<-CreateDimReducObject(embeddings = new_umap, key = "umap_", assay = DefaultAssay(alveolar_epithelium))

# save out files
#write_csv(scrna@meta.data, file = "scrna/submission_materials/20240322_global_ki67trace_metadata.csv")
#saveRDS(alveolar_epithelium, "scrna/submission_materials/20240322_alveolar_epithelium_ki67trace.RDS")
#saveRDS(myeloid_immune, "scrna/submission_materials/20240322_myeloid_immune_ki67trace.RDS")
#saveRDS(capillary_endothelium, "scrna/submission_materials/20240322_capillary_endothelium_ki67trace.RDS")
#saveRDS(scrna, "scrna/submission_materials/20240322_global_ki67trace.RDS")
