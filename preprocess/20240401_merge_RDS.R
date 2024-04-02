# merge individual seurat objects
# 20240401
# Joe Planer

# working directory and other addresses with potentially identifiable information have been redacted
# as they would not be useful in reproducing analyses

# dependencies
library(Seurat)
library(tidyverse)

# setwd
# setwd([redacted])

# read in individual objects
scRNA_113<-read_rds("scRNA-113_H1N1_02to06dpi/Seurat_SoupX/Seurat.RDS")
scRNA_116<-read_rds("scRNA-116_H1N1_07to11dpi/Seurat_SoupX/Seurat.RDS")
scRNA_117<-read_rds("scRNA-117_H1N1_14to19dpi/Seurat_SoupX/Seurat.RDS")
scRNA_118<-read_rds("scRNA-118_H1N1_21to25dpi/Seurat_SoupX/Seurat.RDS")
scRNA_124<-read_rds("scRNA-124_H1N1_02to42dpi/Seurat_SoupX/Seurat.RDS")
scRNA_125<-read_rds("scRNA-125_H1N1_07to42dpi/Seurat_SoupX/Seurat.RDS")
scRNA_126<-read_rds("scRNA-126_H1N1_14to42dpi/Seurat_SoupX/Seurat.RDS")
scRNA_127<-read_rds("scRNA-127_H1N1_21to42dpi/Seurat_SoupX/Seurat.RDS")
scRNA_158<-read_rds("scRNA-158_homeostasis_D42/Seurat_SoupX/Seurat.RDS")
scRNA_162<-read_rds("scRNA-162_H1N1_02to6dpi/Seurat_SoupX/Seurat.RDS")
scRNA_166<-read_rds("scRNA-166_homeostasis_D3/Seurat_SoupX/Seurat.RDS")
scRNA_167<-read_rds("scRNA-167_H1N1_07to11dpi/Seurat_SoupX/Seurat.RDS")
scRNA_170<-read_rds("scRNA-170_H1N1_15to19dpi/Seurat_SoupX/Seurat.RDS")
scRNA_172<-read_rds("scRNA-172_H1N1_21to25dpi/Seurat_SoupX/Seurat.RDS")
scRNA_179<-read_rds("scRNA-179_H1N1_02to42dpi/Seurat_SoupX/Seurat.RDS")
scRNA_180<-read_rds("scRNA-180_H1N1_07to42dpi/Seurat_SoupX/Seurat.RDS")
scRNA_181<-read_rds("scRNA-181_H1N1_15to42dpi/Seurat_SoupX/Seurat.RDS")
scRNA_182<-read_rds("scRNA-182_H1N1_21to42dpi/Seurat_SoupX/Seurat.RDS")
scRNA_236<-read_rds("scRNA-236_H1N1_02to90dpi/Seurat_SoupX/Seurat.RDS")
scRNA_237<-read_rds("scRNA-237_H1N1_07to90dpi/Seurat_SoupX/Seurat.RDS")
scRNA_238<-read_rds("scRNA-238_H1N1_014to90dpi/Seurat_SoupX/Seurat.RDS")
scRNA_239<-read_rds("scRNA-239_H1N1_21_90dpi/Seurat_SoupX/Seurat.RDS")
scRNA_310<-read_rds("EEM-scRNA-310/Seurat_SoupX/Seurat.RDS")
scRNA_311<-read_rds("EEM-scRNA-311/Seurat_SoupX/Seurat.RDS")
scRNA_312<-read_rds("EEM-scRNA-312/Seurat_SoupX/Seurat.RDS")

#

# merge
ki67_combined<-merge(scRNA_113,
                     y= c(scRNA_116,
                          scRNA_117,
                          scRNA_118,
                          scRNA_124,
                          scRNA_125,
                          scRNA_126,
                          scRNA_127,
                          scRNA_158,
                          scRNA_162,
                          scRNA_166,
                          scRNA_167, 
                          scRNA_170,
                          scRNA_172,
                          scRNA_179,
                          scRNA_180,
                          scRNA_181,
                          scRNA_182,
                          scRNA_236,
                          scRNA_237,
                          scRNA_238,
                          scRNA_239,
                          scRNA_310,
                          scRNA_311,
                          scRNA_312),
                     add.cell.ids = c("H1N1_02to06dpi",
                                      "H1N1_07to11dpi",
                                      "H1N1_14to19dpi",
                                      "H1N1_21to25dpi",
                                      "H1N1_02to42dpi",
                                      "H1N1_07to42dpi",
                                      "H1N1_14to42dpi",
                                      "H1N1_21to42dpi",
                                      "homeostasis_d42",
                                      "H1N1_02to06dpi",
                                      "homeostasis_d3",
                                      "H1N1_07to11dpi",
                                      "H1N1_15to19dpi",
                                      "H1N1_21to25dpi",
                                      "H1N1_02to42dpi",
                                      "H1N1_07to42dpi",
                                      "H1N1_15to42dpi",
                                      "H1N1_21to42dpi",
                                      "H1N1_02to90dpi",
                                      "H1N1_07to90dpi",
                                      "H1N1_14to90dpi",
                                      "H1N1_21to90dpi",
                                      "H1N1_07to366dpi",
                                      "H1N1_07to366dpi",
                                      "H1N1_07to366dpi"),
                     project = "H1N1_Ki67")

# save
#write_rds(ki67_combined, [redacted])


### find variable features
# combine objects into list
rds_list<-list(scRNA_113,scRNA_116,scRNA_117,scRNA_118,scRNA_124,
               scRNA_125,scRNA_126,scRNA_127,scRNA_158,scRNA_162,
               scRNA_166,scRNA_167,scRNA_170,scRNA_172,scRNA_179,
               scRNA_180,scRNA_181,scRNA_182,scRNA_236,scRNA_237,
               scRNA_238,scRNA_239,scRNA_310,scRNA_311,scRNA_312)

# find top 4002 variable features
top_features<-SelectIntegrationFeatures(rds_list, nfeatures = 4002, verbose = TRUE)
top_4k_features<-top_features[-which(top_features %in% c("SiteA","SiteB"))]

# write features out in csv
#write_csv(as.data.frame(top_4k_features), [redacted])