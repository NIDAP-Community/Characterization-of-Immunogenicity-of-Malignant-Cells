# # Preprocessing h5 Data ~ required for downstream steps
SO <- source("Preprocessing.R")

# garbage collection
rm(list = setdiff(ls(), "SO"))
gc()

SO <- SO$value

# Harmony batch correction of samples
SO <- source("Harmony.R")

# Classification by ModuleScores ~ Geneset averages
source("ModuleScore.R")

# Append CytoTRACE results for stemlike malignant cells
cytoTRACE_res <- read.csv("files_for_analysis/CCBR_1026_Cyto_meta.csv")
SO$cytotrace_celltype <- cytoTRACE_res$cytotrace_celltype

Idents(SO) <- SO$cytotrace_celltype
SO_malign <- subset(SO, idents = "High_Stem_Malign","Int_Stem_Malign","Low_Stem_Malign")

source("ExpressionHeatmap.R")

source("Color_by_Gene.R")
colorByGene(object = SO_malign,
            gene = c("ALDH1A1","CD24","PROM1","CD44","EPCAM","NANOG",
                     "POU5F1","SOX2","SOX9","KRT19","DCLK1","SALL4"),
            reduction.type = "tsne")

source("Violin_Plots_by_Metadata.R")

# Pathway enrichment across high stem and differentiated malignant cells
source("l2p_pathway.R")

l2p_path(upregulated = T)
l2p_path(upregulated = F)

# Violin plots, replace genes.of.interest with set of interest
set1 <- c("ALDH1A1","CD24","PROM1","CD44","EPCAM","NANOG",
          "POU5F1","SOX2","SOX9","KRT19","DCLK1","SALL4")

set2 <- c("B2M","HLA-A","HLA-B","HLA-C","HLA-DRA","HLA-DRB1","HLA-DQB1","HLA-DQA1")
set3 <- c("CCL2","CCL3","CCL4","CCL5","CCL7","CCL13","CCL14","CCL18","CCL20","CCL23")
set4 <- c("CXCR2","CXCR3","CXCR4","CXCR6","CXCL1","CXCL2","CXCL6","CXCL8","CXCL13","CXCL16")
set5 <- c("IL2RB","IL6ST","IL13","IL16","IL21","IL22","IL26","XCL1","XCL2")
set6 <- c("TNFRSF4","TNFRSF6B","TNFRSF9","TNFRSF12A","TNFRSF14","FASLG","HGF","IFNG","CD27","VEGFA")

violinPlot(object = so, 
           assay = "SCT", 
           slot = "scale.data", 
           group.by = "cytospace_celltype", 
           group.subset = c(), 
           genes.of.interest)

source("GSEA.R")
source("GSEA_plot_path.R")