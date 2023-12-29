library(Seurat)
library(reshape2)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(stringr)
library(ggplot2)
library(tidyverse)

localFilePaths <- paste("/rstudio-files/ccbr-data/users/Jing/CCBR_1026/Characterization-of-Immunogenicity-of-Malignant-Cells/data/",
                        list.files("/rstudio-files/ccbr-data/users/Jing/CCBR_1026/Characterization-of-Immunogenicity-of-Malignant-Cells/data"), sep = "")

subsetRegex = eval(parse(text=gsub('\\[\\]','c()','c()')))
Keep <- TRUE
if (length(subsetRegex) > 0) {
  if (Keep == TRUE){
    for (i in length(subsetRegex)) {
      localFilePaths <- localFilePaths[grepl(subsetRegex[[i]],localFilePaths)]
    }
  } else{
    for (i in length(subsetRegex)) {
      localFilePaths <- localFilePaths[!grepl(subsetRegex[[i]],localFilePaths)]
    }
  }
}

#obj.list <- lapply(localFilePaths, function(x) { return(c(Read10X_h5(x, use.names=FALSE), Read10X_h5(x, use.names=TRUE))) })
obj.list <- lapply(localFilePaths, function(x) {return(Read10X_h5(x, use.names=TRUE)) })

names(obj.list) <- c("X_ICC_18_Adjacent","X_ICC_18_Tumor","X_ICC_20_Tumor","X_ICC_23_Adjacent","X_ICC_23_Tumor","X_ICC_24_Tumor1","X_ICC_24_Tumor2","X_ICC_25_Adjacent")
obj.list <- obj.list[sort(names(obj.list))]

mincells = 3
mingenes = 501
organism = "Human"

mitoch = "^MT-"

seurat_object <- function(i) {
  if (class(obj.list[[i]]) == "dgCMatrix"){
    so.nf <- CreateSeuratObject(counts = obj.list[[i]], assay = "RNA", project=names(obj.list)[[i]], min.cells = 1, min.features = 0)
  } else {
    so.nf <- CreateSeuratObject(counts = obj.list[[i]][1]$`Gene Expression`, assay = "RNA", project=names(obj.list)[[i]], min.cells = 1, min.features = 0)
  }
  so.nf <- NormalizeData(so.nf, normalization.method = "LogNormalize", scale.factor = 10000)
  so.nf[["percent.mt"]] <- PercentageFeatureSet(object = so.nf, pattern = mitoch)
  so.nf$log10GenesPerUMI <- log10(so.nf$nFeature_RNA) / log10(so.nf$nCount_RNA)
  
  #Filtered Seurat Object:
  if (class(obj.list[[i]]) == "dgCMatrix"){
    so <- CreateSeuratObject(counts = obj.list[[i]], assay = "RNA", project=names(obj.list)[[i]], min.cells = mincells, min.features = mingenes)
  } else {
    so <- CreateSeuratObject(counts = obj.list[[i]][1]$`Gene Expression`, assay = "RNA", project=names(obj.list)[[i]], min.cells = mincells, min.features = mingenes)
  }
  so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
  so[["percent.mt"]] <- PercentageFeatureSet(object = so, pattern = mitoch)
  # so <- CellCycleScoring(object = so, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes)
  so$log10GenesPerUMI <- log10(so$nFeature_RNA) / log10(so$nCount_RNA)
  
  so.origcount = dim(so)[2]

  #Start with filtering here:
  so <- subset(so, subset = nCount_RNA > 2001)
  
  maxgenes = 6000
  complexity = 0.25
  MAD_gene <- FALSE
  ngenestdev <- mad(so@meta.data$nFeature_RNA)
  ngenemed <- median(so@meta.data$nFeature_RNA)
  ngenemaxlim <- ngenemed+(3*ngenestdev)
  gl = format(round(ngenemaxlim,0),nsmall=0)
  
  maxmitoch = 20
  
  MAD_mitoch <- FALSE
  mitostdev <- mad(so@meta.data$percent.mt)
  mitomed <- median(so@meta.data$percent.mt)
  mitomaxlim <- mitomed+(3*mitostdev)
  ml = format(round(mitomaxlim,2),nsmall=2)
  
  if (MAD_gene == TRUE & MAD_mitoch == TRUE)       {
    cat(paste0("Gene Count Filter = low:",mingenes," high:",gl),"\n")
    cat(paste0("Mitochondrial Percentage Filter =",ml,"\n"))
    cat(paste0("Complexity Filter =",complexity,"\n"))
    so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < ngenemaxlim & so@meta.data$percent.mt < mitomaxlim & so@meta.data$log10GenesPerUMI > complexity), ]))
    perc.remain = (dim(so)[2]/so.origcount)*100
    perc.remain=formatC(perc.remain,format = "g",digits=3)
    cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
    cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
  } else if (MAD_gene == FALSE & MAD_mitoch == TRUE) {
    cat(paste0("Gene Count Filter = low:", mingenes," high:", maxgenes),"\n")
    cat(paste0("Mitochondrial Percentage Filter =",ml,"\n"))
    so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < maxgenes & so@meta.data$percent.mt < mitomaxlim & so@meta.data$log10GenesPerUMI > complexity), ]))
    perc.remain = (dim(so)[2]/so.origcount)*100
    perc.remain=formatC(perc.remain,format = "g",digits=3)
    cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
    cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
  } else if (MAD_gene == TRUE & MAD_mitoch == FALSE){
    cat(paste0("Gene Count Filter = low:",mingenes," high:",gl),"\n")
    cat(paste0("Mitochondrial Percentage Filter =", maxmitoch,"\n"))
    so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < ngenemaxlim & so@meta.data$nFeature_RNA > mingenes & so@meta.data$percent.mt < maxmitoch & so@meta.data$log10GenesPerUMI > complexity), ]))
    perc.remain = (dim(so)[2]/so.origcount)*100
    perc.remain=formatC(perc.remain,format = "g",digits=3)
    cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
    cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
  } else {
    cat(paste0("Gene Count Filter = low:", mingenes," high:", maxgenes),"\n")
    cat(paste0("Mitochondrial Percentage Filter =", maxmitoch,"\n"))
    so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < maxgenes & so@meta.data$nFeature_RNA > mingenes & so@meta.data$percent.mt < maxmitoch & so@meta.data$log10GenesPerUMI > complexity), ]))
    perc.remain = (dim(so)[2]/so.origcount)*100
    perc.remain=formatC(perc.remain,format = "g",digits=3)
    cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
    cat(paste0("Percent Remaining=" ,perc.remain),"\n\n")
  }
  
  df.m <- melt(so@meta.data)
  df.m$filt <- "filt"
  df.m$filt <- as.factor(df.m$filt)
  df2.m <- melt(so.nf@meta.data)
  df2.m$filt <- "raw"
  df2.m$filt <- as.factor(df2.m$filt)
  
  so2.list <- list(so,so.nf)
  
  return(so2.list)
}

so.list <- lapply(seq_along(obj.list), seurat_object)

so.f.list <- lapply(so.list,function(x) x[[1]])
names(so.f.list) <- sapply(names(obj.list), function(x) gsub("_filtered.h5", "", x))

so.final.list <<- so.f.list

SO = unnamed_2$value

#in case you want to redo this on a merged SO
if (class(SO) =="Seurat") {
  x =list()
  x[[1]] <- SO
  SO <- x
}
vars_to_regress <- c("percent.mt","nCount_RNA")
npcs = 30
linearScale = FALSE

# Linearly scale data without regressing anything.
scale_so <- function(so){
  so <- CellCycleScoring(object = so, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes)
  so$CC.Difference <- so$S.Score - so$G2M.Score
  so <- FindVariableFeatures(object = so, nfeatures = 2000, mean.cutoff = c(0.0125, 8), dispersion.cutoff = c(0.5, 100000), selection.method = "vst")
  all.genes <- rownames(so)
  so <- ScaleData(so,features=all.genes)
  return(so)
}

# Make PCA without regressing anything, and using only SCTransform().
pca_noregress <- function(so) {
  so <- SCTransform(so,do.correct.umi = FALSE,return.only.var.genes = FALSE)
  so <- RunPCA(object = so, features = VariableFeatures(object = so), npcs = npcs)
  return(so)
}

# Make PCA with SCTransform() and optional ScaleData, and do so with
# both regression (if user requests) and on all genes.
pca <- function(so) {
  # If user sets Linear Scaling toggle TRUE, also run ScaleData().
  # Use case: user has legacy project from Seurat 2 and wants to keep
  # methods consistent with pre-SCT Seurat.
  if(linearScale == TRUE) {
    all.genes <- rownames(so)
    if(is.null(vars_to_regress)){     
      so <- so
    }
    else{
      so <- ScaleData(so, features=all.genes, vars.to.regress = vars_to_regress) 
    }
  }
  # Run SCTransform().
  if(is.null(vars_to_regress)){
    so <- so
  }
  else { 
    so <- SCTransform(so,do.correct.umi = TRUE, vars.to.regress = vars_to_regress, return.only.var.genes = FALSE)
  }
  # Make PCA using last transform run, which will always be that from
  # SCTransform().
  so <- RunPCA(object = so, npcs = npcs)
  return(so)
}

# Do transformation with and without regression using SCTransform()
# and ScaleData().
so_scale <- lapply(SO, scale_so) 
SO <- lapply(so_scale, pca) 

#initialize Citeseq functionality as false, 
#later the template will check for a Protein assay and run if it finds it
doCiteSeq <- FALSE

doMergeData <- !FALSE
dat = vector()
integratedata = FALSE

if (length(SO) > 1) {
  for(i in 2:length(SO)){dat=c(dat,SO[[i]]) }
  SO_merge <- merge(SO[[1]], y = dat, add.cell.ids = names(SO), project = "scRNAProject", merge.data = TRUE)
  allgenes <- rownames(SO_merge)
  SO_merge <- ScaleData(SO_merge, assay = "RNA", features=allgenes)
} else {
  SO_merge <- SO[[1]]
  allgenes <- rownames(SO_merge)
  SO_merge <- ScaleData(SO_merge, assay = "RNA", features=allgenes)
}

if (!("orig.ident" %in% colnames(SO_merge@meta.data))) {
  SO_merge@meta.data$orig.ident <- SO_merge@meta.data$orig_ident
}

if ("Protein" %in% names(SO_merge@assays)){
  doCiteSeq <-TRUE
}

npcs = 10
Do_SCTransform = TRUE
vars_to_regress = c("percent.mt","nCount_RNA")

if (Do_SCTransform){
  if(is.null(vars_to_regress)){
    SO_merge <- SCTransform(SO_merge,do.correct.umi = TRUE, return.only.var.genes = FALSE)}
  else{       
    SO_merge <- SCTransform(SO_merge,do.correct.umi = TRUE,vars.to.regress=vars_to_regress,return.only.var.genes = FALSE) 
  }
} else {
  all.genes <- rownames(SO_merge)
  if(is.null(vars_to_regress)){
    SO_merge <- SO_merge
  }
  else{
    SO_merge <- ScaleData(SO_merge, features=all.genes, assay = "RNA", vars.to.regress=vars_to_regress) 
  }
  DefaultAssay(SO_merge) <- "RNA"   
}

if (length(SO)>1) {
  all_features <- lapply(SO, row.names) %>% Reduce(intersect, .)
  if(integratedata==TRUE){
    integ_features <- SelectIntegrationFeatures(object.list = SO, nfeatures = 3000) 
    if(!is.null(SO[[1]]@assays$SCT)){
      SO <- PrepSCTIntegration(object.list = SO, anchor.features = integ_features)
      k.filter <- min(200, min(sapply(SO, ncol)))
      integ_anchors <- FindIntegrationAnchors(object.list = SO, normalization.method = "SCT", k.filter=k.filter, anchor.features = integ_features)
      SO_merge <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT",features.to.integrate = all_features)
      SO_merge <- ScaleData(SO_merge,features=all_features)
    }
    else{
      k.filter <- min(200, min(sapply(SO, ncol)))
      integ_anchors <- FindIntegrationAnchors(object.list = SO, k.filter=k.filter, anchor.features = integ_features)
      SO_merge <- IntegrateData(anchorset = integ_anchors,features.to.integrate = all_features)
      SO_merge <- ScaleData(SO_merge,features=all_features)  
    }}
}


SO_merge <- FindVariableFeatures(object = SO_merge, nfeatures = 2000, mean.cutoff = c(0.0125, 8), dispersion.cutoff = c(0.5, 100000), selection.method = "vst", verbose = FALSE)
SO_merge <- RunPCA(object = SO_merge, npcs = npcs, verbose = FALSE,seed.use = 42)
SO_merge <- RunUMAP(object = SO_merge, reduction = "pca", dims = 1:npcs, seed.use=42)
SO_merge <- RunTSNE(object = SO_merge, reduction = "pca", dim.embed = 2, dims = 1:npcs, seed.use = 1)
SO_merge <- FindNeighbors(SO_merge, dims = 1:npcs)


#check for CITE-seq data and if so, run reductions
if(doCiteSeq) {
  # Moved below integration step. SO_merge is recreated and this information was lost
  SO_merge <- ScaleData(SO_merge, assay = "Protein")
  
  print("finding protein variable features...")
  VariableFeatures(SO_merge,assay="Protein") <- rownames(SO_merge$Protein)
  #Support for partial
  if(all(sapply(seq_along(SO),function(i) "Protein" %in% names(SO[[i]]@assays)))){
    print("running protein pca...")
    SO_merge <- RunPCA(object = SO_merge, assay="Protein",npcs = npcs,verbose = FALSE,reduction.name="protein_pca",seed.use = 42)
    SO_merge <- RunUMAP(object = SO_merge, assay="Protein", features=rownames(SO_merge$Protein), reduction.name="protein_umap",seed.use=42)
    SO_merge <- RunTSNE(object = SO_merge, assay="Protein", features=rownames(SO_merge$Protein),seed.use = 1,reduction.name="protein_tsne",check_duplicates=F)
    SO_merge <- FindNeighbors(SO_merge, assay="Protein",graph.name="Protein_snn",features=rownames(SO_merge$Protein))
  }else{
    doCiteSeq <- FALSE #set to false so we don't cluster protein
  }
  
} else {
  doCiteSeq <- FALSE
}

for (i in seq(0.3,1.2,0.3)){
  SO_merge <- FindClusters(SO_merge, resolution = i, algorithm = 1)
  if(doCiteSeq){
    SO_merge <- FindClusters(SO_merge, graph.name="Protein_snn",resolution = i, algorithm = 1)
  }
}
print("Clustering successful!")

return(SO_merge)
