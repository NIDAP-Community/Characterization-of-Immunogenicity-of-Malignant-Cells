suppressMessages(library(dplyr))
suppressMessages(library(colorspace))
suppressMessages(library(dendsort))
suppressMessages(library(pheatmap))
suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(colorspace))

n <- 2e3
seed=6
set.seed(seed)
ourColorSpace <- colorspace::RGB(runif(n), runif(n), runif(n))
ourColorSpace <- as(ourColorSpace, "LAB")


distinctColorPalette <-function(k=1,seed) {
  currentColorSpace <- ourColorSpace@coords
  # Set iter.max to 20 to avoid convergence warnings.
  set.seed(seed)
  km <- kmeans(currentColorSpace, k, iter.max=20)
  colors <- unname(hex(LAB(km$centers)))
  return(colors)
}   

pal = function (n, h=c(237, 43), c=100, l=c(70, 90), power=1, fixup=TRUE, gamma=NULL, alpha=1, ...) {
  if (n < 1L) 
    return(character(0L))
  h <- rep(h, length.out = 2L)
  c <- c[1L]
  l <- rep(l, length.out = 2L)
  power <- rep(power, length.out = 2L)
  rval <- seq(1, -1, length = n)
  rval <- hex(
    polarLUV(
      L = l[2L] - diff(l) * abs(rval)^power[2L], 
      C = c * abs(rval)^power[1L],
      H = ifelse(rval > 0, h[1L], h[2L])
    ),
    fixup=fixup, ...
  )
  if (!missing(alpha)) {
    alpha <- pmax(pmin(alpha, 1), 0)
    alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
                    width = 2L, upper.case = TRUE)
    rval <- paste(rval, alpha, sep = "")
  }
  return(rval)
}
#Color selections for heatmap:
np0 = pal(100)
np1 = diverge_hcl(100, c=100, l=c(30, 80), power=1)  #Blue to Red
np2 = heat_hcl(100, c=c(80, 30), l=c(30, 90), power=c(1/5, 2))  #Red to Vanilla
np3 = rev(heat_hcl(100, h=c(0, -100), c=c(40, 80), l=c(75, 40), power=1)) #Violet to Pink
np4 = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(100))
np5 = colorRampPalette(c("steelblue","white", "red"))(100) #Steelblue to White to Red

np = list(np0, np1, np2, np3, np4, np5)
names(np) = c("Default","Blue to Red","Red to Vanilla","Violet to Pink","Bu Yl Rd","Bu Wt Rd")


doheatmap <- function(dat, clus, clus2, rn, cn, col) {
  require(pheatmap)
  require(dendsort)
  if (TRUE) {
    tmean.scale = t(scale(t(dat)))
    tmean.scale = tmean.scale[is.finite(rowSums(tmean.scale)),]
  } else {
    tmean.scale = dat
  }
  if(TRUE){
    quantperc <- 0.01
    upperquant <- 1-quantperc
    for(i in 1:nrow(tmean.scale)){
      data <- tmean.scale[i,] 
      dat.quant <- quantile(data,probs=c(quantperc,upperquant))
      data[data > dat.quant[2]] <- dat.quant[2]
      data[data < dat.quant[1]] <- dat.quant[1]
      tmean.scale[i,] <- data    
    }
  }
  col.pal <- np[[col]]
  if (FALSE) {
    col.pal = rev(col.pal)
  }
  # define metrics for clustering
  drows1 <- "euclidean"
  dcols1 <- "euclidean"
  minx = min(tmean.scale)
  maxx = max(tmean.scale)
  
  if (TRUE) {
    breaks = seq(minx, maxx, length=100)
    legbreaks = seq(minx, maxx, length=5)
  } else {
    #absmax = ceiling(max(abs(c(minx, maxx))))
    #breaks = c(-1*absmax, seq(0, 1, length=98), absmax)
    #legbreaks = c(-1*absmax, 0, absmax)
    breaks = seq(-2, 2, length=100)
    legbreaks = seq(-2, 2, length=5)
  }
  breaks = sapply(breaks, signif, 4)
  legbreaks = sapply(legbreaks, signif, 4)
  
  #Run cluster method using 
  hc = hclust(dist(t(tmean.scale)), method="complete")
  hcrow = hclust(dist(tmean.scale), method="complete")
  
  if (clus) {
    sort_hclust <- function(...) as.hclust(rev(dendsort(as.dendrogram(...))))
  } else {
    sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
  }
  if (clus2) {
    rowclus <- sort_hclust(hcrow)
  } else {
    rowclus = FALSE
  }
  #print('sorted the clusters')
  
  if (TRUE) {
    treeheight <- 25
  } else {
    treeheight <- 0
  }
  pathname <- stringr::str_replace_all("Differentially Expressed Genes (High Stem vs Low Stem)", "_", " ") 
  #pathname <- stringr::str_replace_all("T_cell_differentiation_GO_0030217", "_", " ") 
  #pathname <- stringr::str_replace_all("Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell", "_", " ")
  pathname <- stringr::str_wrap(pathname,50)
  
  hm.parameters <- list(
    tmean.scale, 
    color=col.pal,
    legend_breaks=legbreaks,
    # cellwidth=15, 
    # cellheight=10, 
    scale="none",
    treeheight_col=treeheight,
    treeheight_row=treeheight,
    kmeans_k=NA,
    breaks=breaks,
    # height=80,
    fontsize_row=5,
    fontsize_col=10,
    show_rownames=rn, 
    show_colnames=cn,
    main=pathname,
    clustering_method="complete",
    cluster_rows=rowclus, 
    cluster_cols=clus,
    cutree_rows=1,
    clustering_distance_rows=drows1, 
    clustering_distance_cols=dcols1,
    annotation_col = annotation_col,
    annotation_colors = annot_col,
    labels_col = labels_col
  )
  
  mat = t(tmean.scale)
  # print('calculated mat')
  
  callback = function(hc, mat) {
    # print('inside the callback')
    dend=rev(dendsort(as.dendrogram(hc)))
    # print ('reversed the dendsorted hc')
    dend %>% dendextend::rotate(c(1:length(dend))) -> dend
    as.hclust(dend)
  }
  do.call("pheatmap", c(hm.parameters, list(clustering_callback=callback)))
  #pheatmap(hm.parameters)
}

samples_to_include = c("X_ICC_18_Tumor","X_ICC_20_Tumor","X_ICC_23_Tumor","X_ICC_24_Tumor1","X_ICC_24_Tumor2")
samples_to_include <- samples_to_include[samples_to_include != ""]
samples_to_include <- gsub("-","_",samples_to_include)


so <- SO_malign

genes = c("ALDH3A1","SPP1","KRT6","KRT17","COL4A1","COL6A2","PSCA","S100P","FXYD3","KRT19","SLPI","AKR1C2","AGR2","WFDC2","LY6D","TACSTD2","CCL5","NKG7","CXCR4","MT2A","KLRB1","MIR205HG","GNLY","RGS1","CREM","SRGN")
genes = gsub(" ","",genes) 

if(genes[1] != ""){
  genesmiss = setdiff(genes,rownames(so$SCT@scale.data))
  if(length(genesmiss)>0){
    print(paste("missing genes:", genesmiss))
  }
}
genes = genes[genes %in% rownames(so$SCT@scale.data)]

if(!is.null(so@assays$Protein)){
  proteins = c("")
  proteins = gsub(" ","",proteins) 
  if(proteins[1] != ""){
    protmiss = setdiff(proteins,rownames(so$Protein@scale.data))
    if(length(protmiss)>0){
      print(paste("missing proteins:", protmiss))
    }
  }
  proteins = proteins[proteins %in% rownames(so$Protein@scale.data)]
}

df.mat1 = NULL
if(length(genes)>0){
  if(length(genes)==1){
    df.mat1 <- vector(mode="numeric",length=length(so$SCT@scale.data[genes,]))
    df.mat1 <- so$SCT@scale.data[genes,]
  }
  else{
    df.mat1 <- as.matrix(so$SCT@scale.data[genes,])
  }
}

df.mat2 = NULL
if(!is.null(so@assays$Protein)){
  if(length(proteins)>0){
    if(length(proteins)==1){
      df.mat2 <- vector(mode="numeric",length=length(so$Protein@scale.data[proteins,]))
      df.mat2 <- so$Protein@scale.data[proteins,]
      protname <- paste0(proteins,"_Prot")
    }
    else{
      df.mat2 <- as.matrix(so$Protein@scale.data[proteins,])
      protname <- paste0(proteins,"_Prot")
      rownames(df.mat2) <- protname
    }
  }
}

df.mat <- rbind(df.mat1,df.mat2)
if(!is.null(df.mat1)){
  rownames(df.mat)[rownames(df.mat)=="df.mat1"] <- genes
}
if(!is.null(df.mat2)){
  rownames(df.mat)[rownames(df.mat)=="df.mat2"] <- protname
}

df.mat <- df.mat[sort(rownames(df.mat)),]

if(FALSE){
  row.order = c("")  
  row.order = row.order[row.order %in% rownames(df.mat)]
  row.order = c(row.order,setdiff(rownames(df.mat),row.order))
  df.mat <- df.mat[row.order,]
  clusrows <- FALSE
}else{
  clusrows <- TRUE
}

annot <- SO_malign@meta.data
annot$cytotrace_call <- annot$cytotrace_celltype
head(annot)
annot %>% dplyr::filter(orig_ident %in% samples_to_include) -> annot
metadataplot <- c("cytotrace","cytotrace_call","orig_ident")
if(!"Barcode" %in% metadataplot){
  metadataplot = c(metadataplot,"Barcode")
}
annot %>% dplyr::select(metadataplot) -> annot
annot$cytotrace <- as.numeric(annot$cytotrace)
a = dim(annot)[2] - 1

if(FALSE){
  #prot = c("CD4","CD8")
  prot = c()
  if(length(prot) > 0){
    annot1 <- as.matrix(so$Protein@scale.data[prot,])
    if(length(prot)==1){
      annot1 <- annot1[match(annot$Barcode,rownames(annot1))]
      protname <- paste0(prot,"_Prot")
    }
    else{
      annot1 <- annot1[,match(annot$Barcode,colnames(annot1))]
      annot1 <- t(annot1)
      colnames(annot1) = paste0(colnames(annot1),"_Prot")
    }
  }
  #rna = c("CD8A","CD4")
  rna = c()
  if(length(rna)>0){
    annot2 <- as.matrix(so$SCT@scale.data[rna,])
    if(length(rna)==1){
      annot2 <- annot2[match(annot$Barcode,rownames(annot2))]
    }
    else{    
      annot2 <- annot2[,match(annot$Barcode,colnames(annot2))]
      annot2 <- t(annot2)
    }
  }
}


if(exists("annot1")){
  annot <- cbind(annot,annot1)
  colnames(annot)[colnames(annot)=="annot1"] <- protname
}
if(exists("annot2")){
  annot <- cbind(annot,annot2)
  colnames(annot)[colnames(annot)=="annot2"] <- rna
}
print(head(annot))

if(TRUE){
  annot %>% arrange_(.dots=metadataplot) -> annot
  df.mat <- df.mat[,match(annot$Barcode,colnames(df.mat))] 
  df.mat <- df.mat[ , apply(df.mat, 2, function(x) !any(is.na(x)))]
  cluscol = FALSE
}
else{
  cluscol = TRUE
}

#groups=gsub("'\"'","",paste0(groups,collapse=","))

annotation_col = as.data.frame(unclass(annot[,!(names(annot) %in% "Barcode")]))
annotation_col %>% mutate_if(is.logical, as.factor) -> annotation_col
rownames(annotation_col) <- annot$Barcode
if(dim(annot)[2] == 2){
  annottitle = colnames(annot)[1]
  colnames(annotation_col) = annottitle
}
annot_col = list()
groups=colnames(annotation_col)

q=sum(apply(annotation_col[,1:a],2,function(x){length(unique(x))})) 
colors=distinctColorPalette(q,5)

#colors <- c("darkred","greenyellow","darkviolet","black","darkorange","darkorchid","darkturquoise","darkblue","azure","cadetblue","chocolate","deeppink","lavender")
b=1
i=1
nam = NULL
col <- NULL
annot_col <- NULL
for (i in 1:length(groups)){
  nam <- groups[i]
  if(class(annotation_col[,i]) != "numeric"){
    grp <- as.factor(annotation_col[,i])
    c <- b+length(levels(grp))-1
    col = colors[b:c]
    names(col) <- levels(grp)
    assign(nam,col)
    annot_col = append(annot_col,mget(nam))
    b = c+1
    i=i+1
  }
  else{
    grp <- annotation_col[,i]
    np5 = colorRampPalette(c("steelblue","white", "red"))(length(grp))
    col=np5
    #names(col) <- grp
    assign(nam,col)
    annot_col = append(annot_col,mget(nam))
  }
}

print(paste0("The total number of genes in heatmap: ", nrow(df.mat)))

labels_col <- colnames(df.mat)

manually_replace_sample_names = FALSE
if (manually_replace_sample_names) {
  replacements = c("")
  old <- c()
  new <- c()
  for (i in 1:length(replacements)) {
    old[i] <- strsplit(replacements[i], ": ?")[[1]][1]
    new[i] <- strsplit(replacements[i], ": ?")[[1]][2]
  }
  df.relabel <- as.data.frame(cbind(old, new), stringsAsFactors=FALSE)
  labels_col %>% replace(match(df.relabel$old, labels_col), df.relabel$new) -> labels_col
}


p = doheatmap(dat=df.mat, clus=cluscol, clus2=clusrows, rn=TRUE, cn=FALSE, col="Bu Yl Rd")


return(annot)