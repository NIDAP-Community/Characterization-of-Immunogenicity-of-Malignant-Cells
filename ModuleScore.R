library(Seurat)
library(tidyverse)
library(gridExtra)
library(quantmod)
library(grid)
library(data.table)
suppressMessages(library(dplyr))

SO_length <- length(SO@meta.data)

samples = eval(parse(text=gsub('\\[\\]','c()','c("X_ICC_18_Adjacent","X_ICC_18_Tumor","X_ICC_20_Tumor","X_ICC_23_Adjacent","X_ICC_23_Tumor","X_ICC_24_Tumor1","X_ICC_24_Tumor2","X_ICC_25_Adjacent")')))
sample_to_display <- c("X_ICC_18_Adjacent","X_ICC_18_Tumor","X_ICC_20_Tumor","X_ICC_23_Adjacent","X_ICC_23_Tumor","X_ICC_24_Tumor1","X_ICC_24_Tumor2","X_ICC_25_Adjacent")

if (length(samples) == 0) {
  samples = unique(SO@meta.data$sample_name)
}

colnames(SO@meta.data) <- gsub("orig_ident","orig.ident",colnames(SO@meta.data))

if("active.ident" %in% slotNames(SO)){
  sample_name = as.factor(SO@meta.data$orig.ident)
  names(sample_name)=names(SO@active.ident)
  SO@active.ident <- as.factor(vector())
  SO@active.ident <- sample_name
  SO.sub = subset(SO, ident = samples)
} else {
  sample_name = as.factor(SO@meta.data$orig.ident)
  names(sample_name)=names(SO@active.ident)
  SO@active.ident <- as.factor(vector())
  SO@active.ident <- sample_name
  SO.sub = subset(SO, ident = samples)
} 

rm(SO)
SO.sub_clean <- SO.sub

Markers <- unlist(Gene_list)

if (FALSE){
  protein_markers <- Markers[grepl("_prot",Markers)]
  
  protein_orig_markers <- gsub("_prot.*","",protein_markers)
  
  protein_markers_name <- paste(protein_orig_markers,
                                "_prot", sep = "")
  
  i = 0
  protein_array <- list()
  for (i in seq_along(protein_orig_markers)){
    protein_array[[i]] <- SO.sub@assays$Protein[protein_orig_markers[i],]
    rownames(protein_array[[i]]) <- protein_markers_name[i]
  }
  protein_array_comp <- do.call(rbind,protein_array)
  SO.sub@assays$SCT@data <- rbind(SO.sub@assays$SCT@data,protein_array_comp)
}

### Add negative/low identifiers to Module Scores
neg_markers_names <- Markers[grepl("_neg",Markers)]
orig_markers <- gsub("_neg.*","",neg_markers_names)

# Append neg_markers_names to rownames of SO.sub

neg_markers_list <- list() # Create a list for storage and retrieval

# Calculate adjusted expression for negative markers
for (i in seq_along(orig_markers)){
  
  # Format the data so that it can rbinded with SO$SCT@scale.data
  neg_markers_list[[i]] <- t(matrix(max(SO.sub@assays$SCT@data[orig_markers[i],]) - SO.sub@assays$SCT@data[orig_markers[i],]))
  row.names(neg_markers_list[[i]]) <- neg_markers_names[i]
  colnames(neg_markers_list[[i]]) <- colnames(SO.sub@assays$SCT@data)
  
  # Append new Negative/low marker (w Expression Count) to SO slot
  SO.sub@assays$SCT@data <- rbind(SO.sub@assays$SCT@data, neg_markers_list[[i]]) 
}

select_vector <- c("Malignant_Cells","Cholangiocytes","Hepatocytes","B_Cells","T_Cells","NK_Cells","Macrophages","Dendritic_Cells","Fibroblasts","Endothelial_Cells","CD8","CD4")
marker = select(Gene_list, select_vector)
marker.list = as.list(marker)

thres_vec <- c(0.13,0.65,0.07,0.1,0.25,0.12,0.15,0.17,0.08,0.12,0.55,0.5)
if (length(thres_vec) != length(select_vector)){
  if (sum(thres_vec) == 0){
    thres_vec <- rep(0, length(select_vector))
    print("Manual threshold set to zero - outputing preliminary data")
  } else {
    stop("Manual threshold length does not match number of celltypes to analyze - please check manual thresholds")
  }}

figures <- list()
exclude_cells <- c()

i = 0
j = 1

for (i in seq_along(marker.list)) {
  print(names(marker.list[i]))
  present=lapply(marker.list[[i]], function(x) x %in% rownames(SO.sub)) # apply function(x) x %in% rownames(SO.sub) to each element of marker.list
  absentgenes = unlist(marker.list[[i]])[present==FALSE]
  presentgenes = unlist(marker.list[[i]])[present==TRUE]
  print(paste0("Genes not present: ",paste0(absentgenes,collapse=",")))
  print(paste0("Genes present: ",paste0(presentgenes,collapse=",")))
  
  if(length(presentgenes) == 0){
    print(paste0(names(marker.list[i]), " genes were not found in SO and will not be analyzed"))
    exclude_cells[j] <- i
    j = j + 1
  }}  

if (length(exclude_cells) > 0){
  marker.list <- marker.list[-exclude_cells]} else {
    marker.list <- marker.list
  }  

for (i in seq_along(marker.list)) { 
  SO.sub=AddModuleScore(SO.sub,marker.list[i],name = names(marker.list[i]))
  
  m = 0 # m will be plugged into for Seurat Object
  
  m = paste0(names(marker.list[i]),"1")
  SO.sub@meta.data[[m]] <- scales::rescale(SO.sub@meta.data[[m]], to=c(0,1))
  
  clusid = SO.sub@meta.data[[m]] 
  
  d <- density(clusid) # create a density plot for ModScore vs Number of cells
  
  #hist(clusid[!is.na(clusid)], breaks=100, main=m)
  #abline(v=midpt,col="red",lwd=2)
  reduction = "tsne"
  if(reduction=="tsne"){
    p1 <- DimPlot(SO.sub, reduction = "tsne", group.by = "ident") 
  } else if(reduction=="umap"){
    p1 <- DimPlot(SO.sub, reduction = "umap", group.by = "ident")
  } else { 
    p1 <- DimPlot(SO.sub, reduction = "pca", group.by = "ident")
  }
  
  if(reduction=="tsne"){
    clusmat=data.frame(ident=p1$data$ident,umap1=p1$data$tSNE_1,umap2=p1$data$tSNE_2, clusid=as.numeric(SO.sub@meta.data[[m]]))
  } else if(reduction=="umap"){
    clusmat=data.frame(ident=p1$data$ident,umap1=p1$data$UMAP_1,umap2=p1$data$UMAP_2, clusid=as.numeric(SO.sub@meta.data[[m]]))
  } else { 
    clusmat=data.frame(ident=p1$data$ident,umap1=p1$data$PC_1,umap2=p1$data$PC_2, clusid=as.numeric(SO.sub@meta.data[[m]]))
  }
  clusmat <- mutate(clusmat, sample_clusid = clusmat$clusid * grepl(paste(sample_to_display, collapse = "|"), clusmat$ident))
  
  clusmat %>% group_by(clusid) %>% summarise(umap1.mean=mean(umap1), umap2.mean=mean(umap2)) -> umap.pos
  title=as.character(m)
  clusmat %>% dplyr::arrange(clusid) -> clusmat
  
  # Dimension reduction
  g <- ggplot(clusmat, aes(x=umap1, y=umap2)) +
    theme_bw() +
    theme(legend.title=element_blank()) +
    geom_point(aes(colour=sample_clusid),alpha=0.5,shape = 20,size=1) +
    #scale_color_gradient2(low = "blue4", mid = "white", high = "red",
    #          midpoint = midpt[[p]], na.value="grey",limits = c(0, 1)) + 
    scale_color_gradientn(colours = c("blue4","lightgrey", "red"), values = scales::rescale(c(0,thres_vec[i]/2,thres_vec[i],(thres_vec[i]+1)/2,1), limits = c(0, 1))) + guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    ggtitle(paste0("Number of present genes: ", length(presentgenes))) +
    xlab("tsne-1") + ylab("tsne-2")
  
  clusid.df <- data.frame(id=SO.sub@meta.data$orig.ident,score=SO.sub@meta.data[[m]])
  g1 = RidgePlot(SO.sub,features=m,group.by="orig.ident") + theme(legend.position = "none", title = element_blank(), axis.text.x = element_text(size = 6)) + geom_vline(xintercept = thres_vec[i], linetype = "dashed", color = "red3") + scale_x_continuous(breaks = seq(0,1,0.1))
  
  # Violin Plot
  clusid.df <- data.frame(id=SO.sub@meta.data$orig.ident,ModuleScore=SO.sub@meta.data[[m]])
  g2 = ggplot(clusid.df,aes(x=id,y=ModuleScore)) + geom_violin(aes(fill=id)) +  theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.title = element_blank(), panel.background = element_blank(), axis.text.x=element_blank(),legend.text=element_text(size=rel(0.8)),legend.position="top", axis.text.y = element_text(size = 6)) + guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) + geom_hline(yintercept = thres_vec[i], linetype = "dashed", color = "red3") + scale_y_continuous(breaks = seq(0,1,0.1))
  
  # Color gradient density plot
  g3 = ggplot(data.frame(x = d$x, y = d$y), aes(x, y)) + xlab("ModuleScore") + ylab("Density") + geom_line() + 
    geom_segment(aes(xend = d$x, yend = 0, colour = x)) + scale_y_log10() +
    scale_color_gradientn(colours = c("blue4","lightgrey", "red"), values = scales::rescale(c(0,thres_vec[i]/2,thres_vec[i],(thres_vec[i]+1)/2,1), limits = c(0, 1))) + geom_vline(xintercept = thres_vec[i], linetype = "dashed", color = "red3") + geom_vline(xintercept = thres_vec[i], linetype = "dashed", color = "red3") + scale_x_continuous(breaks = seq(0,1,0.1)) + theme(legend.title = element_blank(), axis.text.x = element_text(size = 6))
  
  # Set title for grid.arranged final figure
  
  last.figure = grid.arrange(g,g1,g2,g3, ncol=2, top=textGrob(names(marker.list[i]), gp = gpar(fontsize = 14, fontface = "bold")))
  
  figures[[i]] <- last.figure}

### Housekeeping before classification

## Name threshold vectors for easier retrieval in downstream steps
names(thres_vec) <- select_vector

## Rename meta.data Module Score columns so that it no longer has "1" at its end
names_repl <- substr(colnames(SO.sub@meta.data[(SO_length +1): length(SO.sub@meta.data)]), 1, nchar(colnames(SO.sub@meta.data[(SO_length +1): length(SO.sub@meta.data)]))-1) # nchar returns a vector containing the size of the character vector that was used as an argument. Subtracting 1 from the total size effectively allows you to omit the last character
names_orig <- colnames(SO.sub@meta.data[(SO_length +1): length(SO.sub@meta.data)]) # names to be replaced
setnames(SO.sub@meta.data, names_orig, names_repl) # remove the "1" from the end of the Module Score columns

### New Approach:
General_class <- c("Malignant_Cells","Cholangiocytes","Hepatocytes","B_Cells","T_Cells","NK_Cells","Macrophages","Dendritic_Cells","Fibroblasts","Endothelial_Cells")

#subset the columns of the metadata containing module scores only
SO_Trunc_Metadata_General <- SO.sub@meta.data[General_class] # subsetted metadata 

General_thres_vec <- thres_vec[General_class]

## See if elements in each ModScore column exceeds CellType threshold, set elements below threshold to zero. Keep values of elements above threshold
storage_list_MS_calls <- list()

Predict_Cell_from_ModScore <- function(ModScore_Metadata,thres_vec,rejection){
  
  thres_ls <- list()
  for (i in 1:ncol(ModScore_Metadata)){
    thres_ls[[i]]<- rep(thres_vec[i],nrow(ModScore_Metadata))
  }
  thres_df <- data.frame(matrix(unlist(thres_ls),nrow = nrow(ModScore_Metadata)))
  
  thres_filter <- ModScore_Metadata > thres_df
  ModScore_Metadata_post_thres_filter <- ModScore_Metadata * thres_filter
  
  ## Find column number with highest modscore
  max_col_vector <- max.col(ModScore_Metadata_post_thres_filter)
  
  # If a row contains all zeroes, they will be labeled with unknown
  all_zero_filter <- as.integer(!apply(ModScore_Metadata_post_thres_filter, 1, function(find_zero_rows) all(find_zero_rows == 0)))
  
  # Final filtering: 
  final_filter <- (max_col_vector * all_zero_filter) + 1
  
  # Original names appended to "unknown" classification for cells with ModScores below threshold
  appended_names <- c(rejection, names(ModScore_Metadata))
  
  # Added the names into a Likely_CellType Column
  dupl_data <- ModScore_Metadata
  dupl_data[,"Likely_CellType"] <- appended_names[final_filter]
  return(dupl_data)
}

General_output <- Predict_Cell_from_ModScore(SO_Trunc_Metadata_General,General_thres_vec,rejection = "unknown")
table(General_output$Likely_CellType)

if (TRUE){
  ## Subclass Identification
  Sub_class_storage <- c("T_Cells-CD8","T_Cells-CD4")
  Sub_class_calls <- list()
  
  parent_class <- unique(gsub("(.*)-(.*)","\\1",Sub_class_storage))
  
  for (parent in parent_class){
    Sub_class <- Sub_class_storage[grepl(parent,Sub_class_storage)]
    children_class <- gsub("(.*)-(.*)","\\2",Sub_class)
    
    # Subset out cells predicted to be CD8. CD8_Alloreactive cells will be screened from this population
    parents <- rownames(General_output[General_output$Likely_CellType == parent,])
    
    SO_Trunc_Metadata_parents <- SO.sub@meta.data[parents,] %>% select(children_class)
    
    for (children in children_class){
      
      plot_title <- paste("Density plot for",children,"Module Scores within",           parent,"population", sep = " ")
      adjusted_density_plot <- ggplot(SO_Trunc_Metadata_parents, aes_string(x =                      children)) + geom_density() + ggtitle(plot_title) +                  geom_vline(xintercept = thres_vec[children], linetype = "dashed", color =              "red3")
      figures[[length(figures) + 1]] <- adjusted_density_plot
      
    }
    
    # Create general annotation with barcoded cell and corresponding identity. M1 and M2 identities will be appended to this table after second round of classification
    
    SO_Trunc_Metadata_no_parents <- General_output[!General_output$Likely_CellType == parent,]
    non_parents <- rownames(SO_Trunc_Metadata_no_parents)
    
    gen_annot_table <- data.frame(cells = non_parents, identity = SO_Trunc_Metadata_no_parents$Likely_CellType)
    
    # Repeat Module Score Comparison and Cell Prediction with Macrophage Subset:
    children_thres_vec <- thres_vec[children_class]
    
    Sub_class_calls[[match(parent,parent_class)]] <- Predict_Cell_from_ModScore(SO_Trunc_Metadata_parents,children_thres_vec,rejection = parent) %>% select(Likely_CellType)
  }
  
  ## Updating CellType(s) in metadata with subclass calls
  final_subclass_results <- do.call(rbind,Sub_class_calls)
  
  parent_class_exc_vec <- paste(unique(parent_class),collapse = "|")
  General_output_final <- General_output[!grepl(parent_class_exc_vec,General_output$Likely_CellType),] %>% select(Likely_CellType)
  
  appendable_results <- rbind(General_output_final,final_subclass_results)
  subset_vector <- colnames(SO.sub@meta.data)[!colnames(SO.sub@meta.data) %in% colnames(SO.sub_clean@meta.data)]
  SO.sub_clean@meta.data <- cbind(SO.sub_clean@meta.data,SO.sub@meta.data[subset_vector])
  
  SO.sub_clean@meta.data$Likely_CellType <- appendable_results$Likely_CellType[match(rownames(SO.sub_clean@meta.data),rownames(appendable_results))]} else {
  subset_vector <- colnames(SO.sub@meta.data)[!colnames(SO.sub@meta.data) %in% colnames(SO.sub_clean@meta.data)]
  SO.sub_clean@meta.data <- cbind(SO.sub_clean@meta.data,SO.sub@meta.data[subset_vector])
  SO.sub_clean@meta.data$Likely_CellType <- General_output$Likely_CellType[match(rownames(SO.sub_clean@meta.data),rownames(General_output))]
}

## Set Image Parameters for Module Score figures       
n = ceiling(length(marker.list)^0.5)
m = ceiling(length(marker.list)/n)
imageWidth = 3200*m
imageHeight = 1600*n
dpi = 300

if (imageType == 'png') {
  png(
    filename=graphicsFile,
    width=imageWidth,
    height=imageHeight,
    units="px",
    pointsize=4,
    bg="white",
    res=dpi,
    type="cairo")
} else {
  library(svglite)
  svglite::svglite(
    file=graphicsFile,
    width=round(imageWidth/dpi,digits=2),
    height=round(imageHeight/dpi,digits=2),
    pointsize=1,
    bg="white")
}

Arranged_figures <- do.call("grid.arrange", c(figures))

print(Arranged_figures)

return(RFoundryObject(SO.sub_clean))