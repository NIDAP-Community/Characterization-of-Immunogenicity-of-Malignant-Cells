
library(magrittr)
library(fgsea)
library(dplyr)
library(gtable)
library(grid)
library(gridExtra)
library(data.table)
library(patchwork)    
# INPUT PARAMS ====
## run gsea function

run.gsea <- function(dx, mode, collections, db, collapse, collapse_filter, collapse_threshold) {
  
  # compute gsea stats
  
  ranked = dx$genescores        
  names(ranked) = dx$gene_id
  db$inPathway = sapply(db$gene_symbol, function(x) paste(sort(x[x %in% names(ranked)]), collapse=","))
  
  if (mode == "over all collections") {
    
    set.seed(246642)
    gsea <- fgsea(pathways=collections, stats=ranked, minSize=15, maxSize=500, nperm=5000) 
    gsea$leadingEdge <- sapply(gsea$leadingEdge, function(x) paste(x, collapse=","))
    gsea <- dplyr::inner_join(gsea, select(db, collection, gene_set_name, inPathway), by=c("pathway"="gene_set_name")) %>% dplyr::select(collection, dplyr::everything())         
    
  } else {
    
    included_collections <- setNames(unique(db$collection), unique(db$collection))
    gsea <- lapply(included_collections, function(x) {
      set.seed(246642)
      gsea_collection = fgsea(pathways = collections[names(collections) %in% dplyr::filter(db, collection==x)$gene_set_name], stats=ranked, minSize=15, maxSize=500, nperm=5000)
      gsea_collection$leadingEdge <- sapply(gsea_collection$leadingEdge, function(x) paste(x, collapse=","))
      return(dplyr::inner_join(gsea_collection, select(db, collection, gene_set_name, inPathway) %>% filter(collection == x), by=c("pathway"="gene_set_name")) %>% dplyr::select(collection, dplyr::everything()) )                
    }) %>% dplyr::bind_rows()
    
  }
  
  # collapse function    
  run.collapse <- function(cp.input, pvalue, collections) {
    collapsedPathways <- collapsePathways(as.data.table(cp.input), collections, ranked, pval.threshold=pvalue) # requires the data.table library
    collapsedResults <- data.frame(pathway=names(collapsedPathways$parentPathway), parentPathway=collapsedPathways$parentPathway) 
    collapsedResults$mainPathway <- collapsedPathways$mainPathway[match(collapsedResults$pathway, collapsedPathways$mainPathways)]
    return(collapsedResults)            
  }
  
  if (collapse == "over all collections") {
    
    filter_gsea <- gsea %>% dplyr::select(-collection) %>% dplyr::distinct() %>% dplyr::filter(get(collapse_filter) <= collapse_threshold) %>% dplyr::arrange(pval)
    set.seed(246642)
    collapsedResults = run.collapse(cp.input = filter_gsea, pvalue = 0.05, collections = collections)
    gsea <- gsea %>% dplyr::left_join(collapsedResults, by=c('pathway'='pathway')) %>% dplyr::select(collection, pathway,  parentPathway, mainPathway, dplyr::everything())
    
  } else if (collapse == "within each collection") {
    
    filter_gsea = gsea %>% dplyr::filter(get(collapse_filter) <= collapse_threshold) %>% dplyr::arrange(pval) %>% dplyr::group_by(collection)
    set.seed(246642)
    collapsedResults = dplyr::group_modify(filter_gsea , ~run.collapse(., pvalue = 0.05, collections = collections)) %>% dplyr::ungroup()
    gsea <- gsea %>% dplyr::left_join(collapsedResults, by=c('pathway'='pathway','collection'='collection')) %>% dplyr::select(collection, pathway,  parentPathway, mainPathway, dplyr::everything())
    
  }
  
  return(gsea)
  
}


## pvalue cutoffs table functions

table.pvalue <- function(gsea) {
  
  cuts <- c(-Inf, 1e-04, 0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 1)
  cutsLab <- paste("<",cuts[-1], sep="")
  p= cumsum(table(cut(gsea$pval, breaks = cuts, labels=cutsLab, include.lowest = FALSE, right=TRUE)))
  q= cumsum(table(cut(gsea$padj, breaks = cuts, labels=cutsLab,include.lowest = FALSE, right=TRUE)))
  tab <- data.frame(cutsLab, p, q)
  colnames(tab) <- c("alpha","p-value","*adjusted\np-value")
  rownames(tab) <- NULL
  return(tab)
}       

plot.table <- function(dtab, score) {
  
  title <- textGrob(paste0(unique(dtab$contrast), score),gp=gpar(fontsize=10))
  tab <- as.data.frame.matrix(dtab %>% dplyr::select(-contrast))
  table <- tableGrob(tab, theme = ttheme_default(core=list(fg_params=list(cex=0.9)), colhead = list(fg_params=list(cex = 0.9, parse=FALSE)), rowhead=list(fg_params=list(cex= 0.6))))
  table <- gtable_add_rows(table, heights = grobHeight(title) + unit(2,"line"), pos=0)
  table <- gtable_add_grob(table, list(title), t=c(1), l=c(1), r=ncol(table)) 
  wrap_elements( table)
}

#..dataset

geneset_db = msigDB_v6_2 # gene set collection is pinned (can be changed but assumes presence of columns named: "species", "collection", "gene_set_name", "gene_symbol")
genescore_df = read.csv("~/files_for_analysis/deg_tab_gsea.csv")

#..variables

## ranking
genescore = "avg_logFC_high_vs_low"
genescore_alternative = c()
geneid = "Gene"
contrasts = c()    

## GSEA and pathway collapse
collections_to_include = c("H: hallmark gene sets","CP:REACTOME: Reactome gene sets","CP:KEGG: KEGG gene sets", "custom_geneset") # object Pathway/Collection lookup
FDR_correction_mode = "over all collections"
collapse_mode = "over all collections"
filterBy = "p-value"
filterBy = switch(filterBy, "p-value"="pval", "adjusted p-value"="padj")
filterBy_threshold = 0.05

## output
sortBy = c("collection","pval")
sortDecreasing = FALSE

## image size/resolution
#image: png
png(filename=graphicsFile, width=2500, height=2500, units="px", pointsize=4, bg="white", res=300, type="cairo")         

# DATA HANDLING ====    

#..gene set collection    
geneset_db <- geneset_db %>% SparkR::filter(geneset_db$species=="Human") %>% SparkR::filter(SparkR::`%in%` (geneset_db[["collection"]], collections_to_include)) 
db_unique <- geneset_db %>% SparkR::select("gene_set_name") %>% SparkR::distinct()
db_selected <- geneset_db %>% SparkR::select("collection","gene_set_name") %>% SparkR::distinct()
db_isDuplicated <- SparkR::count(db_selected) > SparkR::count(db_unique)

if (db_isDuplicated) {
  
  db_selected <- SparkR::collect(db_selected)
  within_collection <-  db_selected %>% dplyr::group_by(collection) %>% dplyr::filter(duplicated(gene_set_name)) %>% dplyr::ungroup() %>% dplyr::count()
  
  if (within_collection == 0) {        
    stop("ERROR: duplicated gene set names found in the 'Gene set database' due to overlapping collections selected by the 'Collections to include' parameter")
    
  } else if (within_collection > 0) {        
    between_collection <- Reduce("intersect", split(db_selected$gene_set_name, db_selected$collection)) %>% length()
    
    if (between_collection == 0) {
      stop("ERROR: duplicated gene set names found in the 'Gene set database' due to not unique gene set names within a collection")
      
    } else if (between_collection > 0) {
      "ERROR: duplicated gene set names found in the 'Gene set database' due to overlapping collections selected by the 'Collections to include' parameter and not unique gene set names within a collection"
    } 
  } 
  
} else if (!db_isDuplicated) {    
  geneset_db <- SparkR::collect(geneset_db)     
}

geneset_db <- geneset_db %>% dplyr::group_by(collection, gene_set_name) %>% dplyr::summarize(gene_symbol = as.list(strsplit(paste0(unique(gene_symbol), collapse = " "), " "))) %>% dplyr::ungroup()
geneset_list = geneset_db$gene_symbol
names(geneset_list) = geneset_db$gene_set_name

#..ranking
if (!is.null(genescore_alternative)) { genescore <- genescore_alternative }
rank_columns = "rank"
rank_contrasts = unlist(strsplit(rank_columns, genescore))
if (!is.null(contrasts)) {
  index = match(contrasts, rank_contrasts)
  rank_columns = rank_columns[index]
  rank_contrasts = rank_contrasts[index]
}    
genescore_df <- genescore_df %>% dplyr::select(geneid, rank_columns) %>% tidyr::pivot_longer(!geneid, names_to="contrast", values_to="genescores", values_drop_na=TRUE) %>% dplyr::rename("gene_id"=geneid) %>% dplyr::mutate(contrast=sub(genescore, "", contrast))
genescore_grouped <- dplyr::group_by(genescore_df, contrast)  


# STATISTICAL CALCULATION ====

#..GSEA and pathway

gsea <- dplyr::group_modify( genescore_grouped, ~run.gsea(., db=geneset_db, collections=geneset_list, mode=FDR_correction_mode, collapse=collapse_mode, collapse_filter=filterBy, collapse_threshold=filterBy_threshold) ) 

# OUTPUT ====

#..tables (image)
tab <- dplyr::group_modify(gsea, ~table.pvalue(.x)) %>% dplyr::ungroup()
ltab <- split(tab, tab$contrast) %>% lapply(function(x) plot.table(x, score=genescore) ) %>% wrap_plots()
print(ltab + plot_annotation(title="Cumulative number of significant calls (GSEA)", subtitle = sprintf("*p value adjusted %s by the method of Benjamini and Hochberg (1995)", FDR_correction_mode), tag_levels = 'A', theme=theme(plot.title = element_text(size=20, face='bold', hjust=0.5, margin = margin(t = 0)), plot.subtitle =  element_text(size=11, face='italic', hjust=0.5,  margin = margin(t = 10, b=20)))) )

#..Logs

#cat("\nThe number of tested gene sets per each collection and contrast\n")
#N <- count(gsea,collection); print(N, n=nrow(N))
#if (collapse_mode != "none") {
#    cat("\nThe number of main pathways identified per each collection and contrast\n")
#    N = count(gsea[[2]] %>% dplyr::filter(!is.na(mainPathway)), collection); print(N, n=nrow(N)) 
#}
#cat(sprintf("\nCumulative number of significant calls\np-value adjusted for the false discovery rate %s by the method of Benjamini and Hochberg (1995)\n", FDR_correction_mode))
#tab %>% print(n=nrow(tab))

#..dataset

if (sortDecreasing) { sortBy = sapply( sortBy, function(x) sprintf("desc(%s)", x) ) }

gsea <- gsea %>% dplyr::arrange_(.dots = sortBy) %>% tibble::add_column(geneScore = genescore, .after="contrast") 

return(gsea) 