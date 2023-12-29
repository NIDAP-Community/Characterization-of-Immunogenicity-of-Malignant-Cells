l2p_path <- function(upregulated = T){

suppressMessages(library(l2p))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
#suppressMessages(library(stringr))
suppressMessages(library(RCurl))
suppressMessages(library(RColorBrewer))

deg_res <- read.csv("~/files_for_analysis/KruskalWallis_DEG.csv")

deg_res %>% dplyr::select("Gene", "Rank") -> genesmat
genesmat %>% dplyr::arrange(desc(`Rank`)) -> genesmat
genesmat %>% dplyr::filter(!is.na(`Rank`)) -> genesmat

sort_descending = upregulated
if (sort_descending) {
  genes_to_include = head(genesmat["Gene"], 500)
} else {
  genes_to_include = tail(genesmat["Gene"], 500)
}
genes_to_include <- as.vector(unique(unlist(genes_to_include)))
genes_universe = as.vector(unique(unlist(genesmat["Gene"])))
gene_set_sources_to_include = c("GO","REACTOME","KEGG")
categories_string <- paste(gene_set_sources_to_include, collapse=",")

organism = "Human"
if (organism != "Human") {
  organism_of_interest = "Human"
  orthology_table = SparkR::filter(Ortholog_Map_for_RNA_Seq,     
                                   Ortholog_Map_for_RNA_Seq$Organism==organism_of_interest)
  
  if ("from human"=="from human") {
    orthology_reference_column = "Human_gene_name"
    orthology_conversion_column = "Nonhuman_gene_name"
  } else {
    orthology_reference_column = "Nonhuman_gene_name"
    orthology_conversion_column = "Human_gene_name"
  }
  orthology_table %>% SparkR::withColumnRenamed(orthology_reference_column,           "orthology_reference") %>%
    SparkR::withColumnRenamed(orthology_conversion_column, "orthology_conversion") %>% SparkR::select("orthology_reference", "orthology_conversion") -> orthology_table
  orthology_table <- SparkR::collect(orthology_table)
  orthology_table %>% dplyr::filter(orthology_conversion %in% genes_to_include) %>% dplyr::select(orthology_reference) -> genes_to_include
  genes_to_include <- as.vector(genes_to_include$orthology_reference)
  genes_to_include <- unique(genes_to_include)
  orthology_table %>% dplyr::filter(orthology_conversion %in% genes_universe) %>% dplyr::select(orthology_reference) -> genes_universe
  genes_universe <- as.vector(genes_universe$orthology_reference)
  genes_universe <- unique(genes_universe)
  
}else{
  genes_to_include = genes_to_include
  genes_universe = genes_universe
}

use_built_in_gene_universe = FALSE
if (use_built_in_gene_universe) {
  x <- l2pwcats(genes_to_include, categories_string)
  print("Using built-in gene universe.")
} else {
  x <- l2puwcats(genes_to_include, genes_universe, categories_string)
  print("Using all genes in differential expression analysis as gene universe.")
}
x %>%
  dplyr::arrange(pval) %>% dplyr::mutate(hitsPerc=(pwhitcount/(pwnohitcount+pwhitcount))*100) %>% dplyr::mutate(pathtotal=pwnohitcount+pwhitcount) %>% dplyr::filter(ratio >= 0) %>%
  dplyr::select("pathwayname", "category", "pathwayaccessionidentifier", "pval", "fdr", "pwhitcount", "genesinpathway", "pwnohitcount","pathtotal","hitsPerc", "inputcount", "pwuniverseminuslist","ratio")  %>% dplyr::rename(diff_ratio = ratio) %>% dplyr::filter(pval < 0.05)%>% dplyr::filter(pwhitcount >= 5) -> x
#allgenes <- lapply(x$pathwayaccessionidentifier,function(x) l2pgetgenes4acc(x)) 
#x %>% mutate("allgenes" = paste(allgenes,sep="",collapse=" ")) -> x
print(paste0("Total number of pathways: ", nrow(x)))
goResults <- x
goResults %>% top_n(10, wt=-log(pval)) %>%
  dplyr::arrange(-log(pval)) -> goResults
minp = min(goResults$pval) - 0.1*min(goResults$pval)
maxp = max(goResults$pval) + 0.1*max(goResults$pval)
#print(goResults$pval)
sizemax = ceiling(max(goResults$pwhitcount)/10)*10  
goResults %>% dplyr::mutate(pathwayname2 = stringr::str_replace_all(pathwayname, "_", " ")) -> goResults
goResults %>% dplyr::mutate(pathwayname2 = stringr::str_wrap(pathwayname2,30)) -> goResults

if (FALSE){
  goResults %>% dplyr::mutate(percorder = order(goResults$pval)) -> goResults
  goResults$pathwayname2 <- factor(goResults$pathwayname2, levels = goResults$pathwayname2[goResults$percorder])
  xmin = floor(min(goResults$pval))
  xmax = max(goResults$pval) 
  gplot <- goResults %>% 
    ggplot(aes(x=pval,
               y=pathwayname2, 
               colour=hitsPerc, 
               size=pwhitcount)) +
    geom_point() +
    theme(text = element_text(size=20), legend.position = "right", legend.key.height = unit(1, "cm"),
          axis.title.x = element_text(size = rel(1.2)),axis.title.y = element_text(size = rel(1.2))) +
    xlim(xmin,xmax) +
    expand_limits(colour = seq(minp, maxp, by = 10),
                  size = seq(0, sizemax,by=10)) +
    labs(x="p value", y="GO term", colour="Hits (%)", size="Count") 
  print(gplot)
}else{
  goResults %>% dplyr::mutate(percorder = order(goResults$hitsPerc)) -> goResults
  goResults$pathwayname2 <- factor(goResults$pathwayname2, levels = goResults$pathwayname2[goResults$percorder])
  xmin = floor(min(goResults$hitsPerc)-5)
  xmax = ceiling(max(goResults$hitsPerc)+5) 
  gplot <- goResults %>% 
    ggplot(aes(x=hitsPerc,
               y=pathwayname2, 
               colour=pval, 
               size=pwhitcount)) +
    geom_point() +
    theme_classic() +
    theme(text = element_text(size=20), legend.position = "right", legend.key.height = unit(1, "cm"),
          axis.title.x = element_text(size = rel(1.2)),axis.title.y = element_text(size = rel(1.2))) +
    xlim(xmin,xmax) +
    expand_limits(colour = seq(minp, maxp, by = 10),
                  size = seq(0, sizemax,by=10)) +
    labs(x="Hits (%)", y="GO term", colour="p value", size="Count") 
  print(gplot)
}
return(x)
}