library(Seurat)
library(harmony)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

so = SO

seur.SCT <- so@assays$SCT@scale.data
seur.SCT[1:10, 1:10]
# try small
neigens <- 4000
ngenes <- 2000 
# length(VariableFeatures(so)) 2000, 3/25/21
genes_of_interest <- c("HLA-A","HLA-B","HLA-C","IL6ST","TGFBR2","STAT1","STAT3","CD44","NANOG","SOX2","OCT4","TM4SF4","ANXA4","CDH6","POU5F1")
genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(so@assays$SCT@scale.data)]
VariableFeatures(so) <- c(VariableFeatures(so), genes_of_interest) # can add more genes here (in case they are not present after VariableFeatures)
mvf <- VariableFeatures(so)[1:(ngenes + length(genes_of_interest))] # trusting mvf are in descending order (ngenes + length(genes_of_interest))
# code to add annotation genes to to mvf, from input genelist
# nmissing <- sum(!(genelist %in% mvf))
# mvf <- VariableFeatures(so)[1:(ngenes-nmissing)] # make room for extra genes
# mvf <- unique(c(mvf, genelist)) 
# seur.SCT.obj <- so@assays$SCT
print(dim(so@assays$SCT@scale.data))
seur.SCT <- seur.SCT[mvf,]
seur.SCT <- t(seur.SCT) # want gene columns and cell rows
# need the original loadings and embeddings to compare with my SVD
seur.loads <- so@reductions$pca@feature.loadings
seur.pca <- so@reductions$pca@cell.embeddings
sm.pca <- seur.pca[1:10, 1:10]
sm.seur.ldngs <- seur.loads[1:10, 1:10]
##################################################
# SVD on scaled counts
##################################################
ngenes <- length(mvf)
Sys.time()
pppca <- svd(seur.SCT) # previously propack, hence pppca.  Different argumnts, with svd using defaults
Sys.time() # 27 seconds for 900 genes, 5K cells
ppembed <- pppca$u %*% diag(pppca$d)
pcnames <- vector(mode = "character")
for (i in 1:dim(ppembed)[2])pcnames[i] <- paste("PC", i, sep = "_")
colnames(ppembed) <- pcnames
rownames(ppembed) <- rownames(seur.SCT)
sm.ppembed <- ppembed[1:10, 1:10]
ppldngs <- pppca$v
colnames(ppldngs) <- pcnames
rownames(ppldngs) <- mvf # assuming mvf has been subset to used for calc
sm.ppldngs <- ppldngs[1:10, 1:10]
ppembed[1:10, 1:10]
ppldngs[1:10, 1:10]
# 3/25/21:  these embeddings and loadings match Seur loadings and embeddings approximately, different n pcs?
so@reductions$pca@cell.embeddings <- ppembed
so@reductions$pca@feature.loadings <- ppldngs
so@reductions$pca@stdev <- pppca$d
print(Sys.time())

so <- RunHarmony(so, "orig.ident",
                 do_pca=FALSE,
                 assay.use = "SCT",
                 plot_convergence = TRUE,
                 return_object=TRUE)
head(so@reductions$harmony@cell.embeddings)
head(so@reductions$harmony@feature.loadings)

so <- RunUMAP(so, reduction = "harmony",dims=1:4)
so <- RunTSNE(so, reduction = "harmony",dims=1:4)


#    so <- RunHarmony(so, "orig.ident",
#                 do_pca=FALSE,
#                 plot_convergence = TRUE,
#                 return_object=TRUE)
#    print(Sys.time())
#    so <- RunUMAP(so, reduction = "harmony",dims=1:15)
#    so <- RunTSNE(so, reduction = "harmony",dims=1:15)

sdat <- data.frame(as.vector(so@reductions$tsne@cell.embeddings[,1]),
                   as.vector(so@reductions$tsne@cell.embeddings[,2]),
                   so@meta.data$orig.ident)
names(sdat) <- c("TSNE1","TSNE2","ident")

n <- 2e3
set.seed(10)
ourColorSpace <- colorspace::RGB(runif(n), runif(n), runif(n))
ourColorSpace <- as(ourColorSpace, "LAB")


distinctColorPalette <-function(k=1) {
  currentColorSpace <- ourColorSpace@coords
  # Set iter.max to 20 to avoid convergence warnings.
  set.seed(1)
  km <- kmeans(currentColorSpace, k, iter.max=20)
  colors <- unname(hex(LAB(km$centers)))
  return(colors)
}   

n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
qual_col_pals = qual_col_pals[c(7,6,2,1,8,3,4,5),]
cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


#q=30
#cols=distinctColorPalette(q)


g1 <- ggplot(sdat, aes(x=TSNE1, y=TSNE2)) +
  theme_bw() +
  theme(legend.title=element_blank()) +
  geom_point(aes(colour=ident),size=1) +
  scale_color_manual(values=cols) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position="top",
        panel.background = element_blank(),
        legend.text=element_text(size=rel(0.5))) +
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) 
print(g1)
rtnlist <- list(so, ppldngs)
return(rtnlist[[1]])

