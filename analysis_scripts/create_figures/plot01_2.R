library(ggplot)
library(DropletUtils)
library(scater)
library(EnsDb.Hsapiens.v86)
library(scran)
library(BiocParallel)
library(SingleR)
library(pheatmap)
library(tidyverse)


library(TENxPBMCData)
pbmc <- TENxPBMCData(dataset = "pbmc8k")
rownames(pbmc)


# ------------
# QC
# assume cell calling already performed, see https://support.bioconductor.org/p/123554/#123562
# ------------

unfiltered <- pbmc


is.mito <- grep("MT", rowData(pbmc)$Symbol)
stats <- perCellQCMetrics(pbmc, subsets=list(Mito=is.mito))
pbmc <- addPerCellQC(pbmc, subsets = list(Mito=is.mito))
pbmc <- addPerFeatureQC(pbmc)

#Remove worst cells
initialFilt <- pbmc$sum > 500 & pbmc$detected > 100
length(which(!initialFilt))
pbmc<-pbmc[,initialFilt]


#high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
#pbmc <- pbmc[,!high.mito]



#mitochondrial read fraction and number of unique features (genes) with at least one read
metrics <- as.data.frame(colData(pbmc))
p <- ggplot(metrics, aes(x=detected,y=subsets_Mito_percent)) + geom_point()
p + stat_smooth(method="loess",formula=y~x,size=1,se=F,colour="blue")
p <- ggplot(metrics, aes(x=total_counts,y=subsets_Mito_percent)) + geom_point()
p <- ggplot(metrics, aes(x=detected,y=log10_total_counts)) + geom_point()

# mixture of two different distributions, the "functional" distribution and the "failed" distribution
library(flexmix)
options(scipen = 5)
set.seed(1010)
model<-flexmix(subsets_Mito_percent~detected,data=metrics,k=2)
model@components
slope_1=model@components$Comp.1[[1]]@parameters$coef[2]
intercept_1=model@components$Comp.1[[1]]@parameters$coef[1]
slope_2=model@components$Comp.2[[1]]@parameters$coef[2]
intercept_2=model@components$Comp.2[[1]]@parameters$coef[1]
ggplot(metrics, aes(x=detected,y=subsets_Mito_percent)) + geom_point() + 
  geom_abline(data=metrics,aes(slope=slope_1,intercept=intercept_1)) + 
  geom_abline(data=metrics,aes(slope=slope_2,intercept=intercept_2))

# posterior probability
post <- posterior(model)
metrics$posterior_dist1 <- post[,1]
metrics$posterior_dist2 <- post[,2]
ggplot(metrics, aes(x=detected,y=subsets_Mito_percent,colour=posterior_dist1)) + 
  geom_point() + geom_abline(data=metrics,aes(slope=slope_1,intercept=intercept_1)) + 
  geom_abline(data=metrics,aes(slope=slope_2,intercept=intercept_2))

# throw out cells with low posterior probability of coming from functional distribution
# cutoff = 0.25 
metrics$keep<-metrics$posterior_dist1>=0.25
ggplot(metrics, aes(x=detected,y=subsets_Mito_percent,colour=keep)) + geom_point()
table(metrics$keep)
pbmc <- pbmc[,metrics$keep]
dim(pbmc)


# check data after QC
qcplots <- list()
qcplots[[1]] <- plotColData(pbmc, x="sum", y="subsets_Mito_percent",
) + scale_x_log10()
qcplots[[2]] <-plotColData(pbmc, x="detected", y="subsets_Mito_percent",
) + scale_x_log10()
qcplots[[3]] <-plotColData(pbmc, x="detected", y="sum",
) + scale_y_log10()

do.call(gridExtra::grid.arrange, c(qcplots, ncol=3))
summary(high.mito)



# ------------
# Normalization
# ------------

set.seed(1000)
clusters <- quickCluster(pbmc)
pbmc <- computeSumFactors(pbmc, cluster=clusters)
pbmc <- logNormCounts(pbmc)
summary(sizeFactors(pbmc))

plot(librarySizeFactors(pbmc), sizeFactors(pbmc), pch=16, xlab="Library size factors", ylab="Deconvolution factors", log="xy")


# ------------
# Feature selection assuming near-Poisson variation
# ------------

# Quantifying per-gene variation
set.seed(1001)
pois <- modelGeneVarByPoisson(pbmc)
pois <- pois[order(pois$bio, decreasing=TRUE),]
head(pois)

plot(pois$mean, pois$total, pch=16, xlab="Mean of log-expression", ylab="Variance of log-expression")
curve(metadata(pois)$trend(x), col="dodgerblue", add=TRUE)

# Selecting highly variable genes
# select the top 10% of genes with the highest biological components
dec <- modelGeneVar(pbmc)
top.pbmc <- getTopHVGs(dec, prop=0.1)
str(top.pbmc)

# store hvg in AltExp
pbmc.hvg <- pbmc[top.pbmc,]
altExp(pbmc.hvg, "original") <- pbmc
altExpNames(pbmc.hvg)
# To recover original data: pbmc.original <- altExp(pbmc.hvg, "original", withColData=TRUE)





# ------------
# Dimensionality reductionn with PCA and NMF
# ------------
set.seed(1002) 
pbmc.hvg <- runPCA(pbmc.hvg) 
dim(reducedDim(pbmc.hvg, "PCA"))

#Choose the number of PCs Using the technical noise
set.seed(1003)
denoised.pbmc <- denoisePCA(pbmc.hvg, technical=dec)
ncol(reducedDim(denoised.pbmc))  # 5, elbow point gives 6

# OR: Based on population structure
# takes too long to run
#pcs <- reducedDim(pbmc.hvg)
#choices <- getClusteredPCs(pcs)
#metadata(choices)$chosen

# Set d = 30 for inital analysis
reducedDim(pbmc.hvg, "PCA") <- reducedDim(pbmc.hvg, "PCA")[,1:30]

# NMF  NOT WORKING!!!
set.seed(101001)
nmf.pbmc <- runNMF(pbmc.hvg, ncomponents=10, altexp = NULL)

nmf.out <- reducedDim(nmf.pbmc, "NMF")
nmf.basis <- attr(nmf.out, "basis")
colnames(nmf.out) <- colnames(nmf.basis) <- 1:10

per.cell <- pheatmap::pheatmap(nmf.out, silent=TRUE, 
                               main="By cell", show_rownames=FALSE,
                               color=rev(viridis::magma(100)), cluster_cols=FALSE) 

per.gene <- pheatmap::pheatmap(nmf.basis, silent=TRUE, 
                               main="By gene", cluster_cols=FALSE, show_rownames=FALSE,
                               color=rev(viridis::magma(100)))

gridExtra::grid.arrange(per.cell[[4]], per.gene[[4]], ncol=2)


# Visualization with UMAP
set.seed(1100)
pbmc.hvg <- runUMAP(pbmc.hvg, dimred="PCA")
plotReducedDim(pbmc.hvg, dimred="UMAP") # no available annotation



# ------------
# Clustering (Graph based)
# ------------
library(scran)
g02 <- buildSNNGraph(pbmc.hvg, k=20, use.dimred = 'PCA') # need to further decide what k to use
clust02 <- igraph::cluster_walktrap(g02)$membership
table(clust02)

# colLabels(pbmc.hvg) <- factor(clust)  THIS IS ONLY AVAILABLE IN pbmc 1.9.3, needs BioC-devel
pbmc$label <- factor(clust02) # label the original dataset
pbmc.hvg$label <- factor(clust02)
plotReducedDim(pbmc.hvg, dimred = "UMAP", colour_by="label")  
# plot by TSNE gives error "Error in `rownames<-`(`*tmp*`, value = c("AAACCCAAGCCACCGT-1", "AAACCCAAGGATGGCT-1",  : 
# attempt to set 'rownames' on an object with no dimensions"


#  Assessing cluster separation
ratio <- clusterModularity(g, clust, as.ratio=TRUE)
dim(ratio)
pheatmap(log2(ratio+1), cluster_rows=FALSE, cluster_cols=FALSE,
         color=colorRampPalette(c("white", "blue"))(100))

# Evaluating cluster stability
myClusterFUN <- function(x) {
  g <- buildSNNGraph(x, use.dimred="PCA", type="jaccard")
  igraph::cluster_louvain(g)$membership
}

originals <- myClusterFUN(pbmc.hvg)

set.seed(001001)
coassign <- bootstrapCluster(pbmc.hvg, FUN=myClusterFUN, clusters=originals)
dim(coassign)

pheatmap(coassign, cluster_row=FALSE, cluster_col=FALSE,
         color=rev(viridis::magma(100)))




