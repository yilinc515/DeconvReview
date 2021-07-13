
library(SingleCellExperiment)

load("sce02.RData")
load("sce03.RData")
load("sce04.RData")


#make a copy

func1 <- function(x)
  if (x == "6" || x == "16") {
    "B-cells"
  } else if (x == "4" || x == "17") {
    "Macrophages"
  } else if (x == "9") {
    "Monocytes"
  } else if (x == "8" || x == "3") {
    "Fibroblasts"
  } else if (x == "7") {
    "Endothelial cells"
  } else if (x == "13") {
    "T-cells"
  } else {
    "Malignant cells" # does not keep clustering label
  }
colData(sce)$manual_annotation <- mapply(func1, sce$label)

sce02.shuffled <- sce
label02 <- colData(sce02.shuffled)$manual_annotation
label02 <- unlist(label02)

set.seed(5)
indices02 <-  sample(1:4737, 237, replace=F)
set.seed(5)
shuffled02 <- sample(indices02)

# shuffling label02 
for(i in 1:length(indices02)){
  label02[i] <- colData(sce)$manual_annotation[shuffled02[i]]
}
# set manual annotation
colData(sce02.shuffled)$manual_annotation <- label02




sce02.hvg.shuffled <- sce.hvg
func1 <- function(x)
  if (x == "6" || x == "16") {
    "B-cells"
  } else if (x == "4" || x == "17") {
    "Macrophages"
  } else if (x == "9") {
    "Monocytes"
  } else if (x == "8" || x == "3") {
    "Fibroblasts"
  } else if (x == "7") {
    "Endothelial cells"
  } else if (x == "13") {
    "T-cells"
  } else {
    "Malignant cells" # does not keep clustering label
  }
colData(sce02.hvg.shuffled)$manual_annotation <- mapply(func1, sce.hvg$label)


label02 <- colData(sce02.hvg.shuffled)$manual_annotation
label02 <- unlist(label02)

set.seed(5)
indices02 <-  sample(1:4737, 237, replace=F)
set.seed(5)
shuffled02 <- sample(indices02)

# shuffling label02 
for(i in 1:length(indices02)){
  label02[i] <- colData(sce)$manual_annotation[shuffled02[i]]
}

# set manual annotation
colData(sce02.hvg.shuffled)$manual_annotation <- label02


#_________X03_________________

#make a copy
sce03.shuffled <- sce03

label03 <- colData(sce03)$manual_annotation
label03 <- unlist(label03)

set.seed(6)
indices03 <-  sample(1:950, 48, replace=F)
set.seed(6)
shuffled03 <- sample(indices03)

# shuffling label003
for(i in 1:length(indices03)){
  label03[i] <- colData(sce03)$manual_annotation[shuffled03[i]]
}

# set manual annotation
colData(sce03.shuffled)$manual_annotation <- label03



#_________X04_________________

#make a copy
sce04.shuffled <- sce04

label04 <- colData(sce04)$manual_annotation
label04 <- unlist(label04)

set.seed(7)
indices04 <-  sample(1:4087, 205, replace=F)
# shuffle the selected indices
set.seed(7)
shuffled04 <- sample(indices04)

# shuffling label04
for(i in 1:length(indices04)){
  label04[i] <- colData(sce04)$manual_annotation[shuffled04[i]]
}

# set manual annotation
colData(sce04.shuffled)$manual_annotation <- label04





#--------------------------------Cibersortx-------------------------
library(scater)
library(scran)






# convert sparse matrix to data.frame
sc_ref <- as.data.frame(as.matrix(assays(sce02.hvg.shuffled)$counts))
# set row name to gene symbol
row.names(sc_ref) <- rowData(sce02.hvg.shuffled)$Symbol
# set col name to cell type annotation
colnames(sc_ref) <- sce02.hvg.shuffled$manual_annotation

# write to tab-delimited .txt
write.table(sc_ref, file = "cibersortx_sc02shuffled_ref.txt", quote = FALSE, sep = "\t", col.names = NA)

#-----------------
# Mixture file
#-----------------

library(EnsDb.Hsapiens.v86)


# load bulk data
load("gene_counts_17667X1-3.RData")  #gse

ensembl.ids <- rownames(gse)

# note that gse's Ensemble ID has dot suffix as version number
# in order to unify with sc data, strip the version number in rownames to use just the Ensembl ID
library(stringr)
ensembl.ids <- str_replace(ensembl.ids, pattern = ".[0-9]+$", replacement = "")
# convert ensembl ids to gene symbols
gene.symbols <- mapIds(EnsDb.Hsapiens.v86, keys = ensembl.ids, keytype = "GENEID", column="SYMBOL")

# build bulk reference data frame
bulk_ref <- assays(gse)$counts
row.names(bulk_ref) <- gene.symbols

# unify colnames with scRNA-seq data
colnames(bulk_ref) <- c("16030X3", "16030X2", "16030X4")

# write to tab-delimited .txt
write.table(bulk_ref, file = "cibersortx_bulk.txt", quote = FALSE, sep = "\t", col.names = NA)



#Result
txt<- "16030X3,0.364543405390944,0.0560103628195528,0.000671363890093338,0.109160015937573,0.232876303258855,0.210290961075759,0.0264475876272227,0,0.902930064609702,0.833686023279097
16030X2,0.784344291153795,0.00130177288283927,0.00818138488152475,0,0.0618921475833917,0.144280403498449,0,0,0.839952909302621,0.907987242920674
16030X4,0.610970900669865,0,0,0.0158459116016397,0.12065933209891,0.252523855629585,0,0,0.782496245771328,0.901148903913454"

res.shuffled.cibersortx <- read.table(text = txt, header=F,sep = ",", row.names = 1)
                                      
colnames(res.shuffled.cibersortx) <- c("Malignant cells", "T-cells", "B-cells", "Monocytes", "Macrophages",	"Fibroblasts", "Endothelial cells" ,"P-value", "Correlation", "RMSE")
res.shuffled.cibersortx <- res.shuffled.cibersortx[1:(length(res.shuffled.cibersortx)-3)]





#----------------Bisque------------------------

# ------
# load bulk data
# ------
load("gene_counts_17667X1-3.RData")

# bulk sample name | single-cell sample name
# ------------------------------------------
# 17667X1          | 16030X3
# 17667X2          | 16030X2
# 17667X3          | 16030X4

# unify colnames with scRNA-seq data
colnames(gse) <- c("16030X3", "16030X2", "16030X4")

# note that gse's Ensemble ID has dot suffix as version number
# in order to unify with sc data, strip the version number in rownames to use just the Ensembl ID
library(stringr)
rownames(gse) <- str_replace(rownames(gse), pattern = ".[0-9]+$", replacement = "")

# ------
# convert bulk data to ExpressionSet 
# https://github.com/xuranw/MuSiC/issues/2
# ------
metadata <- data.frame(labelDescription= c("names"), row.names=c("names"))
# gene_exprs.matrix = assays(gse)$counts, pheno.matrix = colData(gse)
# convert DFrame to data.frame
pheno.matrix <- as.data.frame(colData(gse))
# construct ExpressionSet of bulk data
bulk.eset = ExpressionSet(assayData = data.matrix(assays(gse)$counts), phenoData =  new("AnnotatedDataFrame", data = pheno.matrix, varMetadata = metadata))



# ------
# load single-cell data 
# Note that Bisque won't run if all samples have sc and bulk, so X03 is taken out
# ------

# ------
# processing X02
# ------
sce02.copy <- sce02.shuffled
# use Ensemble ID as rownames
rownames(sce02.copy) <- rowData(sce02.copy)$ID
# gene_exprs.matrix = assays(sce.copy)$counts
exprs02.matrix <- assays(sce02.copy)$counts
pheno02.matrix <- colData(sce02.copy)
# ------
# processing X04
# ------
sce04.copy <- sce04.shuffled
# use Ensemble ID as rownames
rownames(sce04.copy) <- rowData(sce04.copy)$ID
# gene_exprs.matrix = assays(sce.copy)$counts
exprs04.matrix <- assays(sce04.copy)$counts
pheno04.matrix <- colData(sce04.copy)


# ------
# combine expression matrix since MuSiC starts with multi-subject scRNA-seq data
# ------
library(Seurat)
# X02+X04
sc.exprs.matrix <- RowMergeSparseMatrices(exprs02.matrix, exprs04.matrix)

# make sure the colnames are unique
colnames(sc.exprs.matrix) <- make.names(colnames(sc.exprs.matrix))

# ------
# combine phenotype matrix (colData)
# ------
# X02+X04
sc.pheno.matrix <- rbind(pheno02.matrix, pheno04.matrix)

# convert DFrame to data.frame
sc.pheno.matrix <- as.data.frame(sc.pheno.matrix)

# ------
# Generate expressionset of sc data
# https://github.com/xuranw/MuSiC/issues/2
# ------

# use make.names() to avoid duplicate row.names
sc.pheno <- data.frame(check.names=F, check.rows=F,
                       stringsAsFactors=F,
                       row.names=make.names(sc.pheno.matrix$Barcode, unique = TRUE), 
                       SubjectName=sc.pheno.matrix$Sample,
                       cellType=sc.pheno.matrix$manual_annotation)

sc.meta <- data.frame(labelDescription=c("SubjectName",
                                         "cellType"),
                      row.names=c("SubjectName",
                                  "cellType"))
sc.eset = ExpressionSet(assayData = data.matrix(sc.exprs.matrix), phenoData =  new("AnnotatedDataFrame", data = sc.pheno, varMetadata = sc.meta))


#----------------------------------------------------------------------
# Run Bisque deconvolution
# Use Reference-based decomposition mode
#----------------------------------------------------------------------
library(Biobase)
library(BisqueRNA)

# without markers: By default, Bisque uses all genes for decomposition.
# may supply a list of genes (such as marker genes) to be used with the markers parameter, see vignette
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset)
# A list is returned with decomposition estimates in slot bulk.props.

ref.based.estimates <- res$bulk.props
knitr::kable(ref.based.estimates, digits=2)
res.shuffled.bisque <- ref.based.estimates





#-----------------------MuSiC-----------------------------

# ------
# load bulk data 
# ------
load("gene_counts_17667X1-3.RData")

# bulk sample name | single-cell sample name
# ------------------------------------------
# 17667X1          | 16030X3
# 17667X2          | 16030X2
# 17667X3          | 16030X4

# unify colnames with scRNA-seq data
colnames(gse) <- c("16030X3", "16030X2", "16030X4")

# note that gse's Ensemble ID has dot suffix as version number
# in order to unify with sc data, strip the version number in rownames to use just the Ensembl ID
library(stringr)
rownames(gse) <- str_replace(rownames(gse), pattern = ".[0-9]+$", replacement = "")

# ------
# convert bulk data to ExpressionSet 
# https://github.com/xuranw/MuSiC/issues/2
# ------
metadata <- data.frame(labelDescription= c("names"), row.names=c("names"))
# gene_exprs.matrix = assays(gse)$counts, pheno.matrix = colData(gse)
# convert DFrame to data.frame
pheno.matrix <- as.data.frame(colData(gse))
# construct ExpressionSet of bulk data
bulk.eset = ExpressionSet(assayData = data.matrix(assays(gse)$counts), phenoData =  new("AnnotatedDataFrame", data = pheno.matrix, varMetadata = metadata))


# ------
# load single-cell data 
# ------

# ------
# processing X02
# ------
sce.copy <- sce02.shuffled
# use Ensemble ID as rownames
rownames(sce.copy) <- rowData(sce.copy)$ID
# gene_exprs.matrix = assays(sce.copy)$counts
exprs02.matrix <- assays(sce.copy)$counts
pheno02.matrix <- colData(sce.copy)

# ------
# processing X03
# ------
sce03.copy <- sce03.shuffled
# use Ensemble ID as rownames
rownames(sce03.copy) <- rowData(sce03.copy)$ID
# gene_exprs.matrix = assays(sce.copy)$counts
exprs03.matrix <- assays(sce03.copy)$counts
pheno03.matrix <- colData(sce03.copy)

# ------
# processing X04
# ------
sce04.copy <- sce04.shuffled
# use Ensemble ID as rownames
rownames(sce04.copy) <- rowData(sce04.copy)$ID
# gene_exprs.matrix = assays(sce.copy)$counts
exprs04.matrix <- assays(sce04.copy)$counts
pheno04.matrix <- colData(sce04.copy)


# ------
# combine expression matrix since MuSiC starts with multi-subject scRNA-seq data
# ------
library(Seurat)
sc.exprs.matrix <- RowMergeSparseMatrices(exprs02.matrix, exprs04.matrix)
#sc.exprs.matrix <- RowMergeSparseMatrices(sc.exprs.matrix, exprs04.matrix)


# ------
# combine phenotype matrix (colData)
# ------
sc.pheno.matrix <- rbind(pheno02.matrix, pheno04.matrix)
#sc.pheno.matrix <- rbind(sc.pheno.matrix, pheno04.matrix)
# convert DFrame to data.frame
sc.pheno.matrix <- as.data.frame(sc.pheno.matrix)

# ------
# Generate expressionset of sc data
# https://github.com/xuranw/MuSiC/issues/2
# ------
sc.metadata <- data.frame(labelDescription= c("Sample", "Barcode", "sum", "detected", "percent_top_50", "percent_top_100", "percent_top_200", "percent_top_500", 
                                              "subsets_Mito_sum", "subsets_Mito_detected", "subsets_Mito_percent", "total", "sizeFactor", "label", "manual_annotation"), 
                          row.names=c("Sample", "Barcode", "sum", "detected", "percent_top_50", "percent_top_100", "percent_top_200", "percent_top_500", 
                                      "subsets_Mito_sum", "subsets_Mito_detected", "subsets_Mito_percent", "total", "sizeFactor", "label", "manual_annotation"))
sc.eset = ExpressionSet(assayData = data.matrix(sc.exprs.matrix), phenoData =  new("AnnotatedDataFrame", data = sc.pheno.matrix, varMetadata = sc.metadata))




#----------------------------------------------------------------------
# RUN MuSiC
# https://xuranw.github.io/MuSiC/articles/MuSiC.html
# see sample anaylsis
#----------------------------------------------------------------------
library(MuSiC)
library(xbioc) # for pVar()
# run!
Est.prop = music_prop(bulk.eset = bulk.eset, sc.eset = sc.eset, clusters = 'manual_annotation', samples = 'Sample')

res.shuffled.music <- Est.prop$Est.prop.weighted







save(res.shuffled.bisque,res.shuffled.music, res.shuffled.cibersortx, file="resShuffled.RData")

#----------------------------SCDC--------------------------

load("resShuffled.RData")
# ------
# load bulk data
# ------
load("gene_counts_17667X1-3.RData")

# bulk sample name | single-cell sample name
# ------------------------------------------
# 17667X1          | 16030X3
# 17667X2          | 16030X2
# 17667X3          | 16030X4

# unify colnames with scRNA-seq data
colnames(gse) <- c("16030X3", "16030X2", "16030X4")

# note that gse's Ensemble ID has dot suffix as version number
# in order to unify with sc data, strip the version number in rownames to use just the Ensembl ID
library(stringr)
rownames(gse) <- str_replace(rownames(gse), pattern = ".[0-9]+$", replacement = "")

# ------
# convert bulk data to ExpressionSet 
# https://github.com/xuranw/MuSiC/issues/2
# ------
metadata <- data.frame(labelDescription= c("names"), row.names=c("names"))
# gene_exprs.matrix = assays(gse)$counts, pheno.matrix = colData(gse)
# convert DFrame to data.frame
pheno.matrix <- as.data.frame(colData(gse))
# construct ExpressionSet of bulk data
bulk.eset = ExpressionSet(assayData = data.matrix(assays(gse)$counts), phenoData =  new("AnnotatedDataFrame", data = pheno.matrix, varMetadata = metadata))



# ------
# load single-cell data 

# ------
# processing X02
# ------
sce.copy <- sce02.shuffled
# use Ensemble ID as rownames
rownames(sce.copy) <- rowData(sce.copy)$ID
# gene_exprs.matrix = assays(sce.copy)$counts
exprs02.matrix <- assays(sce.copy)$counts
pheno02.matrix <- colData(sce.copy)
pheno02.matrix <- as.data.frame(pheno02.matrix)

# use make.names() to avoid duplicate row.names
sc02.pheno <- data.frame(check.names=F, check.rows=F,
                         stringsAsFactors=F,
                         row.names=pheno02.matrix$Barcode, 
                         SubjectName=pheno02.matrix$Sample,
                         cellType=pheno02.matrix$manual_annotation)

sc02.meta <- data.frame(labelDescription=c("SubjectName",
                                           "cellType"),
                        row.names=c("SubjectName",
                                    "cellType"))
sc02.eset = ExpressionSet(assayData = data.matrix(exprs02.matrix), phenoData =  new("AnnotatedDataFrame", data = sc02.pheno, varMetadata = sc02.meta))


# ------
# processing X03
# ------
sce03.copy <- sce03.shuffled
# use Ensemble ID as rownames
rownames(sce03.copy) <- rowData(sce03.copy)$ID
# gene_exprs.matrix = assays(sce.copy)$counts
exprs03.matrix <- assays(sce03.copy)$counts
pheno03.matrix <- colData(sce03.copy)
pheno03.matrix <- as.data.frame(pheno03.matrix)

# use make.names() to avoid duplicate row.names
sc03.pheno <- data.frame(check.names=F, check.rows=F,
                         stringsAsFactors=F,
                         row.names=pheno03.matrix$Barcode, 
                         SubjectName=pheno03.matrix$Sample,
                         cellType=pheno03.matrix$manual_annotation)

sc03.meta <- data.frame(labelDescription=c("SubjectName",
                                           "cellType"),
                        row.names=c("SubjectName",
                                    "cellType"))
sc03.eset = ExpressionSet(assayData = data.matrix(exprs03.matrix), phenoData =  new("AnnotatedDataFrame", data = sc03.pheno, varMetadata = sc03.meta))


# ------
# processing X04
# ------
sce04.copy <- sce04.shuffled
# use Ensemble ID as rownames
rownames(sce04.copy) <- rowData(sce04.copy)$ID
# gene_exprs.matrix = assays(sce.copy)$counts
exprs04.matrix <- assays(sce04.copy)$counts
pheno04.matrix <- colData(sce04.copy)
pheno04.matrix <- as.data.frame(pheno04.matrix)

# use make.names() to avoid duplicate row.names
sc04.pheno <- data.frame(check.names=F, check.rows=F,
                         stringsAsFactors=F,
                         row.names=pheno04.matrix$Barcode, 
                         SubjectName=pheno04.matrix$Sample,
                         cellType=pheno04.matrix$manual_annotation)

sc04.meta <- data.frame(labelDescription=c("SubjectName",
                                           "cellType"),
                        row.names=c("SubjectName",
                                    "cellType"))
sc04.eset = ExpressionSet(assayData = data.matrix(exprs04.matrix), phenoData =  new("AnnotatedDataFrame", data = sc04.pheno, varMetadata = sc04.meta))



# ------
# combine expression matrix since MuSiC starts with multi-subject scRNA-seq data
# ------
library(Seurat)
sc.exprs.matrix <- RowMergeSparseMatrices(exprs02.matrix, exprs04.matrix)
#sc.exprs.matrix <- RowMergeSparseMatrices(sc.exprs.matrix, exprs04.matrix)

# make sure the colnames are unique
colnames(sc.exprs.matrix) <- make.names(colnames(sc.exprs.matrix))

# ------
# combine phenotype matrix (colData)
# ------
sc.pheno.matrix <- rbind(pheno02.matrix, pheno04.matrix)
#sc.pheno.matrix <- rbind(sc.pheno.matrix, pheno04.matrix)

# convert DFrame to data.frame
sc.pheno.matrix <- as.data.frame(sc.pheno.matrix)



# ------
# Generate expressionset of sc data
# https://github.com/xuranw/MuSiC/issues/2
# ------

# use make.names() to avoid duplicate row.names
sc.pheno <- data.frame(check.names=F, check.rows=F,
                       stringsAsFactors=F,
                       row.names=make.names(sc.pheno.matrix$Barcode, unique = TRUE), 
                       SubjectName=sc.pheno.matrix$Sample,
                       cellType=sc.pheno.matrix$manual_annotation)

sc.meta <- data.frame(labelDescription=c("SubjectName",
                                         "cellType"),
                      row.names=c("SubjectName",
                                  "cellType"))
sc.eset = ExpressionSet(assayData = data.matrix(sc.exprs.matrix), phenoData =  new("AnnotatedDataFrame", data = sc.pheno, varMetadata = sc.meta))



#----------------------------------------------------------------------
# RUN SCDC
# https://meichendong.github.io/SCDC/articles/SCDC.html#data-input-1
# see sample anaylsis
#----------------------------------------------------------------------
library(SCDC)
library(xbioc)

# -----
# SCDC Pre-process of scRNA-seq Data
# -----
sc02.qc <- SCDC_qc_ONE(sc02.eset, ct.varname = "cellType", sample = "SubjectName", scsetname = "X02",
                       ct.sub = c("Macrophages", "Monocytes", "Fibroblasts",  "Endothelial cells", "T-cells", "B-cells", "Malignant cells"), qcthreshold = 0.6) 

sc03.qc <- SCDC_qc_ONE(sc03.eset, ct.varname = "cellType", sample = "SubjectName", scsetname = "X03",
                       ct.sub = c("Macrophages", "Monocytes", "Malignant cells"), qcthreshold = 0.6)

load("qc.RData")
sc04.qc <- SCDC_qc_ONE(sc04.eset, ct.varname = "cellType", sample = "SubjectName", scsetname = "X04",
                       ct.sub = c("Monocytes", "Fibroblasts",  "Endothelial cells", "Malignant cells"), qcthreshold = 0.6) 



bulk.ens <- SCDC_ENSEMBLE(bulk.eset = bulk.eset, sc.eset.list = list(X02 = sc02.qc$sc.eset.qc,  X04 = sc04.qc$sc.eset.qc), ct.varname = "cellType",
                          sample = "SubjectName", truep = NULL, ct.sub =  c("Monocytes", "Fibroblasts",  "Endothelial cells", "Malignant cells"), search.length = 0.01, grid.search = T)  



bulk.ens$prop.est
# calculate cell prop from  linear combination of a list of proportions
res.shuffled.scdc <- wt_prop(bulk.ens$w_table[1, 1:2], bulk.ens$prop.only)


save(res.shuffled.scdc, res.shuffled.bisque,res.shuffled.music, res.shuffled.cibersortx, file="resShuffled.RData")





#-------------------BayesPrism---------------------
# ------
# load bulk data
# ------
load("resShuffled.RData")
load("gene_counts_17667X1-3.RData")

# bulk sample name | single-cell sample name
# ------------------------------------------
# 17667X1          | 16030X3
# 17667X2          | 16030X2
# 17667X3          | 16030X4

# unify colnames with scRNA-seq data
colnames(gse) <- c("16030X3", "16030X2", "16030X4")

# note that gse's Ensemble ID has dot suffix as version number
# in order to unify with sc data, strip the version number in rownames to use just the Ensembl ID
library(stringr)
library(scran)
rownames(gse) <- str_replace(rownames(gse), pattern = ".[0-9]+$", replacement = "")

# ------
# convert bulk data to ExpressionSet 
# https://github.com/xuranw/MuSiC/issues/2
# ------
metadata <- data.frame(labelDescription= c("names"), row.names=c("names"))
# gene_exprs.matrix = assays(gse)$counts, pheno.matrix = colData(gse)
# convert DFrame to data.frame
pheno.matrix <- as.data.frame(colData(gse))
# construct ExpressionSet of bulk data
bulk.eset = ExpressionSet(assayData = data.matrix(assays(gse)$counts), phenoData =  new("AnnotatedDataFrame", data = pheno.matrix, varMetadata = metadata))




# ------
# processing X02
# ------
sce.copy <- sce02.shuffled
# use Ensemble ID as rownames
rownames(sce.copy) <- rowData(sce.copy)$ID
# gene_exprs.matrix = assays(sce.copy)$counts
exprs02.matrix <- assays(sce.copy)$counts
pheno02.matrix <- colData(sce.copy)
pheno02.matrix <- as.data.frame(pheno02.matrix)

# use make.names() to avoid duplicate row.names
sc02.pheno <- data.frame(check.names=F, check.rows=F,
                         stringsAsFactors=F,
                         row.names=pheno02.matrix$Barcode, 
                         SubjectName=pheno02.matrix$Sample,
                         cellType=pheno02.matrix$manual_annotation)

sc02.meta <- data.frame(labelDescription=c("SubjectName",
                                           "cellType"),
                        row.names=c("SubjectName",
                                    "cellType"))
sc02.eset = ExpressionSet(assayData = data.matrix(exprs02.matrix), phenoData =  new("AnnotatedDataFrame", data = sc02.pheno, varMetadata = sc02.meta))


# ------
# processing X03
# ------
sce03.copy <- sce03.shuffled
# use Ensemble ID as rownames
rownames(sce03.copy) <- rowData(sce03.copy)$ID
# gene_exprs.matrix = assays(sce.copy)$counts
exprs03.matrix <- assays(sce03.copy)$counts
pheno03.matrix <- colData(sce03.copy)
pheno03.matrix <- as.data.frame(pheno03.matrix)

# use make.names() to avoid duplicate row.names
sc03.pheno <- data.frame(check.names=F, check.rows=F,
                         stringsAsFactors=F,
                         row.names=pheno03.matrix$Barcode, 
                         SubjectName=pheno03.matrix$Sample,
                         cellType=pheno03.matrix$manual_annotation)

sc03.meta <- data.frame(labelDescription=c("SubjectName",
                                           "cellType"),
                        row.names=c("SubjectName",
                                    "cellType"))
sc03.eset = ExpressionSet(assayData = data.matrix(exprs03.matrix), phenoData =  new("AnnotatedDataFrame", data = sc03.pheno, varMetadata = sc03.meta))


# ------
# processing X04
# ------
sce04.copy <- sce04.shuffled
# use Ensemble ID as rownames
rownames(sce04.copy) <- rowData(sce04.copy)$ID
# gene_exprs.matrix = assays(sce.copy)$counts
exprs04.matrix <- assays(sce04.copy)$counts
pheno04.matrix <- colData(sce04.copy)
pheno04.matrix <- as.data.frame(pheno04.matrix)

# use make.names() to avoid duplicate row.names
sc04.pheno <- data.frame(check.names=F, check.rows=F,
                         stringsAsFactors=F,
                         row.names=pheno04.matrix$Barcode, 
                         SubjectName=pheno04.matrix$Sample,
                         cellType=pheno04.matrix$manual_annotation)

sc04.meta <- data.frame(labelDescription=c("SubjectName",
                                           "cellType"),
                        row.names=c("SubjectName",
                                    "cellType"))
sc04.eset = ExpressionSet(assayData = data.matrix(exprs04.matrix), phenoData =  new("AnnotatedDataFrame", data = sc04.pheno, varMetadata = sc04.meta))



# ------
# combine expression matrix since MuSiC starts with multi-subject scRNA-seq data
# ------
library(Seurat)
#sc.exprs.matrix <- RowMergeSparseMatrices(exprs02.matrix, exprs03.matrix)
sc.exprs.matrix <- RowMergeSparseMatrices(exprs02.matrix, exprs04.matrix)
#sc.exprs.matrix <- RowMergeSparseMatrices(exprs03.matrix, exprs04.matrix)
#sc.exprs.matrix <- RowMergeSparseMatrices(sc.exprs.matrix, exprs04.matrix)

# make sure the colnames are unique
colnames(sc.exprs.matrix) <- make.names(colnames(sc.exprs.matrix))

# ------
# combine phenotype matrix (colData)
# ------
#sc.pheno.matrix <- rbind(pheno02.matrix, pheno03.matrix)
sc.pheno.matrix <- rbind(pheno02.matrix, pheno04.matrix)
#sc.pheno.matrix <- rbind(pheno03.matrix, pheno04.matrix)
#sc.pheno.matrix <- rbind(sc.pheno.matrix, pheno04.matrix)

# convert DFrame to data.frame
sc.pheno.matrix <- as.data.frame(sc.pheno.matrix)



# ------
# Generate expressionset of sc data
# https://github.com/xuranw/MuSiC/issues/2
# ------

# use make.names() to avoid duplicate row.names
sc.pheno <- data.frame(check.names=F, check.rows=F,
                       stringsAsFactors=F,
                       row.names=make.names(sc.pheno.matrix$Barcode, unique = TRUE), 
                       SubjectName=sc.pheno.matrix$Sample,
                       cellType=sc.pheno.matrix$manual_annotation)

sc.meta <- data.frame(labelDescription=c("SubjectName",
                                         "cellType"),
                      row.names=c("SubjectName",
                                  "cellType"))
sc.eset = ExpressionSet(assayData = data.matrix(sc.exprs.matrix), phenoData =  new("AnnotatedDataFrame", data = sc.pheno, varMetadata = sc.meta))






#----------------------------------------------------------------------
# RUN BayesPrism
# https://github.com/Danko-Lab/TED/blob/master/vignette.pdf
# see sample anaylsis
#----------------------------------------------------------------------
library(TED)

#  Each row is a cell (input.type="scRNA") or a cell type (input.type="GEP"); 
#  each column is a gene. rownames of ref.dat are cell IDs (input.type="scRNA") or cell type names(input.type="GEP");
# colnames of ref.dat are gene IDs.


# denote ambiguous as tumor
sc.pheno.matrix$manual_annotation[grepl("Ambiguous", sc.pheno.matrix$manual_annotation)] <- "Malignant cells"
pheno03.matrix$manual_annotation[grepl("Ambiguous", pheno03.matrix$manual_annotation)] <- "Malignant cells"
# Use X02 and X04 as a reference
test.ted <-  run.Ted(ref.dat= t(exprs(sc.eset)), 
                         X=t(exprs(bulk.eset)),
                         cell.type.labels=sc.pheno.matrix$manual_annotation,
                         cell.subtype.labels= NULL,
                         tum.key=c("Malignant cells"),
                         input.type="scRNA",
                         n.cores=45,
                         pdf.name="BayesPrismTest")
#save.image(test.ted, "TEDData03.RData")

# To see results:
#load("BPresults.RData")

#--------
# Datasets : test02.ted, test03.ted, test04.ted, testAll.ted
#  save(test02.ted, test03.ted, test04.ted, test0203.ted, testAll.ted, file = "BPresults.RData")
#--------

# #the correlation heatmap of tumor expressions across bulk samples
test.ted$res$cor.mat
# final updated cell type fraction estimate
res.shuffled.bayesprism <- test.ted$res$ final.gibbs.theta


save(res.shuffled.bayesprism, res.shuffled.scdc, res.shuffled.bisque,res.shuffled.music, res.shuffled.cibersortx, file="resShuffled.RData")










