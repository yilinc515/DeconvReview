
load("sce02.RData")
load("sce03.RData")
load("sce04.RData")


#make a copy
sce02.shuffled <- sce

label02 <- colData(sce)$manual_annotation
label02 <- as.data.frame(label02)

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



#_________X03_________________

#make a copy
sce03.shuffled <- sce03

label03 <- colData(sce03)$manual_annotation
label03 <- as.data.frame(label03)

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
label04 <- as.data.frame(label04)

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
sc_ref <- as.data.frame(as.matrix(assays(sce02.shuffled)$counts))
# set row name to gene symbol
row.names(sc_ref) <- rowData(sce02.shuffled)$Symbol
# set col name to cell type annotation
colnames(sc_ref) <- sce02.shuffled$manual_annotation

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
load("sce02.RData")
load("sce03.RData")
load("sce04.RData")

# ------
# processing X02
# ------
sce.copy <- sce
# use Ensemble ID as rownames
rownames(sce.copy) <- rowData(sce.copy)$ID
# gene_exprs.matrix = assays(sce.copy)$counts
exprs02.matrix <- assays(sce.copy)$counts
pheno02.matrix <- colData(sce.copy)

# ------
# processing X03
# ------
sce03.copy <- sce03
# use Ensemble ID as rownames
rownames(sce03.copy) <- rowData(sce03.copy)$ID
# gene_exprs.matrix = assays(sce.copy)$counts
exprs03.matrix <- assays(sce03.copy)$counts
pheno03.matrix <- colData(sce03.copy)

# ------
# processing X04
# ------
sce04.copy <- sce04
# use Ensemble ID as rownames
rownames(sce04.copy) <- rowData(sce04.copy)$ID
# gene_exprs.matrix = assays(sce.copy)$counts
exprs04.matrix <- assays(sce04.copy)$counts
pheno04.matrix <- colData(sce04.copy)


# ------
# combine expression matrix since MuSiC starts with multi-subject scRNA-seq data
# ------
library(Seurat)
sc.exprs.matrix <- RowMergeSparseMatrices(exprs02.matrix, exprs03.matrix)
sc.exprs.matrix <- RowMergeSparseMatrices(sc.exprs.matrix, exprs04.matrix)


# ------
# combine phenotype matrix (colData)
# ------
sc.pheno.matrix <- rbind(pheno02.matrix, pheno03.matrix)
sc.pheno.matrix <- rbind(sc.pheno.matrix, pheno04.matrix)
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






