library(ggplot2)
library(reshape2)
library(tidyverse)
library(ggsci)
library(plyr)
library(SingleCellExperiment)


res <- load("res.RData")
load("sce02.RData")
load("sce03.RData")
load("sce04.RData")


txt<- "16030X3	0	0.360681564095981	0.0391082676450189	0.385662817321777	0.180065340160833	0.0333372666135224	0.00114474416286758	0	0.868540697970453	0.769336960751569
16030X2	0.329775579645059	0.448623421687672	0	0.190401637174546	0.0138969466601463	0	0.0173024148325771	0	0.722939847297778	0.844779642006333
16030X4	0	0.592459399290005	0	0.284067420004284	0.123473180705711	0	0	0	0.794266805033956	0.820486205235996"

res.cibersortx <- read.table(text = txt, row.names = 1)
colnames(res.cibersortx) <- c("Malignant cells", "Fibroblasts", "T-cells", 	"Macrophages", "Monocytes", "Endothelial cells", "B-cells", "P-value", "Correlation", "RMSE")
res.cibersortx <- res.cibersortx[1:(length(res.cibersortx)-3)]

res.bisque <- t(res.bisque)


res.music <- as.data.frame.matrix(res.music)
res.ted <- as.data.frame.matrix(res.ted)
res.bisque <- as.data.frame.matrix(res.bisque)
res.scdc <- as.data.frame.matrix(res.scdc)






#------------X03----------------
#-------------------------------
X03 <- rbind.fill(res.music['16030X3', ], res.scdc['16030X3', ], 
                  res.bisque['16030X3', ], res.ted['16030X3', ], res.cibersortx['16030X3', ])


row.names(X03) <- c("MuSiC", "SCDC", "Bisque", "BayesPrism", "Cibersortx")

X03 <- t(X03)
X03 <- as.data.frame(X03)


prop02 <- table(colData(sce)$manual_annotation)
prop02 <- prop02/sum(prop02)
prop02 <- prop02[match(rownames(X03), names(prop02))]
X03$Reference02 <- prop02

prop04 <- table(colData(sce04)$manual_annotation)
prop04 <- prop04/sum(prop04)
prop04 <- prop04[match(rownames(X03), names(prop04))]
X03$Reference04 <- prop04



dat <- X03
dat$base <- rownames(X03)
dat <- melt(dat)

colnames(dat) <- c("Cell Types", "Methods", "Proportions")

dat$`Cell Types` <- factor(dat$`Cell Types`, 
                           levels = c("Malignant cells", "Fibroblasts",      
                                      "Monocytes", "Endothelial cells",
                                      "B-cells", "Macrophages",      
                                      "T-cells"))

# plot
plot1.X03 <- ggplot(dat, aes(fill = `Cell Types`, y = Proportions, x = Methods)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) + 
  scale_fill_npg()



#---------------------------------------------
#---------------X02-------------------------------------
X02 <- rbind.fill(res.music['16030X2', ], res.scdc['16030X2', ], 
                  res.bisque['16030X2', ], res.ted['16030X2', ], res.cibersortx['16030X2', ])



row.names(X02) <- c("MuSiC", "SCDC", "Bisque", "BayesPrism", "Cibersortx")

X02 <- t(X02)
X02 <- as.data.frame(X02)




prop02 <- table(colData(sce)$manual_annotation)
prop02 <- prop02/sum(prop02)
prop02 <- prop02[match(rownames(X02), names(prop02))]
X02$Reference02 <- prop02

prop04 <- table(colData(sce04)$manual_annotation)
prop04 <- prop04/sum(prop04)
prop04 <- prop04[match(rownames(X02), names(prop04))]
X02$Reference04 <- prop04


library(ggplot2)
library(reshape2)
library(tidyverse)

dat2 <- X02
dat2$base <- rownames(X02)
dat2 <- melt(dat2)


colnames(dat2) <- c("Cell Types", "Methods", "Proportions")

dat2$`Cell Types` <- factor(dat2$`Cell Types`, 
                            levels = c("Malignant cells", "Fibroblasts",      
                                       "Monocytes", "Endothelial cells",
                                       "B-cells", "Macrophages",      
                                       "T-cells"))

# plot
plot1.X02 <- ggplot(dat2, aes(fill = `Cell Types`, y = Proportions, x = Methods)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) + 
  scale_fill_npg()




#----------------------------------------------
#-----------------------X04----------------------
#-----------------------------

X04 <- rbind.fill(res.music['16030X4', ], res.scdc['16030X4', ], 
                  res.bisque['16030X4', ], res.ted['16030X4', ], res.cibersortx['16030X4', ])


row.names(X04) <- c("MuSiC", "SCDC", "Bisque", "BayesPrism", "Cibersortx")

X04 <- t(X04)
X04 <- as.data.frame(X04)


prop02 <- table(colData(sce)$manual_annotation)
prop02 <- prop02/sum(prop02)
prop02 <- prop02[match(rownames(X04), names(prop02))]
X04$Reference02 <- prop02

prop04 <- table(colData(sce04)$manual_annotation)
prop04 <- prop04/sum(prop04)
prop04 <- prop04[match(rownames(X04), names(prop04))]
X04$Reference04 <- prop04



dat3 <- X04
dat3$base <- rownames(X04)
dat3 <- melt(dat3)


colnames(dat3) <- c("Cell Types", "Methods", "Proportions")

dat3$`Cell Types` <- factor(dat3$`Cell Types`, 
                            levels = c("Malignant cells", "Fibroblasts",      
                                       "Monocytes", "Endothelial cells",
                                       "B-cells", "Macrophages",      
                                       "T-cells"))

# plot
plot1.X04 <- ggplot(dat3, aes(fill = `Cell Types`, y = Proportions, x = Methods)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) + 
  scale_fill_npg()







#-------------
# Side by side
#-------------
library(gdata)
Reference <- data.frame(matrix(ncol = 7, nrow = 2))
rownames(Reference) <- c("Reference02", "Reference04")
colnames(Reference) <- c("Malignant.cells", "Fibroblasts",      
                         "Monocytes", "Endothelial.cells",
                         "B.cells", "Macrophages",      
                         "T.cells")

X02 <- rbind.fill(res.music['16030X2', ], res.scdc['16030X2', ], 
                  res.bisque['16030X2', ], res.ted['16030X2', ], res.cibersortx['16030X2', ])

row.names(X02) <- c("MuSiC", "SCDC", "Bisque", "BayesPrism", "Cibersortx")
colnames(X02) <- c("Malignant.cells", "Fibroblasts",      
                   "Monocytes", "Endothelial.cells",
                   "B.cells", "Macrophages",      
                   "T.cells")



X03 <- rbind.fill(res.music['16030X3', ], res.scdc['16030X3', ], 
                  res.bisque['16030X3', ], res.ted['16030X3', ], res.cibersortx['16030X3', ])
row.names(X03) <- c("MuSiC", "SCDC", "Bisque", "BayesPrism", "Cibersortx")
colnames(X03) <- c("Malignant.cells", "Fibroblasts",      
                   "Monocytes", "Endothelial.cells",
                   "B.cells", "Macrophages",      
                   "T.cells")



X04 <- rbind.fill(res.music['16030X4', ], res.scdc['16030X4', ], 
                  res.bisque['16030X4', ], res.ted['16030X4', ], res.cibersortx['16030X4', ])
row.names(X04) <- c("MuSiC", "SCDC", "Bisque", "BayesPrism", "Cibersortx")
colnames(X04) <- c("Malignant.cells", "Fibroblasts",      
                   "Monocytes", "Endothelial.cells",
                   "B.cells", "Macrophages",      
                   "T.cells")

combined <- combine(Reference, X02, X03, X04)

combined <- t(combined)
combined <- as.data.frame(combined)


prop02 <- table(colData(sce)$manual_annotation)
prop02 <- prop02/sum(prop02)
names(prop02) <- c("B.cells", "Endothelial.cells","Fibroblasts",      
                   "Macrophages", "Malignant.cells", "Monocytes",        
                   "T.cells")
prop02 <- prop02[match(rownames(combined), names(prop02))]
combined$Reference02 <- prop02

prop04 <- table(colData(sce04)$manual_annotation)
prop04 <- prop04/sum(prop04)
names(prop04) <- c("Endothelial.cells", "Fibroblasts", "Malignant.cells",  
                   "Monocytes" )
prop04 <- prop04[match(rownames(combined), names(prop04))]
combined$Reference04 <- prop04


combined[8, 1:2] <- c("Reference02", "Reference04")


dat4 <- combined
source <- dat4[8,]
dat4 <- dat4[1:7, ]
dat4$base <- rownames(combined)[1:7]

dat4 <- melt(dat4, id = c("base"))
dat4$source <- rep(source,each =7)


colnames(dat4) <- c( "Cell Types", "Methods",  "Proportions", "Samples")

dat4$`Cell Types` <- factor(dat4$`Cell Types`, 
                            levels = c("Malignant.cells", "Fibroblasts",      
                                       "Monocytes", "Endothelial.cells",
                                       "B.cells", "Macrophages",      
                                       "T.cells"))

dat4$Methods[15:119] <- gsub('[[:digit:]]+', '', dat4$Methods[15:119])
dat4$Samples[1:14] <- gsub('[[:digit:]]+', '', dat4$Samples[1:14])



dat4$Samples <- as.factor(as.character(unlist(dat4$Samples)))
dat4$Samples[1:7] <- rep("X02", 7)
dat4$Samples[8:14] <- rep("X04", 7)

#add reference03
prop03 <- table(colData(sce03)$manual_annotation)
prop03 <- prop03/sum(prop03)
prop03[4] <- prop03[4]+prop03[5]+prop03[6]
prop03 <- prop03[1:4]
prop03[5:7] <- NA         
names(prop03) <- c("Macrophages", "Monocytes", "T.cells", "Malignant.cells",
                   "B.cells", "Fibroblasts", "Endothelial.cells")
prop03 <- prop03[c(4, 6, 2, 7, 5, 1, 3)]

prop03 <- data.frame('Cell Types' = names(prop03), 
                     "Methods" = rep("Reference03", 7),
                     "Proportions" = matrix(prop03),
                     "Samples" = rep("X03", 7))
colnames(prop03) <- c("Cell Types", "Methods", "Proportions", "Samples")

dat5 <- rbind(prop03, dat4)
dat5$`Cell Types` <- factor(dat5$`Cell Types`, 
                            levels = c("Malignant.cells", "Fibroblasts",      
                                       "Monocytes", "Endothelial.cells",
                                       "B.cells", "Macrophages",      
                                       "T.cells"))

dat5$Methods  <- factor(dat5$Methods, 
                        levels = c("Reference02", "Reference03", "Reference04",      
                                   "Cibersortx", "BayesPrism",
                                   "Bisque", "MuSiC",      
                                   "SCDC"))



# plot
plot1.comb <- ggplot(dat5, aes(fill = `Cell Types`, y = Proportions, x = Methods)) + 
  geom_bar(stat = "identity") + 
  facet_grid(. ~ Samples, scales = "free_x") + 
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) + scale_fill_npg() 



