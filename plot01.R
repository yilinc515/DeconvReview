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

X02 <- rbind.fill(res.music['16030X2', ], res.scdc['16030X2', ], 
                  res.bisque['16030X2', ], res.ted['16030X2', ], res.cibersortx['16030X2', ])

row.names(X02) <- c("MuSiC", "SCDC", "Bisque", "BayesPrism", "Cibersortx")

X03 <- rbind.fill(res.music['16030X3', ], res.scdc['16030X3', ], 
                  res.bisque['16030X3', ], res.ted['16030X3', ], res.cibersortx['16030X3', ])
row.names(X03) <- c("MuSiC", "SCDC", "Bisque", "BayesPrism", "Cibersortx")

X04 <- rbind.fill(res.music['16030X4', ], res.scdc['16030X4', ], 
                  res.bisque['16030X4', ], res.ted['16030X4', ], res.cibersortx['16030X4', ])
row.names(X04) <- c("MuSiC", "SCDC", "Bisque", "BayesPrism", "Cibersortx")


combined <- combine(X02, X03, X04)



combined <- t(combined)
combined <- as.data.frame(combined)






dat4 <- combined
dat4$base <- rownames(combined)
dat4 <- melt(dat4)


colnames(dat4) <- c("Samples", "Methods",  "Cell Types", "Proportions")

dat4$`Cell Types` <- factor(dat4$`Cell Types`, 
                            levels = c("Malignant.cells", "Fibroblasts",      
                                       "Monocytes", "Endothelial.cells",
                                       "B.cells", "Macrophages",      
                                       "T.cells"))
library(stringi)
dat4$Methods <- stri_sub(dat4$Methods,1, -2)

# plot
plot1.comb <- ggplot(dat4, aes(fill = `Cell Types`, y = Proportions, x = Methods)) + 
  geom_bar(stat = "identity") + 
  facet_grid(. ~ Samples) + 
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) + scale_fill_npg() 
  


