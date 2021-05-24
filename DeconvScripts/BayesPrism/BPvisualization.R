library(TED)
# To see results:
load("TEDresults.RData")
load("TED.RData")
load("TED_repeat.RData")

#--------
# Datasets : test02.ted, test03.ted, test04.ted, testAll.ted
#  save(test02.ted, test03.ted, test04.ted, test0203.ted, test0204.ted, test0304.ted, testAll.ted, file = "BPresults.RData")
#--------

# #the correlation heatmap of tumor expressions across bulk samples
test.ted$res$cor.mat
# final updated cell type fraction estimate
test.ted$res$ final.gibbs.theta

# ----------------------------------
# BP results 
# data visualization

result02 <- test02.ted$res$ final.gibbs.theta
result03 <- test03.ted$res$ final.gibbs.theta
result04 <- test04.ted$res$ final.gibbs.theta
result0203 <- test0203.ted$res$ final.gibbs.theta
result0304 <- test0304.ted$res$ final.gibbs.theta
result0204 <- test0204.ted$res$ final.gibbs.theta

result0203_2 <- test0203_2.ted$res$ final.gibbs.theta
result0304_2 <- test0304_2.ted$res$ final.gibbs.theta
result0204_2 <- test0204_2.ted$res$ final.gibbs.theta

result0203_3 <- test0203_3.ted$res$ final.gibbs.theta
result0304_3 <- test0304_3.ted$res$ final.gibbs.theta
result0204_3 <- test0204_3.ted$res$ final.gibbs.theta



#single cell dataset proportion
# -----X02------- --
sc02.prop <- pheno02.matrix$manual_annotation
sc02.prop[grepl("tumor", sc02.prop)] <- "tumor"
sc02.prop <- prop.table(table(sc02.prop))


# -----X03------- --
sc03.prop <- pheno03.matrix$manual_annotation
sc03.prop[grepl("tumor", sc03.prop)] <- "tumor"
sc03.prop <- prop.table(table(sc03.prop))

# -----X04------- --
sc04.prop <- pheno04.matrix$manual_annotation
sc04.prop[grepl("tumor", sc04.prop)] <- "tumor"
sc04.prop <- prop.table(table(sc04.prop))



save(sc02.prop, sc03.prop, sc04.prop, result02, result03, result04, result0204, result0304, result0203, file = "TEDresults.RData")



# --------data visualization     X02--------------------------------------------
load("TEDresults.RData")



sc02.prop <- as.matrix(sc02.prop)
colnames(sc02.prop) <- c("actual")
rownames(sc02.prop)[rownames(sc02.prop) == "tumor"] <- "Tumor"


r02.x02 <- as.matrix(result02[2, ])
colnames(r02.x02) <- c("ref02")
r03.x02 <- as.matrix(result03[2, ])
colnames(r03.x02) <- c("ref03")
r04.x02 <- as.matrix(result04[2, ])
colnames(r04.x02) <- c("ref04")

r0203.x02 <- as.matrix(result0203[2, ])
colnames(r0203.x02) <- c("ref0203")
r0204.x02 <- as.matrix(result0204[2, ])
colnames(r0204.x02) <- c("ref0204")
r0304.x02 <- as.matrix(result0304[2, ])
colnames(r0304.x02) <- c("ref0304")


# combine into 1 dataframe
res <- Reduce(function(a,b){
  ans <- merge(a,b,by="row.names",all=T)
  row.names(ans) <- ans[,"Row.names"]
  ans[,!names(ans) %in% "Row.names"]
}, list(sc02.prop, r02.x02, r03.x02, r04.x02, r0203.x02, r0204.x02, r0304.x02))


library(tidyverse)
library(ggplot2)
library(reshape2)
res02 <- rownames_to_column(res, var = "rowname")

df1 <- melt(res02, "rowname")

g1 <- ggplot(df1, aes(x = rowname, y=value)) +
  geom_bar(aes(fill=variable),stat="identity", position ="dodge") + 
  theme_bw()+ 
  theme(axis.text.x = element_text(angle=-40, hjust=.1)) + 
  labs(x= "Cell Types", y = "Proportions", title = "Sample X02")

g2 <- ggplot(subset(df1, rowname %in% c("Tumor")), aes(x = rowname, y=value)) +
  geom_bar(aes(fill=variable),stat="identity", position ="dodge") + 
  theme_bw()+ 
  theme(axis.text.x = element_text(angle=-40, hjust=.1)) + 
  coord_cartesian(ylim=c(0.75, 1)) +
  labs(x= "Cell Types", y = "Proportions", title = "Sample X02 - Tumor Proportions")  
  



# --------data visualization     X03--------------------------------------------


sc03.prop <- as.matrix(sc03.prop)
colnames(sc03.prop) <- c("actual")
rownames(sc03.prop)[rownames(sc03.prop) == "tumor"] <- "Tumor"


r02.x03 <- as.matrix(result02[1, ])
colnames(r02.x03) <- c("ref02")
r03.x03 <- as.matrix(result03[1, ])
colnames(r03.x03) <- c("ref03")
r04.x03 <- as.matrix(result04[1, ])
colnames(r04.x03) <- c("ref04")

r0203.x03 <- as.matrix(result0203[1, ])
colnames(r0203.x03) <- c("ref0203")
r0204.x03 <- as.matrix(result0204[1, ])
colnames(r0204.x03) <- c("ref0204")
r0304.x03 <- as.matrix(result0304[1, ])
colnames(r0304.x03) <- c("ref0304")


# combine into 1 dataframe
res <- Reduce(function(a,b){
  ans <- merge(a,b,by="row.names",all=T)
  row.names(ans) <- ans[,"Row.names"]
  ans[,!names(ans) %in% "Row.names"]
}, list(sc03.prop, r02.x03, r03.x03, r04.x03, r0203.x03, r0204.x03, r0304.x03))


library(tidyverse)
library(ggplot2)
library(reshape2)
res03 <- rownames_to_column(res, var = "rowname")

df2 <- melt(res03, "rowname")

g1 <- ggplot(df2, aes(x = rowname, y=value)) +
  geom_bar(aes(fill=variable),stat="identity", position ="dodge") + 
  theme_bw()+ 
  theme(axis.text.x = element_text(angle=-40, hjust=.1)) + 
  labs(x= "Cell Types", y = "Proportions", title = "Sample X03")

g2 <- ggplot(subset(df2, rowname %in% c("Tumor")), aes(x = rowname, y=value)) +
  geom_bar(aes(fill=variable),stat="identity", position ="dodge") + 
  theme_bw()+ 
  theme(axis.text.x = element_text(angle=-40, hjust=.1)) + 
  coord_cartesian(ylim=c(0.3, 1)) +
  labs(x= "Cell Types", y = "Proportions", title = "Sample X03 - Tumor Proportions") 






# --------data visualization     X04--------------------------------------------



sc04.prop <- as.matrix(sc04.prop)
colnames(sc04.prop) <- c("actual")
rownames(sc04.prop)[rownames(sc04.prop) == "tumor"] <- "Tumor"


r02.x04 <- as.matrix(result02[3, ])
colnames(r02.x04) <- c("ref02")
r03.x04 <- as.matrix(result03[3, ])
colnames(r03.x04) <- c("ref03")
r04.x04 <- as.matrix(result04[3, ])
colnames(r04.x04) <- c("ref04")

r0203.x04 <- as.matrix(result0203[3, ])
colnames(r0203.x04) <- c("ref0203")
r0204.x04 <- as.matrix(result0204[3, ])
colnames(r0204.x04) <- c("ref0204")
r0304.x04 <- as.matrix(result0304[3, ])
colnames(r0304.x04) <- c("ref0304")


# combine into 1 dataframe
res <- Reduce(function(a,b){
  ans <- merge(a,b,by="row.names",all=T)
  row.names(ans) <- ans[,"Row.names"]
  ans[,!names(ans) %in% "Row.names"]
}, list(sc02.prop, r02.x04, r03.x04, r04.x04, r0203.x04, r0204.x04, r0304.x04))


library(tidyverse)
library(ggplot2)
library(reshape2)
res03 <- rownames_to_column(res, var = "rowname")

df3 <- melt(res03, "rowname")

g1 <- ggplot(df3, aes(x = rowname, y=value)) +
  geom_bar(aes(fill=variable),stat="identity", position ="dodge") + 
  theme_bw()+ 
  theme(axis.text.x = element_text(angle=-40, hjust=.1)) + 
  labs(x= "Cell Types", y = "Proportions", title = "Sample X04")

g2 <- ggplot(subset(df3, rowname %in% c("Tumor")), aes(x = rowname, y=value)) +
  geom_bar(aes(fill=variable),stat="identity", position ="dodge") + 
  theme_bw()+ 
  theme(axis.text.x = element_text(angle=-40, hjust=.1)) + 
  coord_cartesian(ylim=c(0.5, 1)) +
  labs(x= "Cell Types", y = "Proportions", title = "Sample X04 - Tumor Proportions")  













#-------------------Compare three trials for X02----------------------
sc02.prop <- as.matrix(sc02.prop)
colnames(sc02.prop) <- c("actual")
rownames(sc02.prop)[rownames(sc02.prop) == "tumor"] <- "Tumor"



r0304.x02 <- as.matrix(result0304[2, ])
colnames(r0304_1.x02) <- c("run 1")

r0304.x02 <- as.matrix(result0304[2, ])
colnames(r0304_1.x02) <- c("run 1")
r0203_3.x02 <- as.matrix(result0203_3[2, ])
colnames(r0203_3.x02) <- c("run 3")


# combine into 1 dataframe
res <- Reduce(function(a,b){
  ans <- merge(a,b,by="row.names",all=T)
  row.names(ans) <- ans[,"Row.names"]
  ans[,!names(ans) %in% "Row.names"]
}, list(sc02.prop, r0203.x02, r0203_2.x02, r0203_3.x02))



library(tidyverse)
library(ggplot2)
library(reshape2)
res02 <- rownames_to_column(res, var = "rowname")

df1 <- melt(res02, "rowname")

g1 <- ggplot(df1, aes(x = rowname, y=value)) +
  geom_bar(aes(fill=variable),stat="identity", position ="dodge") + 
  theme_bw()+ 
  theme(axis.text.x = element_text(angle=-40, hjust=.1)) + 
  labs(x= "Cell Types", y = "Proportions", title = "Sample X02")

g2 <- ggplot(subset(df1, rowname %in% c("Tumor")), aes(x = rowname, y=value)) +
  geom_bar(aes(fill=variable),stat="identity", position ="dodge") + 
  theme_bw()+ 
  theme(axis.text.x = element_text(angle=-40, hjust=.1)) + 
  coord_cartesian(ylim=c(0.75, 1)) +
  labs(x= "Cell Types", y = "Proportions", title = "Sample X02 - Tumor Proportions")  




#-------------------Compare three trials for X03 ----------------------
sc03.prop <- as.matrix(sc03.prop)
colnames(sc03.prop) <- c("actual")
rownames(sc03.prop)[rownames(sc03.prop) == "tumor"] <- "Tumor"

r03.x03 <- as.matrix(result03[1, ])
colnames(r03.x03) <- c("ref03")

r0304.x03 <- as.matrix(result0304[1, ])
colnames(r0304.x03) <- c("0304-run 1")

r0304_2.x03 <- as.matrix(result0304_2[1, ])
colnames(r0304_2.x03) <- c("0304-run 2")

r0304_3.x03 <- as.matrix(result0304_3[1, ])
colnames(r0304_3.x03) <- c("0304-run 3")


# combine into 1 dataframe
res <- Reduce(function(a,b){
  ans <- merge(a,b,by="row.names",all=T)
  row.names(ans) <- ans[,"Row.names"]
  ans[,!names(ans) %in% "Row.names"]
}, list(sc03.prop, r03.x03, r0304.x03, r0304_2.x03, r0304_3.x03))



library(tidyverse)
library(ggplot2)
library(reshape2)
res02 <- rownames_to_column(res, var = "rowname")

df1 <- melt(res02, "rowname")

g1 <- ggplot(df1, aes(x = rowname, y=value)) +
  geom_bar(aes(fill=variable),stat="identity", position ="dodge") + 
  theme_bw()+ 
  theme(axis.text.x = element_text(angle=-40, hjust=.1)) + 
  labs(x= "Cell Types", y = "Proportions", title = "Sample X03 - 3 repeats of Ref0304")

g2 <- ggplot(subset(df1, rowname %in% c("Tumor")), aes(x = rowname, y=value)) +
  geom_bar(aes(fill=variable),stat="identity", position ="dodge") + 
  theme_bw()+ 
  theme(axis.text.x = element_text(angle=-40, hjust=.1)) + 
  coord_cartesian(ylim=c(0.75, 1)) +
  labs(x= "Cell Types", y = "Proportions", title = "Sample X03 - Tumor Proportions")  












r0204.x02 <- as.matrix(result0204[2, ])
colnames(r0204.x02) <- c("ref0204")
r0304.x02 <- as.matrix(result0304[2, ])
colnames(r0304.x02) <- c("ref0304")


# combine into 1 dataframe
res <- Reduce(function(a,b){
  ans <- merge(a,b,by="row.names",all=T)
  row.names(ans) <- ans[,"Row.names"]
  ans[,!names(ans) %in% "Row.names"]
}, list(sc02.prop, r02.x02, r03.x02, r04.x02, r0203.x02, r0204.x02, r0304.x02))


library(tidyverse)
library(ggplot2)
library(reshape2)
res02 <- rownames_to_column(res, var = "rowname")

df1 <- melt(res02, "rowname")

g1 <- ggplot(df1, aes(x = rowname, y=value)) +
  geom_bar(aes(fill=variable),stat="identity", position ="dodge") + 
  theme_bw()+ 
  theme(axis.text.x = element_text(angle=-40, hjust=.1)) + 
  labs(x= "Cell Types", y = "Proportions", title = "Sample X02")

g2 <- ggplot(subset(df1, rowname %in% c("Tumor")), aes(x = rowname, y=value)) +
  geom_bar(aes(fill=variable),stat="identity", position ="dodge") + 
  theme_bw()+ 
  theme(axis.text.x = element_text(angle=-40, hjust=.1)) + 
  coord_cartesian(ylim=c(0.75, 1)) +
  labs(x= "Cell Types", y = "Proportions", title = "Sample X02 - Tumor Proportions")  








