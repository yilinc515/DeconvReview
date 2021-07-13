load("resShuffled.RData")
dat <- readRDS("plot1data.RDS")

res.shuffled.bisque <- t(res.shuffled.bisque)
res.shuffled.music <- as.data.frame.matrix(res.shuffled.music)
res.shuffled.bayesprism <- as.data.frame.matrix(res.shuffled.bayesprism)
res.shuffled.bisque <- as.data.frame.matrix(res.shuffled.bisque)
res.shuffled.scdc <- as.data.frame.matrix(res.shuffled.scdc)
res.shuffled.cibersortx <- as.data.frame.matrix(res.shuffled.cibersortx)
