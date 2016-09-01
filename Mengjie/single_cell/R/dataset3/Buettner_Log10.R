#######################################################################################################
## Buettner_Log10.R
## Author: Meilei
#######################################################################################################
library(ggplot2)
library(dplyr)
library(reshape2)
library(RUnit)
library(sva)

source("R/helper.R")
source("R/nullirwsva.R")
load(file = "rdata/Buettner_Log10.RData")


# visualize the dataset ---------------------------------------------------

# Par1: non cell cycle gene

dim(part1)
# [1]  251 2795

part1.df = data.frame(t(part1))
dim(part1.df)

mpart1.df =  melt(part1.df)
colnames(mpart1.df ) = c("gene","cell","expression")

gp1 = ggplot(data = mpart1.df , aes(x = cell, y = gene, fill = expression)) +
  labs(x = "Cell", y = "Gene", fill = "Value", 
       title = "Heatmap of Expression value for Non Cell Cycle Gene in Buettner Log10 Data") +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

# Par2: cell cycle gene

dim(part2)
# [1] 182 629

part2.df = t(data.frame(part2))
dim(part2.df)

mpart2.df =  melt(part2.df)
colnames(mpart2.df ) = c("gene","cell","expression")

gp2 = ggplot(data = mpart1.df , aes(x = cell, y = gene, fill = expression)) +
  labs(x = "Cell", y = "Gene", fill = "Value", 
       title = "Heatmap of Expression value for Cell Cycle Gene in Buettner Log10 Data") +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

# Spikein

dim(Spikein)
# [1] 182  92

Spikein1 = data.frame(Spikein) 
checkEquals(rownames(Spikein1), rownames(part1)) 

Spikein2 = t(Spikein1)
mSpikein = melt(Spikein2)
colnames(mSpikein) = c("spikein","cell","expression")

gs = ggplot(data = mSpikein, aes(x = cell, y = spikein, fill = expression)) +
  labs(x = "Cell", y = "Spikein", fill = "Value", title = "Heatmap of Spikein for Buettner Log10 Data") +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())


# Surrogate Variable Analysis  --------------------------------------------

# Part1

p1.mod = model.matrix(~1, data = data.frame(t(part1.df)))
p1.n.sv = num.sv(part1.df, mod = p1.mod)
p1.sv.obj = nullirwsva(dat = part1.df, n.sv = p1.n.sv)
table(p1.sv.obj$pprob.gam)
#    1 
# 2795 
p1.sv = p1.sv.obj$sv
Id <- diag(dim(part1.df)[2])
p1.resid <- part1.df %*% (Id - p1.sv %*% solve(t(p1.sv) %*% p1.sv) %*% 
                                t(p1.sv))

# Part2

p2.mod = model.matrix(~1, data = data.frame(t(part2.df)))
p2.n.sv = num.sv(part2.df, mod = p2.mod)
p2.sv.obj = nullirwsva(dat = part2.df, n.sv = p2.n.sv)
table(p2.sv.obj$pprob.gam)
#    1 
# 2795 
p2.sv = p2.sv.obj$sv
Id <- diag(dim(part2.df)[2])
p2.resid <- part2.df %*% (Id - p2.sv %*% solve(t(p2.sv) %*% p2.sv) %*% 
                            t(p2.sv))

# Spikein Matrix

spi.mod = model.matrix(~1, data = Spikein1)
spi.n.sv = num.sv(Spikein2, mod = spi.mod, B = 50)
spi.sv.obj = nullirwsva(dat = Spikein2, n.sv = spi.n.sv)
table(spi.sv.obj$pprob.gam)
#  1 
# 17  
spi.sv = spi.sv.obj$sv
Id <- diag(dim(Spikein2)[2])
spi.resid <- Spikein2 %*% (Id - spi.sv %*% solve(t(spi.sv) %*% spi.sv) %*% 
                             t(spi.sv))


# Principal Component Analysis --------------------------------------------

# part1

p1.svd <- svd(part1.df)
p1.d = data.frame(eigen = p1.svd$d, index = c(1:dim(part1.df)[2]))
p1.pc <- p1.svd$v[,1:(p1.n.sv+1)]
pairs(p1.pc)
p1.pc.df = data.frame(p1.pc)
colnames(p1.pc.df) = paste0("PC", 1:(p1.n.sv+1))

gp1.scr = ggplot(data = p1.d, aes(x = index, y = eigen) ) +
  geom_point(col = "blue") +
  geom_line(col = "red") +
  labs(x = "Principal Components", y = "Eigenvalue", 
       title = "Scree plot of Expression Matrix for Non Cell Cycle Gene in Buettner Log10 Data")

p1.resid.svd <- svd(p1.resid)
p1.resid.d <- data.frame(eigen = p1.resid.svd$d, index = c(1:length(p1.resid.svd$d)))
p1.resid.pc <- p1.resid.svd$v[,1:5]
pairs(p1.resid.pc)

# part2

p2.svd <- svd(part2.df)
p2.d = data.frame(eigen = p2.svd$d, index = c(1:dim(part2.df)[2]))
p2.pc <- p2.svd$v[,1:(p2.n.sv+1)]
pairs(p2.pc)
p2.pc.df = data.frame(p2.pc)
colnames(p2.pc.df) = paste0("PC", 1:(p2.n.sv+1))

gp2.scr = ggplot(data = p2.d, aes(x = index, y = eigen) ) +
  geom_point(col = "blue") +
  geom_line(col = "red") +
  labs(x = "Principal Components", y = "Eigenvalue", 
       title = "Scree plot of Expression Matrix for Cell Cycle Gene in Buettner Log10 Data")

p2.resid.svd <- svd(p1.resid)
p2.resid.d <- data.frame(eigen = p2.resid.svd$d, index = c(1:length(p2.resid.svd$d)))
p2.resid.pc <- p2.resid.svd$v[,1:5]
pairs(p2.resid.pc)


# Spikein Matrix

spi.svd <- svd(Spikein2)
spi.d = data.frame(eigen = spi.svd$d, index = c(1:dim(Spikein2)[1]))
spi.pc = spi.svd$v[,1:(spi.n.sv+1)]
pairs(spi.pc)

gs.scr = ggplot(data = spi.d, aes(x = index, y = eigen) ) +
  geom_point(col = "blue") +
  geom_line(col = "red") +
  labs(x = "Principal Components", y = "Eigenvalue", 
       title = "Scree plot of Expression Matrix for Buettner Log10 Data")



# summarize the result ----------------------------------------------------

pdf("figures/Buetter_Log10.pdf")
print(gp1)
print(gp2)
print(gs)
print(gp1.scr)
print(gp2.scr)
print(gs.scr)
pairs(p1.pc, main = "Scatterplot of Principal Components For Expression Matrix for Non Cell Cycle Gene in Buettner Log10 Data")
pairs(p2.pc, main = "Scatterplot of Principal Components For Expression Matrix for Cell Cycle Gene in Buettner Log10 Data")
pairs(spi.pc, main = "Scatterplot of Principal Components For Spikein Matrix in Buettner Log10 Data")
dev.off()
