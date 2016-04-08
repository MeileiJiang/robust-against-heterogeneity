#######################################################################################################
## Grun_reads_cutoff.R
## Author: Meilei
#######################################################################################################
library(ggplot2)
library(dplyr)
library(reshape2)
library(RUnit)
library(sva)

source("R/helper.R")
source("R/nullirwsva.R")
load(file = "rdata/Grun_reads_cutoff.RData")


# visualize the dataset ---------------------------------------------------

# Expression

dim(Expression)
# [1]  251 2795

Expression1 = data.frame(Expression)
Expression2 = t(Expression1)
dim(Expression2)

mExpression =  melt(Expression2)
colnames(mExpression) = c("gene","cell","expression")

ge = ggplot(data = mExpression, aes(x = cell, y = gene, fill = expression)) +
  labs(x = "Cell", y = "Gene", fill = "Value", title = "Heatmap of Expression value for Grun reads cutoff") +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

ge

cellMean = colMeans(Expression2)
cellRank = names(sort(cellMean)) 

sortExpression2 = Expression2[, cellRank]
msortExpression2 = melt(sortExpression2)
colnames(msortExpression2) = c("gene","cell","expression")

gse = ggplot(data = msortExpression2, aes(x = cell, y = gene, fill = expression)) +
  labs(x = "Cell", y = "Gene", fill = "Value", title = "Heatmap of Expression value for Grun reads cutoff \n sorted by column mean") +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

gse




# Spikein

dim(Spikein)
# [1] 251  17

Spikein1 = data.frame(Spikein) 
checkEquals(rownames(Spikein1), rownames(Expression1)) 

Spikein2 = t(Spikein1)
mSpikein = melt(Spikein2)
colnames(mSpikein) = c("spikein","cell","expression")

gs = ggplot(data = mSpikein, aes(x = cell, y = spikein, fill = expression)) +
  labs(x = "Cell", y = "Spikein", fill = "Value", title = "Heatmap of Spikein for Grun reads cutoff") +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

gs

sortSpikein2 = Spikein2[, cellRank]
msortSpikein2 = melt(sortSpikein2)
colnames(msortSpikein2) = c("gene","cell","expression")

gss = ggplot(data = msortSpikein2, aes(x = cell, y = gene, fill = expression)) +
  labs(x = "Cell", y = "Gene", fill = "Value", title = "Heatmap of  Spikein for Grun reads cutoff \n sorted by column mean") +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

gss

# Surrogate Variable Analysis  --------------------------------------------

# Expression Matrix

exp.mod = model.matrix(~1, data = Expression1)
exp.n.sv = num.sv(Expression2, mod = exp.mod, B = 50)
exp.sv.obj = nullirwsva(dat = Expression2, n.sv = exp.n.sv)
table(exp.sv.obj$pprob.gam)
#    1 
# 2795 
exp.sv = exp.sv.obj$sv
Id <- diag(dim(Expression2)[2])
exp.resid <- Expression2 %*% (Id - exp.sv %*% solve(t(exp.sv) %*% exp.sv) %*% 
                    t(exp.sv))

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

# Expression Matrix

exp.svd <- svd(Expression2)
exp.d = data.frame(eigen = exp.svd$d, index = c(1:dim(Expression2)[2]))
exp.pc <- exp.svd$v[,1:(exp.n.sv+1)]
pairs(exp.pc)
exp.pc.df = data.frame(exp.pc)
colnames(exp.pc.df) = paste0("PC", 1:(exp.n.sv+1))

ge.scr = ggplot(data = exp.d, aes(x = index, y = eigen) ) +
  geom_point(col = "blue") +
  geom_line(col = "red") +
  labs(x = "Principal Components", y = "Eigenvalue", 
       title = "Scree plot of Expression Matrix")

exp.resid.svd <- svd(exp.resid)
exp.resid.d <- data.frame(eigen = exp.resid.svd$d, index = c(1:length(exp.resid.svd$d)))
exp.resid.pc <- exp.resid.svd$v[,1:5]
pairs(exp.resid.pc)

# Spikein Matrix

spi.svd <- svd(Spikein2)
spi.d = data.frame(eigen = spi.svd$d, index = c(1:dim(Spikein2)[1]))
spi.pc = spi.svd$v[,1:(spi.n.sv+1)]
pairs(spi.pc)

gs.scr = ggplot(data = spi.d, aes(x = index, y = eigen) ) +
  geom_point(col = "blue") +
  geom_line(col = "red") +
  labs(x = "Principal Components", y = "Eigenvalue", 
       title = "Scree plot of Expression Matrix")



# summarize the result ----------------------------------------------------

pdf("figures/Grun_reads_cutoff.pdf")
print(ge)
print(gs)
print(ge.scr)
print(gs.scr)
pairs(exp.pc, main = "Scatterplot of Principal Components For Expression Matrix")
pairs(spi.pc, main = "Scatterplot of Principal Components For Spikein Matrix")
dev.off()

