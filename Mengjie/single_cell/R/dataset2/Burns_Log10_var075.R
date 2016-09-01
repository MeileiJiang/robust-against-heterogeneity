#######################################################################################################
## Burns_Log10_var075.R
## Author: Meilei
#######################################################################################################
library(dplyr)
library(ggplot2)
library(reshape2)
library(sva)

load("rdata/Burns_Log10_var075.RData")
source("R/helper.R")
source("R/nullirwsva.R")
source('~/researchspace/robust-against-heterogeneity/Surrogate Variable Analysis/pcafuns/R/getPcaResult.R')
source('~/researchspace/robust-against-heterogeneity/Surrogate Variable Analysis/pcafuns/R/pcaScreePlot.R')

dim(Control)
# [1] 307 282
dim(Expression)
# [1]  307 3846
Control1 = data.frame(Control)
Control2 = t(Control1)
Expression1 = data.frame(Expression)
Expression2 = t(Expression1)


con.mod = model.matrix(~1, data = data.frame(Control2))
ep.mod = model.matrix(~1, data = data.frame(Expression2))

# Visualize ---------------------------------------------------------------

mControl2 = melt(Control2)
colnames(mControl2) = c("gene","cell","expression")
gc = ggplot(data = mControl2, aes(x = cell, y = gene, fill = expression)) +
  labs(x = "Cell", y = "Gene", fill = "Value",
       title = "Heatmap of Expression Value for Cell Cycle Genes In Burns Log10 var075 Dataset") +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

mExpression2 = melt(Expression2)
colnames(mExpression2) = c("gene","cell","expression")
ge = ggplot(data = mExpression2, aes(x = cell, y = gene, fill = expression)) +
  labs(x = "Cell", y = "Gene", fill = "Value",
       title = "Heatmap of Expression Value for Non Cell Cycle Genes In Burns Log10 var075 Dataset") +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())


# surrogate variable analysis ---------------------------------------------

mod = model.matrix(~1, data = Control1)
con.n.sv = num.sv(Control2, mod = mod, method="be")
exp.n.sv = num.sv(Expression2, mod = mod, method="be")

con.svobj = nullirwsva(dat = Control2, n.sv = con.n.sv)
exp.svobj = nullirwsva(dat = Expression2, n.sv = exp.n.sv)



# pca ---------------------------------------------------------------------

# Contral 
con.svd <- svd(Control2)
con.d = data.frame(eigen = con.svd$d, index = c(1:dim(Control2)[1]))
con.pc <- con.svd$v[,1:(con.n.sv+1)]
pairs(con.pc)
con.pc.df = data.frame(con.pc)
colnames(con.pc.df) = paste0("PC", 1:(con.n.sv+1))

gc.scr = ggplot(data = con.d, aes(x = index, y = eigen) ) +
  geom_point(col = "blue") +
  geom_line(col = "red") +
  labs(x = "Principal Components", y = "Eigenvalue",
       title = "Scree plot of Expression Value for Cell Cycle Genes In Burns Log10 var075 Dataset")

# Expression
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
        title = "Scree plot of Expression Value for Non Cell Cycle Genes In Burns Log10 var075 Dataset")


# summarize the result ----------------------------------------------------

pdf("figures/Brunr_Log10_var075.pdf")
print(gc)
print(ge)
print(gc.scr)
print(ge.scr)
pairs(con.pc, main = "Scatterplot of Principal Components For Cell Cycle Genes In Burns Log10 var075 Datasetx ")
pairs(exp.pc, main = "Scatterplot of Principal Components For Non Cell Cycle Genes In Burns Log10 var075 Dataset ")
dev.off()
