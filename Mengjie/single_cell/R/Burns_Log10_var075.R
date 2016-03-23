#######################################################################################################
## Burns_Log10_var075.R
## Author: Meilei
#######################################################################################################
library(dplyr)
library(ggplot2)
library(reshape2)
library(sva)

load("rdata/Burns_Log10_var075.RData")
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

mod = model.matrix(~1, data = Expression2)

# Visualize ---------------------------------------------------------------

mControl2 = melt(Control2)
colnames(mControl2) = c("gene","cell","expression")
ggplot(data = mControl2, aes(x = cell, y = gene, fill = expression)) +
  labs(x = "Cell", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

mExpression2 = melt(Expression2)
colnames(mExpression2) = c("gene","cell","expression")
ggplot(data = mExpression2, aes(x = cell, y = gene, fill = expression)) +
  labs(x = "Cell", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())


# surrogate variable analysis ---------------------------------------------

mod = model.matrix(~1, data = Control1)
n.sv1 = num.sv(Control2, mod = mod, method="be")
n.sv2 = num.sv(Expression2, mod = mod, method="be")

svobj1 = sva(Control2, mod = mod, mod0 = mod, n.sv = n.sv1)
svobj2 = sva(Expression2, mod = mod, mod0 = mod, n.sv = n.sv2)



# pca ---------------------------------------------------------------------

svdobj1 = svd(Control2)
plot(svdobj1$d, type = "l")

pcs1 = svdobj1$v[,1:5]
pairs(pcs1)

pca1 = getPcaResult(Control2, varNames = colnames(Control2), scale=F, center = F)
pcaScreePlot(pcaResult = pca1, title="Scree Plot")
svdobj2 = svd(Expression2)
plot(svdobj2$d, type = "l")

pcs2 = svdobj2$v[,1:10]
pairs(pcs2)
