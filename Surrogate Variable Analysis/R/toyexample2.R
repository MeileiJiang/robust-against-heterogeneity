############################################################################################################
## toyexample2.R
## This is 2-d toy example. There are two classes and each class is a Gaussian mixture. The Gaussian
## mixture will be consider as heterogeneity. SVA will be applied to alliviate this problem. SVA is expected
## to be good at regression rather than classification.
## Author: Meilei
###########################################################################################################
library(dplyr)
library(ggplot2)
library(reshape2)
# Each class have 40 samples and each sample has 100 features.
# Feature 1 is constant signal. For Class 1, X1 = 4; for Class 2, X1 = -4.
# Feature 2 is batch signal. For Batch 1, X2 = 4; for Batch 2, X2 = -4
# Other features are random noise.

X1 = c(rep(4, 40), rep(-4, 40))
X2 = c(rep(c(4, -4, 4, -4), rep(20, 4)))
XMat = matrix(rnorm(98 * 80, 0, 1), ncol = 98)
colnames(XMat) = paste0('X', 3:100)
M = data.frame(X1, X2, XMat)

y = c(rep('A', 40), rep('B', 40))
batch = c(rep(c(1, 0, 1, 0), rep(20, 4)))

train = data.frame(y, M, batch)


# Visualize the training data ---------------------------------------------

Edata = t(M)
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()


# save the data set -------------------------------------------------------

save(train,file = "rdata/constToyexample.RData")
