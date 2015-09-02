############################################################################################################
## unbalanced_scale_heterogeneity_example.R
## This is 2-d toy example. There are two classes and each class is a Gaussian mixture. The Gaussian
## mixture will be consider as heterogeneity. SVA will be applied to alliviate this problem. SVA is expected
## to be good at regression rather than classification. Especially, the batch effect is unbalanced.
## Now the batch effect is unbalanced.
## Author: Meilei
###########################################################################################################
library(dplyr)
library(ggplot2)
library(reshape2)

# first set of unbalanced scale heterogeity -------------------------------

# Each class have 40 samples and each sample has 100 features.
# Feature 1 is constant signal. For Class 1, X1 = 4; for Class 2, X1 = -4.
# Feature 2 is batch effect. For Batch 1, X2 ~ N(0, 0.01); for Batch 2, X2 ~ N(0,2).70% of Class 1 from Batch 1
# and 20 % of Class 2 from Batch 1. The others from Batch 2. 
# Other features are random noise.

X1 = c(rep(4, 40), rep(-4, 40))
X2 = c(rnorm(28, 0, 0.1), rnorm(12, 0, 2), rnorm(8, 0, 0.1), rnorm(32, 0, 2))
XMat = matrix(rnorm(98 * 80, 0, 1), ncol = 98)
colnames(XMat) = paste0('X', 3:100)
M = data.frame(X1, X2, XMat)

y = c(rep('Class1', 40), rep('Class2', 40))
batch = c(rep(c('Batch1', 'Batch2', 'Batch1', 'Batch2'), c(28, 12, 8, 32)))

train = data.frame(y, M, batch)


# Visualize the training data ---------------------------------------------

Edata = t(M)
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()

# Second set of unbalanced scale heterogeity -------------------------------

# Each class have 40 samples and each sample has 100 features.
# Feature 1 is constant signal. For Class 1, X1 = 4; for Class 2, X1 = -4.
# Feature 2 is batch effect. For Batch 1, X2 ~ N(0, 0.01); for Batch 2, X2 ~ N(0,2).33% of Class 1 from Batch 1
# and 67 % of Class 2 from Batch 1. The others from Batch 2. 
# Other features are random noise.

X1 = c(rep(4, 40), rep(-4, 40))
X2 = c(rnorm(13, 0, 0.1), rnorm(27, 0, 2), rnorm(27, 0, 0.1), rnorm(13, 0, 2))
XMat = matrix(rnorm(98 * 80, 0, 1), ncol = 98)
colnames(XMat) = paste0('X', 3:100)
M = data.frame(X1, X2, XMat)

y = c(rep('Class1', 40), rep('Class2', 40))
batch = c(rep(c('Batch1', 'Batch2', 'Batch1', 'Batch2'), c(13, 27, 27, 13)))

train2 = data.frame(y, M, batch)


# Visualize the training data ---------------------------------------------

Edata = t(M)
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()

# Third set of unbalanced scale heterogeity -------------------------------

# Each class have 40 samples and each sample has 100 features.
# Feature 1 is constant signal. For Class 1, X1 = 4; for Class 2, X1 = -4.
# Feature 2 is batch effect. For Batch 1, X2 ~ N(0, 0.01); for Batch 2, X2 ~ N(0,2).10% of Class 1 from Batch 1
# and 90% of Class 2 from Batch 1. The others from Batch 2. 
# Other features are random noise.

X1 = c(rep(4, 40), rep(-4, 40))
X2 = c(rnorm(4, 0, 0.1), rnorm(36, 0, 2), rnorm(36, 0, 0.1), rnorm(4, 0, 2))
XMat = matrix(rnorm(98 * 80, 0, 1), ncol = 98)
colnames(XMat) = paste0('X', 3:100)
M = data.frame(X1, X2, XMat)

y = c(rep('Class1', 40), rep('Class2', 40))
batch = c(rep(c('Batch1', 'Batch2', 'Batch1', 'Batch2'), c(4, 36, 36, 4)))

train3 = data.frame(y, M, batch)


# Visualize the training data ---------------------------------------------

Edata = t(M)
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()

# save the data set -------------------------------------------------------

save(train, train2, train3, file = "rdata/unbalanced_scale_hetero_example.RData")
