############################################################################################################
## unbalanced_mean_heterogeneity_example.R
## This is 2-d toy example. There are two classes and each class is a Gaussian mixture. The Gaussian
## mixture will be consider as heterogeneity. SVA will be applied to alliviate this problem. SVA is expected
## to be good at regression rather than classification.
## Moreover, our analysis shows that the detection of SV is better at the case of random effects than the 
## case of fixed effects. 
## Author: Meilei
###########################################################################################################

library(dplyr)
library(ggplot2)
library(reshape2)

# Each class have 40 samples and each sample has 100 features.
# Feature 1 is primary variable effects signal. For Class 1, X1 ~ N(4, 1); for Class 2, X1 ~ N(-4, 1).
# Feature 3 is batch effect. For Batch 1, X2 ~ N(2, 1); for Batch 2, X2 ~ N(-2, 1). pi_1 of Class 1 from Batch 1
# Feature 2 is affected by both primary variable and batch effect.
# and pi_2 % of Class 2 from Batch 1. The others from Batch 2. 
# Other features are random noise.
n = 80 # number of features
m = 97 # number of noise feature
# Balanced Case -----------------------------------------------------------

# pi_1 = pi_2 = 0.5
pi_1 = 0.5; pi_2 = 0.5;
n_11 = 40*pi_1; n_12 = 40 - n_11; n_21 = 40*pi_2; n_22 = 40 - n_21

X1 = c(rnorm(40, 2, 1), rnorm(40, -2, 1))
X3 = c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1))
X2 = X1 + X3 + rnorm(n, 0, 1)
#X2 = c(rep(2, n_11), rep(-2, n_12), rep(2, n_21), rep(-2, n_22))
XMat = matrix(rnorm(m * 80, 0, 1), ncol = m)
colnames(XMat) = paste0('X', 4:(m + 3))
M = data.frame(X1, X2, X3, XMat)

y = c(rep('Class1', 40), rep('Class2', 40))
batch = rep(c('Batch1', 'Batch2', 'Batch1', 'Batch2'), c(n_11, n_12, n_21, n_22) ) 


train0 = data.frame(y, M, batch, pi_1, pi_2)

# Visualize the training data 

Edata = t(M)
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()



# Unbalanced Case 1 -------------------------------------------------------

# In this case: pi_1 = pi_2, pi_1 + pi_2 < 1
# The propotion of each class in the Batch 1 are the same.
# The sample numbers of two batches are different.

## Case 1.1
pi_1 = 0.4; pi_2 = 0.4
  
n_11 = 40*pi_1; n_12 = 40 - n_11; n_21 = 40*pi_2; n_22 = 40 - n_21
X3 = c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1))
X2 = X1 + X3 + rnorm(n, 0, 1)
#X2 = c(rep(2, n_11), rep(-2, n_12), rep(2, n_21), rep(-2, n_22))
M = data.frame(X1, X2, X3, XMat)

batch = rep(c('Batch1', 'Batch2', 'Batch1', 'Batch2'), c(n_11, n_12, n_21, n_22) ) 
train1.1 = data.frame(y, M, batch, pi_1, pi_2)

## Visualize the training data 

Edata = t(M)
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()

## Case 1.2
pi_1 = 0.1; pi_2 = 0.1

n_11 = 40*pi_1; n_12 = 40 - n_11; n_21 = 40*pi_2; n_22 = 40 - n_21
X3 = c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1))
X2 = X1 + X3 + rnorm(n, 0, 1)
M = data.frame(X1, X2, X3, XMat)

batch = rep(c('Batch1', 'Batch2', 'Batch1', 'Batch2'), c(n_11, n_12, n_21, n_22) ) 
train1.2 = data.frame(y, M, batch, pi_1, pi_2)

## Visualize the training data 

Edata = t(M)
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()


# Unbalanced Case 2 -------------------------------------------------------

# pi_1 > pi_2, pi_1 + pi_2 = 1.
# The propotion of each class in the Batch 1 are different, 
# The sample numbers of two batches are the same.

## Case 2.1
pi_1 = 0.6; pi_2 = 0.4

n_11 = 40*pi_1; n_12 = 40 - n_11; n_21 = 40*pi_2; n_22 = 40 - n_21
X3 = c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1))
X2 = X1 + X3 + rnorm(n, 0, 1)
M = data.frame(X1, X2, X3, XMat)

batch = rep(c('Batch1', 'Batch2', 'Batch1', 'Batch2'), c(n_11, n_12, n_21, n_22) ) 
train2.1 = data.frame(y, M, batch, pi_1, pi_2)

## Visualize the training data 

Edata = t(M)
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()

## Case 2.2
pi_1 = 0.9; pi_2 = 0.1

n_11 = 40*pi_1; n_12 = 40 - n_11; n_21 = 40*pi_2; n_22 = 40 - n_21
X3 = c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1))
X2 = X1 + X3 + rnorm(n, 0, 1)
M = data.frame(X1, X2, X3, XMat)

batch = rep(c('Batch1', 'Batch2', 'Batch1', 'Batch2'), c(n_11, n_12, n_21, n_22) ) 
train2.2 = data.frame(y, M, batch, pi_1, pi_2)

## Visualize the training data 

Edata = t(M)
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()

# Unbalanced Case 3 -------------------------------------------------------

# pi_1 > pi_2, pi_1 + pi_2 < 1.
# The propotion of each class in the Batch 1 are different. 
# The sample numbers of two batches are different.

## Case 3.1
pi_1 = 0.55; pi_2 = 0.35

n_11 = 40*pi_1; n_12 = 40 - n_11; n_21 = 40*pi_2; n_22 = 40 - n_21
X3 = c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1))
X2 = X1 + X3 + rnorm(n, 0, 1)
M = data.frame(X1, X2, X3, XMat)

batch = rep(c('Batch1', 'Batch2', 'Batch1', 'Batch2'), c(n_11, n_12, n_21, n_22) ) 
train3.1 = data.frame(y, M, batch, pi_1, pi_2)

## Visualize the training data 

Edata = t(M)
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()

## Case 3.2
pi_1 = 0.4; pi_2 = 0.1

n_11 = 40*pi_1; n_12 = 40 - n_11; n_21 = 40*pi_2; n_22 = 40 - n_21
X3 = c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1))
X2 = X1 + X3 + rnorm(n, 0, 1)
M = data.frame(X1, X2, X3, XMat)

batch = rep(c('Batch1', 'Batch2', 'Batch1', 'Batch2'), c(n_11, n_12, n_21, n_22) ) 
train3.2 = data.frame(y, M, batch, pi_1, pi_2)

## Visualize the training data 

Edata = t(M)
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()


# None batch effect case --------------------------------------------------

# pi_1 = pi_2 = 1
# All samples come from the same batch

pi_1 = 1; pi_2 = 1

n_11 = 40*pi_1; n_12 = 40 - n_11; n_21 = 40*pi_2; n_22 = 40 - n_21
X3 = c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1))
X2 = X1 + X3 + rnorm(n, 0, 1)
M = data.frame(X1, X2, X3, XMat)

batch = rep(c('Batch1', 'Batch2', 'Batch1', 'Batch2'), c(n_11, n_12, n_21, n_22) ) 
train4 = data.frame(y, M, batch, pi_1, pi_2)

## Visualize the training data 

Edata = t(M)
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()

# save the data set -------------------------------------------------------

save(train0, train1.1, train1.2, train2.1, train2.2, 
     train3.1, train3.2, train4, file = "rdata/random_effects_mean_hetero_example.RData")
