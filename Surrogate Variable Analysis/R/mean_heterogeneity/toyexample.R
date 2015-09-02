############################################################################################################
## toyexample.R
## This is 2-d toy example. There are two classes and each class is a Gaussian mixture. The Gaussian
## mixture will be consider as heterogeneity. SVA will be applied to alliviate this problem. SVA is expected
## to be good at regression rather than classification.
## Author: Meilei
###########################################################################################################

# Class1: X ~ pi_1 N([-4, 4, rep(0, 98)], I_100 ) + (1 - pi_1) N([4, 4, rep(0, 98)], I_100)
# Class2: X ~ pi_2 N([-4, -4, rep(0, 98)], I_100 ) + (1 - pi_2) N([4, -4, rep(0, 98)], I_100)

library(dplyr)
library(ggplot2)
library(mnormt)
library(reshape2)
# function to generate Gaussion mixture -----------------------------------

toy = function(n, pi, mu1, mu2, cov1, cov2){
  m = length(mu1)
  x = matrix(rep(NA, m * n), ncol = m)
  k = rbinom(n, 1, pi)
  for(i in 1:n){
    x[i, ] = (1 - k[i]) * rmnorm(mean = mu1, varcov = cov1) + k[i] * rmnorm(mean = mu2, varcov = cov2)
  }
  list(data = x, label = k)
}


# generate the data -------------------------------------------------------

# 3-d training set --------------------------------------------------------


n = 100
pi = 0.5
mu1 = c(4, -4, rep(0, 1))
mu2 = c(4, 4, rep(0, 1))
mu3 = c(-4, -4, rep(0, 1))
mu4 = c(-4, 4, rep(0, 1))
cov1 = diag(rep(1,3))
cov2 = diag(rep(1,3))
cov3 = diag(rep(1,3))
cov4 = diag(rep(1,3))

M1 = toy(n, pi, mu1, mu2, cov1, cov2)
M2 = toy(n, pi, mu3, mu4, cov3, cov4)
M = rbind(M1$data, M2$data)
label = c(M1$label, M2$label)
train3d = data.frame(y = c(rep("A",n), rep("B", n)), M, batch = as.factor(label))


# 10-d training set -------------------------------------------------------

n = 100
pi = 0.5
mu1 = c(2, -2, rep(0, 8))
mu2 = c(2, 2, rep(0, 8))
mu3 = c(-2, -2, rep(0, 8))
mu4 = c(-2, 2, rep(0, 8))
cov1 = diag(rep(1,10))
cov2 = diag(rep(1,10))
cov3 = diag(rep(1,10))
cov4 = diag(rep(1,10))

M1 = toy(n, pi, mu1, mu2, cov1, cov2)
M2 = toy(n, pi, mu3, mu4, cov3, cov4)
M = rbind(M1$data, M2$data)
label = c(M1$label, M2$label)
train10d = data.frame(y = c(rep("A",n), rep("B", n)), M, batch = as.factor(label))


# 100-training set --------------------------------------------------------
n = 100
pi = 0.5
mu1 = c(2, -2, rep(0, 98))
mu2 = c(2, 2, rep(0, 98))
mu3 = c(-2, -2, rep(0, 98))
mu4 = c(-2, 2, rep(0, 98))
cov1 = diag(rep(1,100))
cov2 = diag(rep(1,100))
cov3 = diag(rep(1,100))
cov4 = diag(rep(1,100))

M1 = toy(n, pi, mu1, mu2, cov1, cov2)
M2 = toy(n, pi, mu3, mu4, cov3, cov4)
M = rbind(M1$data, M2$data)
label = c(M1$label, M2$label)
train100d = data.frame(y = c(rep("A",n), rep("B", n)), M, batch = as.factor(label))
# visualize data ----------------------------------------------------------

ggplot(data = train, aes(x = X1, y = X2, color = y)) + geom_point(aes(shape = label))
ggplot(data = train, aes(x = X1, y = X3, color = y)) + geom_point(aes(shape = label))
ggplot(data = train, aes(x = X2, y = X3, color = y)) + geom_point(aes(shape = label))
ggplot(data = train, aes(x = X3, y = X4, color = y)) + geom_point(aes(shape = label))

Edata = t(M)
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()

# save the simulation data ------------------------------------------------

save(train3d, train10d, train100d, file = "rdata/toyexample.RData")


