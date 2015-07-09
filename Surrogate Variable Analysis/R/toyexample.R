############################################################################################################
## toyexample.R
## This is 2-d toy example. There are two classes and each class is a Gaussian mixture. The Gaussian
## mixture will be consider as heterogeneity. SVA will be applied to alliviate this problem. SVA is expected
## to be good at regression rather than classification.
## Author: Meilei
###########################################################################################################

# Class1: X ~ pi_1 N([-4, 4], I_2 ) + (1 - pi_1) N([4, 4], I_2)
# Class2: X ~ pi_2 N([-4, -4], I_2 ) + (1 - pi_2) N([4, -4], I_2)

library(dplyr)
library(ggplot2)
library(mnormt)

# function to generate Gaussion mixture -----------------------------------

toy = function(n, pi, mu1, mu2, cov1, cov2){
  x = matrix(rep(NA, 2 * n), ncol = 2)
  k = rbinom(n, 1, pi)
  for(i in 1:n){
    x[i, ] = (1 - k[i]) * rmnorm(mean = mu1, varcov = cov1) + k[i] * rmnorm(mean = mu2, varcov = cov2)
  }
  list(data = x, label = k)
}


# generate the data -------------------------------------------------------


n = 200
pi = 0.5
mu1 = c(-4, 4)
mu2 = c(4, 4)
mu3 = c(-4, -4)
mu4 = c(4, -4)
cov1 = diag(c(1,1))
cov2 = diag(c(1,1))
cov3 = diag(c(1,1))
cov4 = diag(c(1,1))

M1 = toy(n, pi, mu1, mu2, cov1, cov2)
M2 = toy(n, pi, mu3, mu4, cov3, cov4)
M = rbind(M1$data, M2$data)
label = c(M1$label, M2$label)
train = data.frame(y = c(rep(1,200), rep(-1, 200)), M, label = label)


# visualize data ----------------------------------------------------------

ggplot(data = train, aes(x = X1, y = X2, color = factor(y))) + geom_point(aes(shape = factor(label)))


# save the simulation data ------------------------------------------------

save(train, file = "rdata/toyexample1.RData")


