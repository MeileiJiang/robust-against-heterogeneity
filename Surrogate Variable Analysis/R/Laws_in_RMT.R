###############################################################################
## Laws_in_RMT.R
## Famous laws in the Random Matrix Theory: 
## 1. Semicircle law. 
## 2. Marcenko-Pastur law. 
## Author: Meilei
###############################################################################
library(ggplot2)
library(dplyr)
library(mnormt)

A = matrix(rnorm(1000*1000), ncol = 1000)


# semicircle law ----------------------------------------------------------


X = (A + t(A))/(2*sqrt(1000))

eigen.X = eigen(X)
eg.x = data.frame(egv = eigen.X$values)
hist(eg.x$egv)
ggplot(data = eg.x, aes(x = egv)) + 
  geom_density(col = "red") + 
  geom_histogram( aes(y = ..density..))

