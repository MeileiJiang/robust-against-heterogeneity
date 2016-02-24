#############################################################################################
## 3-simplex.R
## Standard normal random vectors of size n = 3 tend to be a 3-simplex. 
## Verify the idea proposed by Hall, Marron and Neeman (2005) about geometry representation.
## Author: Meilei
############################################################################################
library(ggplot2)
library(dplyr)


# generate 3 normal random vectors ----------------------------------------

v1 = rnorm(1000)
v2 = rnorm(1000)
v3 = rnorm(1000)

sum(v1^2)
sum(v2^2)
sum(v3^2)

sum((v1-v2)^2)
