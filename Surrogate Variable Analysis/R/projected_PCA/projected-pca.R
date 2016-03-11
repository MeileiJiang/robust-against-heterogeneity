################################################################################
## Projected-PCA.R
## Applied Projected-PCA to adjust batch effects when we know the batch labels.
## Author: Meilei
## Date: March 2016
################################################################################
library(ggplot2)
library(dplyr)
library(reshape2)
library(sva)

source('R/mean_heterogeneity/fixed_effects_mean_hetero_sva/Generate_mean_hetero_data.R')
source('R/modifiedSVA/newirwsva.R')
source('R/modifiedSVA/helper.R')

data = Generate_mean_hetero_data(pi_1= 0.1, pi_2 = 0.9, c1 = 10, c2 = 10, t = c(0,0,5000,5000), noise = 10)

M <- data %>% select(-y, -batch, -pi_1, -pi_2, -c1, -c2)
Edata <- t(M)
colnames(Edata) <- paste0("sample", c(1: dim(Edata)[2]))

# # visualize the dataset
# Mdata = melt(Edata)
# 
# gg0 = ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
#   labs(x = "Sample", y = "Gene", fill = "Value", title = 'Original Simulation Data Set') +
#   geom_tile() + 
#   scale_fill_gradient2(limits=c(-8, 8)) +
#   theme(plot.title = element_text(size = 20, face = "bold"), axis.ticks = element_blank(), 
#         axis.title=element_text(size=16, face="bold"), axis.text.x = element_blank(),
#         legend.title = element_text(size=14, face="bold"), legend.text = element_text(face="bold"))
# print(gg0)

## full model 
mod <- model.matrix(~ y, data = data)
## null model
mod0 <- model.matrix(~ 1, data = data)
## consturct the hat matrix and get residual matrix R
HatB = Edata %*% mod %*% solve(t(mod) %*% mod) 
R = Edata - HatB %*% t(mod)

# get the PC1 from R
pc1 = svd(R)$v[,1]
plot(pc1)
## get the batch vector and class vector
bv = 2 - as.numeric(unlist(data %>% select(batch)))
y = 2 - as.numeric(unlist(data %>% select(y)))

cor(pc1, bv)

## estimate the number of surrogate variable
n.sv <- num.sv(Edata, mod, method="be", B = 100)
## estimate the surrogate varialbe
# IRW-SVA
svobj <- sva(Edata, mod, mod0, n.sv= n.sv) 
# revised IRW-SVA
svobj2 <- newirwsva(Edata, mod, mod0, n.sv= n.sv) 

plot(svobj$pprob.gam)
plot(svobj$pprob.b)

plot(svobj2$pprob.gam)
plot(svobj2$pprob.b)

## get the surrogate variable
sv <- svobj$sv
sv2 <- svobj2$sv
plot(sv)
plot(sv2)

cor(sv, bv)
cor(sv2, bv)
