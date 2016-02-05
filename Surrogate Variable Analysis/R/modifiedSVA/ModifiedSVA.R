###########################################################################
## ModifiedSVA.R
## This file is written to understand the behavior of ModifiedSVA
## Author: Meilei
###########################################################################
library(ggplot2)
library(dplyr)
library(reshape2)
library(sva)

source('R/mean_heterogeneity/fixed_effects_mean_hetero_sva/Generate_mean_hetero_data.R')
source('R/modifiedSVA/newirwsva.R')
source('R/modifiedSVA/helper.R')

# Privous research shows that IRW-SVA could be problematic at the case when the data has no genes only
# affected by batches.

n = 80; # Total number of samples
c1 = 40; # Total number of samples in the Class 1
c2 = 40; # Total number of samples in the Class 2
t = c(20, 20, 20, 40) # Settings of features

tempdata <- Generate_mean_hetero_data(pi_1 = 0.5, pi_2 = 0.5, t = t)
M <- tempdata %>% select(-y, -batch, -pi_1, -pi_2, -c1, -c2)
# Expression data
Edata <- t(M)
colnames(Edata) <- paste0("sample", c(1: dim(Edata)[2]))

# visualize the dataset
Mdata = melt(Edata)

gg0 = ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value", title = 'Original Simulation Data Set') +
  geom_tile() + 
  scale_fill_gradient2(limits=c(-8, 8)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

pdf(file = 'figures/Modified_SVA/simulate0.pdf', width = 8, height = 8)
print(gg0)
dev.off()

## full model 
mod <- model.matrix(~ y, data = tempdata)
## null model
mod0 <- model.matrix(~ 1, data = tempdata)
## consturct the hat matrix and get residual matrix R
HatB = Edata %*% mod %*% solve(t(mod) %*% mod) 
R = Edata - HatB %*% t(mod)

MR = melt(R)

ggplot(MR,aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value", 
       title = 'Sample data after removing the design matrix') +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

# get the PC1 from R
pc1 = svd(R)$v[,1]
plot(pc1)
## get the batch vector and class vector
bv = 2 - as.numeric(unlist(tempdata %>% select(batch)))
y = 2 - as.numeric(unlist(tempdata %>% select(y)))



## estimate the number of surrogate variable
n.sv <- num.sv(Edata, mod, method="be", B = 100)
## estimate the surrogate varialbe
# IRW-SVA
svobj <- sva(Edata, mod, mod0, n.sv= n.sv) 
# revised IRW-SVA
svobj2 <- newirwsva(Edata, mod, mod0, n.sv= n.sv) 

# Visualize the posterior probability of gamma and b

plot(svobj$pprob.gam)
plot(svobj$pprob.b)

plot(svobj2$pprob.gam)
plot(svobj2$pprob.b)

genes = data.frame(sva.pprob.gam = svobj$pprob.gam, sva.pprob.b = svobj$pprob.b, 
                   newsva.pprob.gam = svobj2$pprob.gam, newsva.pprob.b = svobj2$pprob.b, 
                   type = rep(c('Type A: b','Type C: b & gam', 'Type B: gam', 'Type D'), t), index = c(1:100))
rgenes = melt(genes, id.vars = c('type','index'))
rgenes$variable = factor(rgenes$variable, levels = c('sva.pprob.gam','newsva.pprob.gam', 'sva.pprob.b', 'newsva.pprob.b')) 
gg1 = ggplot(data = rgenes, aes(x = index, y = value, col = type)) +
  facet_wrap(~variable) +
  geom_point() +
  labs(x= 'Gene', y='probablity', col = 'Gene Type',
       title = 'Posteriori probability for gamma and b') +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())


pdf(file = 'figures/Modified_SVA/pprop4.pdf', width = 8, height = 8)
print(gg1)
dev.off()
## get the surrogate variable
sv <- svobj$sv
sv2 <- svobj2$sv

# vectors = data.frame('batch' = bv, 'pc1' = pc1, 'sv' = sv, 'new-sv' = sv2, 
#                      sample = c(1:80) )

vectors = data.frame('batch' = bv, 'sv' = sv, 'new-sv' = sv2, 
                     sample = c(1:80) )

mvectors = melt(vectors, id.vars = c('batch','sample'))

gg2 = ggplot(data = mvectors, aes(x = sample, y = value, col = as.factor(batch)) ) +
  facet_wrap(~variable) +
  geom_point() +
  labs(x= 'Sample', y='Value', col = 'Batch',
       title = 'Comparison of Surrogate Variables') +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())


pdf(file = 'figures/Modified_SVA/vector4.pdf', width = 8, height = 8)
print(gg2)
dev.off()  

## Remove the batch effect
fsvobj = fsva(Edata, mod, svobj, Edata)

fsvobj2 = fsva(Edata, mod, svobj2, Edata)

Edb = fsvobj$db
Edb2 = fsvobj2$db

Mdb = melt(Edb)
Mdb2 = melt(Edb2)

gg3 = ggplot(Mdb, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",title = "Remove Batch Effect through IRW-SVA") +
  geom_tile() + 
  scale_fill_gradient2(limits=c(-8, 8)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

gg4 = ggplot(Mdb2, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",title = "Remove Batch Effect through modified SVA") +
  geom_tile() + 
  scale_fill_gradient2(limits=c(-8, 8)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

pdf(file = 'figures/Modified_SVA/sva0.pdf', width = 8, height = 8)
print(gg3)
dev.off()

pdf(file = 'figures/Modified_SVA/new_sva0.pdf', width = 8, height = 8)
print(gg4)
dev.off()
