###########################################################################
## ModifiedSVA_Case1.R
## This file is written to understand the behavior of ModifiedSVA
## Case 1: The data set does not contains genes with batch label not with class label signal.
## Author: Meilei
###########################################################################
library(ggplot2)
library(dplyr)
library(reshape2)
library(sva)

source('R/modifiedSVA/Generate_mean_hetero_data_new.R')
source('R/modifiedSVA/newirwsva.R')
source('R/modifiedSVA/helper.R')

# Privous research shows that IRW-SVA could be problematic at the case when the data has no genes only
# affected by batches.

n = 80; # Total number of samples
c1 = 40; # Total number of samples in the Class 1
c2 = 40; # Total number of samples in the Class 2
t = c(0, 0, 50, 50) # Settings of features

tempdata <- Generate_mean_hetero_data(pi_1 = 0.5, pi_2 = 0.5, t = t)
M <- tempdata %>% select(-y, -batch, -pi_1, -pi_2, -c1, -c2)

# save simulation data in Case 1 ------------------------------------------

save(tempdata, file = "rdata/simulate_data_case_1.RData")
# Expression data
Edata <- t(M)
colnames(Edata) <- paste0("sample", c(1: dim(Edata)[2]))

# visualize the dataset
Mdata = melt(Edata)

gg0 = ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value", title = 'Original Simulation Data Set') +
  geom_tile() + 
  scale_fill_gradient2(limits=c(-8, 8)) +
  theme(plot.title = element_text(size = 20, face = "bold"), axis.ticks = element_blank(), 
        axis.title=element_text(size=16, face="bold"), axis.text.x = element_blank(),
        legend.title = element_text(size=14, face="bold"), legend.text = element_text(face="bold"))


pdf(file = 'figures/Modified_SVA/simulate1.pdf', width = 8, height = 8)
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

genes.pprob.gam = data.frame(irw_sva.pprob.gam = svobj$pprob.gam, new_sva.pprob.gam = svobj2$pprob.gam, 
                             type = rep(c('Type A: b', 'Type B: gam', 'Type C: b & gam', 'Type D: noise'), t), index = c(1:100))
genes.pprob.b = data.frame(irw_sva.pprob.b = svobj$pprob.b, new_sva.pprob.b = svobj2$pprob.b, 
                           type = rep(c('Type A: b', 'Type B: gam', 'Type C: b & gam', 'Type D: noise'), t), index = c(1:100))
mgenes.pprob.gam = melt(genes.pprob.gam, id.vars = c('type','index'))
mgenes.pprob.gam$variable = factor(mgenes.pprob.gam$variable, levels = c('irw_sva.pprob.gam', 'new_sva.pprob.gam')) 
mgenes.pprob.b = melt(genes.pprob.b, id.vars = c('type','index'))
mgenes.pprob.b$variable = factor(mgenes.pprob.b$variable, levels = c('irw_sva.pprob.b', 'new_sva.pprob.b')) 

gg1.1 = ggplot(data = mgenes.pprob.gam, aes(x = index, y = value, col = type)) +
  facet_wrap(~variable) +
  geom_point() +
  labs(x= 'Gene', y='probablity', col = 'Gene Type:',
       title = expression(paste('Posterior probability of ', gamma != 0, ' for each gene'))) +
  theme(plot.title = element_text(size = 20, face = "bold"), axis.ticks = element_blank(), 
        axis.title=element_text(size=16, face="bold"), axis.text.x = element_blank(), 
        strip.text.x = element_text(size = 16),
        legend.title = element_text(size=14, face="bold"), legend.text = element_text(size=12, face="bold"),
        legend.position = 'bottom')


pdf(file = 'figures/Modified_SVA/pprop1_1.pdf', width = 7, height = 6)
print(gg1.1)
dev.off()

gg1.2 = ggplot(data = mgenes.pprob.b, aes(x = index, y = value, col = type)) +
  facet_wrap(~variable) +
  geom_point() +
  labs(x= 'Gene', y='probablity', col = 'Gene Type:',
       title = expression(paste('Posterior probability of ', b != 0, ' for each gene'))) +
  theme(plot.title = element_text(size = 20, face = "bold"), axis.ticks = element_blank(), 
        axis.title=element_text(size=16, face="bold"), axis.text.x = element_blank(), 
        strip.text.x = element_text(size = 16),
        legend.title = element_text(size=14, face="bold"), legend.text = element_text(size=12, face="bold"),
        legend.position = 'bottom')



pdf(file = 'figures/Modified_SVA/pprop1_2.pdf', width = 7, height = 6)
print(gg1.2)
dev.off()
## get the surrogate variable
sv <- svobj$sv
sv2 <- svobj2$sv

# vectors = data.frame('batch' = bv, 'pc1' = pc1, 'sv' = sv, 'new-sv' = sv2, 
#                      sample = c(1:80) )

vectors = data.frame('batch' = bv, 'sv' = sv, 'new_sv' = sv2, 
                     sample = c(1:80) )

mvectors = melt(vectors, id.vars = c('batch','sample'))

gg2 = ggplot(data = mvectors, aes(x = sample, y = value, col = as.factor(batch)) ) +
  facet_wrap(~variable) +
  geom_point() +
  labs(x= 'Sample', y='Value', col = 'Batch:',
       title = 'Comparison of Surrogate Variables') +
  scale_colour_manual(values = c('orange','blue')) +
  theme(plot.title = element_text(size = 20, face = "bold"), axis.ticks = element_blank(), 
        axis.title=element_text(size=16, face="bold"), axis.text.x = element_blank(), 
        strip.text.x = element_text(size = 16),
        legend.title = element_text(size=14, face="bold"), legend.text = element_text(size=12, face="bold"),
        legend.position = 'bottom')


pdf(file = 'figures/Modified_SVA/vector1.pdf', width = 7, height = 6)
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
  theme(plot.title = element_text(size = 20, face = "bold"), axis.ticks = element_blank(), 
        axis.title=element_text(size=16, face="bold"), axis.text.x = element_blank(),
        legend.title = element_text(size=14, face="bold"), legend.text = element_text(face="bold"))

gg4 = ggplot(Mdb2, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",title = "Remove Batch Effect through modified SVA") +
  geom_tile() + 
  scale_fill_gradient2(limits=c(-8, 8)) +
  theme(plot.title = element_text(size = 20, face = "bold"), axis.ticks = element_blank(), 
        axis.title=element_text(size=16, face="bold"), axis.text.x = element_blank(),
        legend.title = element_text(size=14, face="bold"), legend.text = element_text(face="bold"))


pdf(file = 'figures/Modified_SVA/sva1.pdf', width = 8, height = 8)
print(gg3)
dev.off()

pdf(file = 'figures/Modified_SVA/new_sva1.pdf', width = 8, height = 8)
print(gg4)
dev.off()
