#################################################################################
## unbalanced_mean_hetero_SV_detect_case2.R
## Try to understand the behaviour 
## Data has unbalanced 
##
## Meilei
#################################################################################
library(dplyr)
library(ggplot2)
library(sva)
library(reshape2)

# data processing ---------------------------------------------------------

load("rdata/mean_hetero_example.RData")
source('~/researchspace/robust-against-heterogeneity/Surrogate Variable Analysis/pcafuns/R/getPcaResult.R')

n = 80 # number of samples

# analysis first data set -------------------------------------------------

M1 <- train2.1 %>% select(-y, -batch, -pi_1, -pi_2)
pi_1 <- unique(train2.1$pi_1); pi_2 <- unique(train2.1$pi_2)
p_1 <- (pi_1 + pi_2)/2; p_2 <- 1 - p_1


# Expression data
Edata1 = t(M1)
colnames(Edata1) = paste0("sample", c(1: dim(Edata1)[2]))
Mdata1 = melt(Edata1)

title_ubg1.1 <- paste0("Simulation data \n pi_1 = ", pi_1, 
                       ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)
title_ubg2.1 <- paste0("Estimation of primary effects before adjustment \n pi_1 = ", pi_1, 
                       ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)
title_ubg3.1 <- paste0("Estimation of primary effects after adjustment \n pi_1 = ", pi_1, 
                       ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)
title_ubg4.1 <- paste0("Simulation data after ajustment \n pi_1 = ", pi_1, 
                       ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)
title_ubg5.1 <- paste0("Residual Matrix After Removing Primary Effect \n pi_1 = ", pi_1, 
                       ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)
title_ubg6.1 <- paste0("Eigenvalues Boxplot \n pi_1 = ", pi_1, 
                       ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)
title_ubg7.1 <- paste0("Screeplot \n pi_1 = ", pi_1, 
                       ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)
title_ubg8.1 <- paste0("Angles between PCs and Batch Effect \n blue dashed line: angle(SV, Batch Effect) \n pi_1 = ", pi_1, 
                       ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)

ubg1.1 = ggplot(Mdata1, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value", title = title_ubg1.1) +
  geom_tile() + 
  scale_fill_gradient2(limits=c(-10, 10)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())
print(ubg1.1)
# primary variable of interest
Y1 <- train2.1 %>% select(y)

# make model matrix -------------------------------------------------------
## full model 
mod1 = model.matrix(~ y, data = train2.1)
## null model
mod10 = model.matrix(~ 1, data = train2.1)


# estimate the number of latent factors that need to be estimated ---------
n.sv10 = num.sv(Edata1, mod1, method="leek")

n.sv1 = num.sv(Edata1, mod1, method="be", B = 1000)

print(paste0("The number of surrogate variable: ", n.sv1))

if(n.sv1 == 0){
  print("Mannually setting the number of surrogate variable as 1 to apply the sva algorithm")
  svobj1 = sva(t(M1),mod1, mod10,n.sv= 1, B = 5)
} else{
  svobj1 = sva(t(M1),mod1, mod10,n.sv= n.sv1, B = 5)
}


# Remove the batch effects
fsvobj1 = fsva(Edata1, mod1, svobj1, Edata1)

Edb1 = fsvobj1$db

# Estimate class effects on features from original data matrix

HatB1 = Edata1 %*% mod1 %*% solve(t(mod1) %*% mod1) 

Effect11 = HatB1 %*% matrix(c(1,1,0,1), nrow = 2, byrow = T)
rownames(Effect11) = rownames(Edata1)
colnames(Effect11) = c('yClass1', 'yClass2')
Effect11_df = data.frame(Effect11, feature = fe, x = factor(paste0('X',1:100), levels = paste0('X',1:100)))

MEff11_df = melt(Effect11_df)
ubg2.1 = ggplot(data = MEff11_df, aes(x = x, y = value, col = feature, shape = variable)) + 
  geom_point() + 
  scale_y_continuous(limits = c(-4.5, 4.5), breaks=-4:4) + 
  scale_shape_manual(values = c(4, 16)) + 
  labs(x = 'Gene',title = title_ubg2.1) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())
print(ubg2.1)
# Estimate class effects on features from original data matrix after removing the effects of surrogate variable 

HatB12 = Edb1 %*% mod1 %*% solve(t(mod1) %*% mod1) 

Effect12 = HatB12 %*% matrix(c(1,1,0,1), nrow = 2, byrow = T)
rownames(Effect12) = rownames(Edb1)
colnames(Effect12) = c('yClass1', 'yClass2')
Effect12_df = data.frame(Effect12, feature = fe, x = factor(paste0('X',1:100), levels = paste0('X',1:100)))

MEff12_df = melt(Effect12_df)
ubg3.1 = ggplot(data = MEff12_df, aes(x = x, y = value, col = feature, shape = variable)) + 
  geom_point() + 
  scale_y_continuous(limits = c(-4.5, 4.5), breaks=-4:4) + 
  scale_shape_manual(values = c(4, 16)) +
  labs(x = 'Gene',title = title_ubg3.1) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())
print(ubg3.1)

# visualize the result

Mdb1 = melt(Edb1)

ubg4.1 = ggplot(Mdb1, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",title = title_ubg4.1) +
  geom_tile() + 
  scale_fill_gradient2(limits=c(-10, 10)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

print(ubg4.1)
# dig into the analysis ---------------------------------------------------

R1 = Edata1 - HatB1 %*% t(mod1)
MR1 = melt(R1)
ubg5.1 =  ggplot(MR1, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",title = title_ubg5.1) +
  geom_tile() + 
  scale_fill_gradient2(limits=c(-10, 10)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())
print(ubg5.1)
R.pc1 = getPcaResult(R1, varNames = colnames(R1), scale=F, center = F)

R.pv1 =  data.frame(Var1 = "Origin", Var2 = rownames(R.pc1$varDf), value = R.pc1$varDf[,2])

# use permutation test to find the unusual large eigenvalues --------------

B = 1000
n1 = dim(R.pc1$varDf)[1]; n2 = dim(R1)[2]
pvMat = matrix(nrow = B, ncol = n1)
for(k in 1:B){
  # make permutation of each row independently
  newE = t(apply(R1, 1, sample, replace = FALSE))
  #  newE = apply(newE0, 2, sample, replace = FALSE)
  # refit the model to get the new residual matrix
  newR = newE - newE %*% mod1 %*% solve(t(mod1) %*% mod1) %*% t(mod1)
  colnames(newR) = c(1:n1)
  # do pca on newR and take out the proprotion variance vector
  newR.pc = getPcaResult(newR, varNames = colnames(newR), scale=F, center = F)
  pvMat[k, 1:n1] = newR.pc$varDf[,2]
}

colnames(pvMat) = rownames(R.pc1$varDf)
rownames(pvMat) = paste0("Run", 1:B)

EpvMat1 = melt(pvMat)
pv_stat1 = EpvMat1 %>% 
  group_by(Var2) %>%
  summarise(Median = median(value), Q900 = quantile(value, .9), Q000 = quantile(value, .0)) 

ubg6.1 = ggplot(data = pv_stat1, aes(x = Var2, y = Median) )+ 
  geom_point(col = "blue") +
  geom_errorbar(aes(ymax = Q900, ymin = Q000)) + 
  geom_point(data = R.pv1, aes(x = Var2, y = value), col = "red") +
  labs(x= 'PC', y = 'Porpotion of variance', title = title_ubg6.1) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

ubg7.1 = ggplot(data = EpvMat1, aes(x = Var2, y = value, group = Var1)) + 
  geom_line(col = "blue") +
  geom_line(data = R.pv1, aes(x = Var2, y = value), col = "red") +
  geom_point(data = R.pv1, aes(x = Var2, y = value), col = "red", size = 2) +
  labs(x= 'PC', y = 'Porpotion of variance', title = title_ubg7.1) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

print(ubg6.1)
print(ubg7.1)
# Analysis of the angle between eigenvector and batch vector

bv1 = as.numeric(unlist(train2.1 %>% select(batch)))
y1 = as.numeric(unlist(Y1))
acos(abs(cor(y1, bv1)))/pi * 180 
# [1] 78.46304
acos(abs(cor(svobj1$sv, bv1)))/pi * 180
# [1] 5.841545

R.dir1 = R.pc1$dirDf
R.cor1 = data.frame(PC = colnames(R.dir1), angle = acos(abs(cor(R.dir1, bv1)))/pi * 180 )

ubg8.1 = ggplot(data = R.cor1, aes(x = PC, y = angle))+ 
  geom_point() + 
  geom_hline(yintercep = acos(abs(cor(svobj1$sv, bv1)))/pi * 180, col = "blue", linetype = "dashed") +
  labs(x = "PC", y = "Angle", title = title_ubg8.1) + 
  scale_y_continuous(limits = c(0, 90), breaks = seq(0, 90, by = 30)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

print(ubg8.1)

# save the plots
pdf(file = 'figures/mean_heterogeneity/unbalance_case2/ubg1.1.pdf', height  = 10, width = 8)
print(ubg1.1)
dev.off()

pdf(file = 'figures/mean_heterogeneity/unbalance_case2/ubg2.1.pdf')
print(ubg2.1)
dev.off()

pdf(file = 'figures/mean_heterogeneity/unbalance_case2/ubg3.1.pdf')
print(ubg3.1)
dev.off()

pdf(file = 'figures/mean_heterogeneity/unbalance_case2/ubg4.1.pdf', height  = 10, width = 8)
print(ubg4.1)
dev.off()

pdf(file = 'figures/mean_heterogeneity/unbalance_case2/ubg5.1.pdf', height  = 10, width = 8)
print(ubg5.1)
dev.off()

pdf(file = 'figures/mean_heterogeneity/unbalance_case2/ubg6.1.pdf')
print(ubg6.1)
dev.off()

pdf(file = 'figures/mean_heterogeneity/unbalance_case2/ubg7.1.pdf')
print(ubg7.1)
dev.off()

pdf(file = 'figures/mean_heterogeneity/unbalance_case2/ubg8.1.pdf')
print(ubg8.1)
dev.off()
#############################################################################
# analysis second data set -------------------------------------------------

M2 <- train2.2 %>% select(-y, -batch, -pi_1, -pi_2)
pi_1 <- unique(train2.2$pi_1); pi_2 <- unique(train2.2$pi_2)
p_1 <- (pi_1 + pi_2)/2; p_2 <- 1 - p_1
# Expression data
Edata2 = t(M2)
colnames(Edata2) = paste0("sample", c(1: dim(Edata2)[2]))
Mdata2 = melt(Edata2)

title_ubg1.2 <- paste0("Simulation data \n pi_1 = ", pi_1, 
                       ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)
title_ubg2.2 <- paste0("Estimation of primary effects before adjustment \n pi_1 = ", pi_1, 
                       ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)
title_ubg3.2 <- paste0("Estimation of primary effects after adjustment \n pi_1 = ", pi_1, 
                       ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)
title_ubg4.2 <- paste0("Simulation data after ajustment \n pi_1 = ", pi_1, 
                       ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)
title_ubg5.2 <- paste0("Residual Matrix After Removing Primary Effect \n pi_1 = ", pi_1, 
                       ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)
title_ubg6.2 <- paste0("Eigenvalues Boxplot \n pi_1 = ", pi_1, 
                       ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)
title_ubg7.2 <- paste0("Screeplot \n pi_1 = ", pi_1, 
                       ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)
title_ubg8.2 <- paste0("Angles between PCs and Batch Effect \n blue dashed line: angle(SV, Batch Effect) \n pi_1 = ", pi_1, 
                       ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)

ubg1.2 = ggplot(Mdata2, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value", title = title_ubg1.2) +
  geom_tile() + 
  scale_fill_gradient2(limits=c(-10, 10)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())
print(ubg1.2)
# primary variable of interest
Y2 <- train2.2 %>% select(y)

# make model matrix -------------------------------------------------------
## full model 
mod2 = model.matrix(~ y, data = train2.2)
## null model
mod20 = model.matrix(~ 1, data = train2.2)


# estimate the number of latent factors that need to be estimated ---------
n.sv20 = num.sv(Edata2, mod2, method="leek")

n.sv2 = num.sv(Edata2, mod2, method="be", B = 1000)

print(paste0("The number of surrogate variable: ", n.sv2))

if(n.sv2 == 0){
  print("Mannually setting the number of surrogate variable as 1 to apply the sva algorithm")
  svobj2 = sva(t(M2), mod2, mod20, n.sv= 1, B = 5)
} else{
  svobj2 = sva(t(M2), mod2, mod20, n.sv= n.sv2, B = 5)
}


# Remove the batch effects
fsvobj2 = fsva(Edata2, mod2, svobj2, Edata2)

Edb2 = fsvobj2$db

# Estimate class effects on features from original data matrix

HatB2 = Edata2 %*% mod2 %*% solve(t(mod2) %*% mod2) 

Effect21 = HatB2 %*% matrix(c(1,1,0,1), nrow = 2, byrow = T)
rownames(Effect21) = rownames(Edata2)
colnames(Effect21) = c('yClass1', 'yClass2')
Effect21_df = data.frame(Effect21, feature = fe, x = factor(paste0('X',1:100), levels = paste0('X',1:100)))

MEff21_df = melt(Effect21_df)
ubg2.2 = ggplot(data = MEff21_df, aes(x = x, y = value, col = feature, shape = variable)) + 
  geom_point() + 
  scale_y_continuous(limits = c(-4.5, 4.5), breaks=-4:4) + 
  scale_shape_manual(values = c(4, 16)) + 
  labs(x = 'Gene',title = title_ubg2.2) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())
print(ubg2.2)
# Estimate class effects on features from original data matrix after removing the effects of surrogate variable 

HatB22 = Edb2 %*% mod2 %*% solve(t(mod2) %*% mod2) 

Effect22 = HatB22 %*% matrix(c(1,1,0,1), nrow = 2, byrow = T)
rownames(Effect22) = rownames(Edb2)
colnames(Effect22) = c('yClass1', 'yClass2')
Effect22_df = data.frame(Effect22, feature = fe, x = factor(paste0('X',1:100), levels = paste0('X',1:100)))

MEff22_df = melt(Effect22_df)
ubg3.2 = ggplot(data = MEff22_df, aes(x = x, y = value, col = feature, shape = variable)) + 
  geom_point() + 
  scale_y_continuous(limits = c(-4.5, 4.5), breaks=-4:4) + 
  scale_shape_manual(values = c(4, 16)) +
  labs(x = 'Gene',title = title_ubg3.2) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())
print(ubg3.2)

# visualize the result

Mdb2 = melt(Edb2)

ubg4.2 = ggplot(Mdb2, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",title = title_ubg4.2) +
  geom_tile() + 
  scale_fill_gradient2(limits=c(-10, 10)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

print(ubg4.2)
# dig into the analysis ---------------------------------------------------

R2 = Edata2 - HatB2 %*% t(mod2)
MR2 = melt(R2)
ubg5.2 =  ggplot(MR2, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",title = title_ubg5.2) +
  geom_tile() + 
  scale_fill_gradient2(limits=c(-10, 10)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())
print(ubg5.2)
R.pc2 = getPcaResult(R2, varNames = colnames(R2), scale=F, center = F)

R.pv2 =  data.frame(Var1 = "Origin", Var2 = rownames(R.pc2$varDf), value = R.pc2$varDf[,2])

# use permutation test to find the unusual large eigenvalues --------------

B = 1000
n1 = dim(R.pc2$varDf)[1]; n2 = dim(R2)[2]
pvMat = matrix(nrow = B, ncol = n1)
for(k in 1:B){
  # make permutation of each row independently
  newE = t(apply(R2, 1, sample, replace = FALSE))
  #  newE = apply(newE0, 2, sample, replace = FALSE)
  # refit the model to get the new residual matrix
  newR = newE - newE %*% mod2 %*% solve(t(mod2) %*% mod2) %*% t(mod2)
  colnames(newR) = c(1:n1)
  # do pca on newR and take out the proprotion variance vector
  newR.pc = getPcaResult(newR, varNames = colnames(newR), scale=F, center = F)
  pvMat[k, 1:n1] = newR.pc$varDf[,2]
}

colnames(pvMat) = rownames(R.pc2$varDf)
rownames(pvMat) = paste0("Run", 1:B)

EpvMat2 = melt(pvMat)
pv_stat2 = EpvMat2 %>% 
  group_by(Var2) %>%
  summarise(Median = median(value), Q900 = quantile(value, .9), Q000 = quantile(value, .0)) 

ubg6.2 = ggplot(data = pv_stat2, aes(x = Var2, y = Median) )+ 
  geom_point(col = "blue") +
  geom_errorbar(aes(ymax = Q900, ymin = Q000)) + 
  geom_point(data = R.pv2, aes(x = Var2, y = value), col = "red") +
  labs(x= 'PC', y = 'Porpotion of variance', title = title_ubg6.2) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

ubg7.2 = ggplot(data = EpvMat2, aes(x = Var2, y = value, group = Var1)) + 
  geom_line(col = "blue") +
  geom_line(data = R.pv2, aes(x = Var2, y = value), col = "red") +
  geom_point(data = R.pv2, aes(x = Var2, y = value), col = "red", size = 2) +
  labs(x= 'PC', y = 'Porpotion of variance', title = title_ubg7.2) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

print(ubg6.2)
print(ubg7.2)
# Analysis of the angle between eigenvector and batch vector

bv2 = as.numeric(unlist(train2.2 %>% select(batch)))
y2 = as.numeric(unlist(Y2))
acos(abs(cor(y2, bv2))) /pi * 180
# [1] 36.8699
acos(abs(cor(svobj2$sv, bv2))) /pi * 180
# [1] 82.6872
R.dir2 = R.pc2$dirDf
R.cor2 = data.frame(PC = colnames(R.dir2), angle = acos(abs(cor(R.dir2, bv2)))/pi * 180 )

ubg8.2 = ggplot(data = R.cor2, aes(x = PC, y = angle))+ 
  geom_point() + 
  geom_hline(yintercep = acos(abs(cor(svobj2$sv, bv2)))/pi * 180, col = "blue", linetype = "dashed") +
  labs(x = "PC", y = "Angle", title = title_ubg8.2) + 
  scale_y_continuous(limits = c(0, 90), breaks = seq(0, 90, by = 30)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

print(ubg8.2)

# save the plots
pdf(file = 'figures/mean_heterogeneity/unbalance_case2/ubg1.2.pdf', height  = 10, width = 8)
print(ubg1.2)
dev.off()

pdf(file = 'figures/mean_heterogeneity/unbalance_case2/ubg2.2.pdf')
print(ubg2.2)
dev.off()

pdf(file = 'figures/mean_heterogeneity/unbalance_case2/ubg3.2.pdf')
print(ubg3.2)
dev.off()

pdf(file = 'figures/mean_heterogeneity/unbalance_case2/ubg4.2.pdf', height  = 10, width = 8)
print(ubg4.2)
dev.off()

pdf(file = 'figures/mean_heterogeneity/unbalance_case2/ubg5.2.pdf', height  = 10, width = 8)
print(ubg5.2)
dev.off()

pdf(file = 'figures/mean_heterogeneity/unbalance_case2/ubg6.2.pdf')
print(ubg6.2)
dev.off()

pdf(file = 'figures/mean_heterogeneity/unbalance_case2/ubg7.2.pdf')
print(ubg7.2)
dev.off()

pdf(file = 'figures/mean_heterogeneity/unbalance_case2/ubg8.2.pdf')
print(ubg8.2)
dev.off()
