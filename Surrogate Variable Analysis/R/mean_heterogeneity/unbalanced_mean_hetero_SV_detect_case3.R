#################################################################################
## unbalanced_mean_hetero_SV_detect_case3.R
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
source('~/researchspace/robust-against-heterogeneity/Surrogate Variable Analysis/pcafuns/R/pcaScreePlot.R')
source('~/researchspace/robust-against-heterogeneity/Surrogate Variable Analysis/pcafuns/R/getPcaResult.R')

# analysis first data set -------------------------------------------------

M1 <- train3.1 %>% select(-y, -batch)
# Expression data
Edata1 = t(M1)
colnames(Edata1) = paste0("sample", c(1: dim(Edata1)[2]))
Mdata1 = melt(Edata1)

ggplot(Mdata1, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()

# primary variable of interest
Y1 <- train3.1 %>% select(y)

# make model matrix -------------------------------------------------------
## full model 
mod1 = model.matrix(~ y, data = train3.1)
## null model
mod10 = model.matrix(~ 1, data = train3.1)


# estimate the number of latent factors that need to be estimated ---------
n.sv1 = num.sv(Edata1, mod1, method="be", B = 1000)
n.sv1
# [1] 2
# it succesfully estimate the number of surogate variable

# estimate the surrogate variables
svobj1 = sva(t(M1), mod1, mod10, n.sv= n.sv1, B = 5)


# Remove the batch effects
fsvobj1 = fsva(Edata1, mod1, svobj1, Edata1)

Edb1 = fsvobj1$db

# visualize the result
trainSV1 <- data.frame(Y1, t(Edb1), batch = train3.1$batch)

ggplot(data = trainSV1, aes(x = as.factor(batch), y = X2)) + geom_boxplot()

ggplot(data = train3.1, aes(x = as.factor(batch), y = X2)) + geom_boxplot()

ggplot(data = trainSV1, aes(x = X2, group = batch, col = batch)) + 
  geom_density(linetype = "dashed") + 
  geom_density(data = train3.1, aes(x = X2, group = batch, col = batch))

Mdb1 = melt(Edb1)
ggplot(Mdb1, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",title = "Remove Effect of 1 Surrogate Variables") +
  geom_tile() + 
  scale_fill_gradient2()

# dig into the analysis ---------------------------------------------------

# estimate the coefficient of basis matrix

H1 <- mod1 %*% solve(t(mod1) %*% mod1) %*% t(mod1)
R1 <- Edata1 - t(H1 %*% t(Edata1))
MR1 = melt(R1)
ggplot(MR1, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",title = "Residual Matrix After Removing Primary Effect") +
  geom_tile() + 
  scale_fill_gradient2()

R.pc1 = getPcaResult(R1, varNames = colnames(R1), scale=F, center = F)

R.pv1 =  data.frame(Var1 = "Origin", Var2 = rownames(R.pc1$varDf), value = R.pc1$varDf[,2])



# use permutation test to find the unusual large eigenvalues --------------
B = 1000
n1 = dim(R.pc1$varDf)[1]; n2 = dim(R1)[2]
pvMat = matrix(nrow = B, ncol = n1)
tempR = R1
for(k in 1:B){
  # make permutation of each row independently
  newE = t(apply(tempR, 1, sample, replace = FALSE))
  # refit the model to get the new residual matrix
  newR = newE - newE %*% mod1 %*% solve(t(mod1) %*% mod1) %*% t(mod1)
  colnames(newR) = c(1:n1)
  # do pca on newR and take out the proprotion variance vector
  newR.pc = getPcaResult(newR, varNames = colnames(newR), scale=F, center = F)
  pvMat[k, 1:n1] = newR.pc$varDf[,2]
  temoR = newR
}

colnames(pvMat) = rownames(R.pc1$varDf)
rownames(pvMat) = paste0("Run", 1:B)

EpvMat1 = melt(pvMat)
pv_stat1 = EpvMat1 %>% 
  group_by(Var2) %>%
  summarise(Median = median(value), Q900 = quantile(value, .9), Q000 = quantile(value, .0)) 

ggplot(data = pv_stat1, aes(x = Var2, y = Median) )+ 
  geom_point(col = "blue") +
  geom_errorbar(aes(ymax = Q900, ymin = Q000)) + 
  geom_point(data = R.pv1, aes(x = Var2, y = value), col = "red")

ggplot(data = EpvMat1, aes(x = Var2, y = value, group = Var1)) + 
  geom_line(col = "blue") +
  geom_line(data = R.pv1, aes(x = Var2, y = value), col = "red") +
  geom_point(data = R.pv1, aes(x = Var2, y = value), col = "red", size = 2)

# Analysis of the angle between eigenvector and batch vector

bv1 = as.numeric(unlist(train3.1 %>% select(batch)))
y1 = as.numeric(unlist(Y1))
cor(y1, bv1)
# [1] 0.1005038
cor(svobj1$sv, bv1)
#              [,1]
# [1,]  0.843334076
# [2,] -0.005835153
R.dir1 = R.pc1$dirDf
R.cor1 = data.frame(PC = colnames(R.dir1), angle = acos(cor(R.dir1, bv1))/pi * 180 )

ggplot(data = R.cor1, aes(x = PC, y = angle))+ 
  geom_point() + 
  geom_hline(yintercep = acos(cor(svobj1$sv, bv1))/pi * 180, col = "blue", linetype = "dashed") +
  geom_hline(yintercep = 90, col = "red", linetype = "dashed") + 
  labs(x = "PC", y = "Angle", 
       title = "Angle between PCs and Batch Effect \n blue dashed line: angle(SV, Batch Effect)") + 
  scale_y_continuous(limits = c(0, 180), breaks = seq(0, 180, by = 30))


# analysis second data set -------------------------------------------------

M2 <- train3.2 %>% select(-y, -batch)
# Expression data
Edata2 = t(M2)
colnames(Edata2) = paste0("sample", c(1: dim(Edata2)[2]))
Mdata2 = melt(Edata2)

ggplot(Mdata2, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()

# primary variable of interest
Y2 <- train3.2 %>% select(y)

# make model matrix -------------------------------------------------------
## full model 
mod2 = model.matrix(~ y, data = train3.2)
## null model
mod20 = model.matrix(~ 1, data = train3.2)


# estimate the number of latent factors that need to be estimated ---------
n.sv2 = num.sv(Edata2, mod2, method="be", B = 1000)

# it fails to estimate the number of surogate variable

# estimate the surrogate variables
svobj2 = sva(t(M2), mod2, mod20, n.sv= 1, B = 5)




# Remove the batch effects
fsvobj2 = fsva(Edata2, mod2, svobj2, Edata2)

Edb2 = fsvobj2$db

# visualize the result
trainSV2 <- data.frame(Y2, t(Edb2), batch = train3.2$batch)

ggplot(data = trainSV2, aes(x = as.factor(batch), y = X2)) + geom_boxplot()

ggplot(data = train3.2, aes(x = as.factor(batch), y = X2)) + geom_boxplot()

ggplot(data = trainSV2, aes(x = X2, group = batch, col = batch)) + 
  geom_density(linetype = "dashed") + 
  geom_density(data = train3.2, aes(x = X2, group = batch, col = batch))

Mdb2 = melt(Edb2)
ggplot(Mdb2, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",title = "1 Surrogate Variables") +
  geom_tile() + 
  scale_fill_gradient2()

# dig into the analysis ---------------------------------------------------

# estimate the coefficient of basis matrix
HatB2 = Edata2 %*% mod2 %*% solve(t(mod2) %*% mod2) 
R2 = Edata2 - HatB2 %*% t(mod2)

MR2 = melt(R2)
ggplot(MR2, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",title = "Residual Matrix After Removing Primary Effect") +
  geom_tile() + 
  scale_fill_gradient2()

R.pc2 = getPcaResult(R2, varNames = colnames(R2), center = F, scale = F)

R.pv2 =  data.frame(Var1 = "Origin", Var2 = rownames(R.pc2$varDf), value = R.pc2$varDf[,2])

pcaScreePlot(pcaResult = R.pc2, title="Scree Plot")

# use permutation test to find the unusual large eigenvalues --------------
B = 1000
n1 = dim(R.pc2$varDf)[1]; n2 = dim(R2)[2]
pvMat = matrix(nrow = B, ncol = n1)
tempR = R2
for(k in 1:B){
  # make permutation of each row independently
  newE <- t(apply(tempR, 1, sample, replace = FALSE))
  # refit the model to get the new residual matrix
  newR = newE - newE %*% mod2 %*% solve(t(mod2) %*% mod2) %*% t(mod2)
  colnames(newR) = c(1:n1)  
  # do pca on newR and take out the proprotion variance vector
  newR.pc = getPcaResult(newR, varNames = colnames(newR), center = F, scale = F)
  pvMat[k, 1:n1] = newR.pc$varDf[,2]
  temoR = newR
}

colnames(pvMat) = rownames(R.pc2$varDf)
rownames(pvMat) = paste0("Run", 1:B)

EpvMat2 = melt(pvMat)
pv_stat2 = EpvMat2 %>% 
  group_by(Var2) %>%
  summarise(Median = median(value), Q900 = quantile(value, .9), Q000 = quantile(value, .0)) 

ggplot(data = pv_stat2, aes(x = Var2, y = Median) )+ 
  geom_point(col = "blue") +
  geom_errorbar(aes(ymax = Q900, ymin = Q000)) + 
  geom_point(data = R.pv2, aes(x = Var2, y = value), col = "red")

ggplot(data = EpvMat2, aes(x = Var2, y = value, group = Var1)) + 
  geom_line(col = "blue") +
  geom_line(data = R.pv2, aes(x = Var2, y = value), col = "red") +
  geom_point(data = R.pv2, aes(x = Var2, y = value), col = "red", size = 2)

# Analysis of the angle between eigenvector and batch vector

bv2 = as.numeric(unlist(train3.2 %>% select(batch)))
y2 = as.numeric(unlist(Y2))
cor(y2, bv2)
# [1] 0.3464102
cor(svobj2$sv, bv2)
# [1] 0.2584565
R.dir2 = R.pc2$dirDf
R.cor2 = data.frame(PC = colnames(R.dir2), angle = acos(cor(R.dir2, bv2))/pi * 180 )

ggplot(data = R.cor2, aes(x = PC, y = angle))+ 
  geom_point() + 
  geom_hline(yintercep = acos(cor(svobj2$sv, bv2))/pi * 180, col = "blue", linetype = "dashed") +
  geom_hline(yintercep = 90, col = "red", linetype = "dashed") + 
  labs(x = "PC", y = "Angle", 
       title = "Angle between PCs and Batch Effect \n blue dashed line: angle(SV, Batch Effect)") + 
  scale_y_continuous(limits = c(0, 180), breaks = seq(0, 180, by = 30))

