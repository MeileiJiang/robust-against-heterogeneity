#################################################################################
## SV_detect5.R
## Try to understand the details of permutation of Residual matrix.
## Data has unbalanced scale heterogeneity 
##
## Meilei
#################################################################################
library(dplyr)
library(ggplot2)
library(sva)
library(reshape2)

# data processing ---------------------------------------------------------

load("rdata/unbalanced_scale_hetero_example.RData")
source('~/researchspace/robust-against-heterogeneity/Surrogate Variable Analysis/pcafuns/R/pcaScreePlot.R')
source('~/researchspace/robust-against-heterogeneity/Surrogate Variable Analysis/pcafuns/R/getPcaResult.R')

# analysis first data set -------------------------------------------------

M <- train %>% select(-y, -batch)
# Expression data
Edata = t(M)
colnames(Edata) = paste0("sample", c(1: dim(Edata)[2]))
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()

# primary variable of interest
Y <- train %>% select(y)

# make model matrix -------------------------------------------------------
## full model 
mod = model.matrix(~ y, data = train)
## null model
mod0 = model.matrix(~ 1, data = train)


# estimate the number of latent factors that need to be estimated ---------
n.sv = num.sv(Edata,mod,method="be")

# it fails to estimate the number of surogate variable

# estimate the surrogate variables
svobj = sva(t(M),mod,mod0,n.sv= 1, B = 5)


# Remove the batch effects
fsvobj = fsva(Edata, mod, svobj, Edata)

Edb = fsvobj$db

# visualize the result
trainSV <- data.frame(Y, t(Edb), batch = train$batch)
ggplot(data = trainSV, aes(x = X2, y = X1, group = y , color = y)) +
  geom_point(aes(shape = factor(batch))) 

ggplot(data = trainSV, aes(x = as.factor(batch), y = X2)) + geom_boxplot()

ggplot(data = train, aes(x = as.factor(batch), y = X2)) + geom_boxplot()

Mdb = melt(Edb)
ggplot(Mdb, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",title = "1 Surrogate Variables") +
  geom_tile() + 
  scale_fill_gradient2()


Mdb.0 = melt(Edb[1:10,])
ggplot(Mdb.0, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",title = "1 Surrogate Variables") +
  geom_tile() + 
  scale_fill_gradient2()

# dig into the analysis ---------------------------------------------------

# estimate the coefficient of basis matrix
HatB = Edata %*% mod %*% solve(t(mod) %*% mod) 
R = Edata - HatB %*% t(mod)

R.pc = getPcaResult(R, varNames = colnames(R))
names(R.pc)

R.pc$npc
R.pv =  data.frame(Var1 = "Origin", Var2 = rownames(R.pc$varDf), value = R.pc$varDf[,2])

pcaScreePlot(pcaResult = R.pc, title="Scree Plot")

# use permutation test to find the unusual large eigenvalues --------------
B = 1000
pvMat = matrix(nrow = B, ncol = dim(R.pc$varDf)[1])
tempR = R
for(k in 1:B){
  # make permutation of each row independently
  newE = matrix(nrow = dim(R.pc$varDf)[1],ncol = dim(R)[2])
  for(i in 1:dim(R.pc$varDf)[1]){
    rank = sample(dim(R)[2])
    newE[i,1:dim(R)[2]] = tempR[i,rank]
  }
  # refit the model to get the new residual matrix
  newR = newE - newE %*% mod %*% solve(t(mod) %*% mod) %*% t(mod)
  
  # do pca on newR and take out the proprotion variance vector
  newR.pc = getPcaResult(newR, varNames = colnames(newR))
  pvMat[k, 1:dim(R.pc$varDf)[1]] = newR.pc$varDf[,2]
  temoR = newR
}

colnames(pvMat) = rownames(R.pc$varDf)
rownames(pvMat) = paste0("Run", 1:B)

EpvMat = melt(pvMat)
pv_stat = EpvMat %>% 
  group_by(Var2) %>%
  summarise(Median = median(value), Q975 = quantile(value, .975), Q025 = quantile(value, .025)) 

ggplot(data = pv_stat, aes(x = Var2, y = Median) )+ 
  geom_point(col = "blue") +
  geom_errorbar(aes(ymax = Q975, ymin = Q025)) + 
  geom_point(data = R.pv, aes(x = Var2, y = value), col = "red")

ggplot(data = EpvMat, aes(x = Var2, y = value, group = Var1)) + 
  geom_line(col = "blue") +
  geom_line(data = R.pv, aes(x = Var2, y = value), col = "red") +
  geom_point(data = R.pv, aes(x = Var2, y = value), col = "red", size = 2)


# analysis second data set -------------------------------------------------

M <- train2 %>% select(-y, -batch)
# Expression data
Edata = t(M)
colnames(Edata) = paste0("sample", c(1: dim(Edata)[2]))
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()

# primary variable of interest
Y <- train2 %>% select(y)

# make model matrix -------------------------------------------------------
## full model 
mod = model.matrix(~ y, data = train)
## null model
mod0 = model.matrix(~ 1, data = train)


# estimate the number of latent factors that need to be estimated ---------
n.sv = num.sv(Edata,mod,method="be")

# it fails to estimate the number of surogate variable

# estimate the surrogate variables
svobj = sva(t(M),mod,mod0,n.sv= 1, B = 5)




# Remove the batch effects
fsvobj = fsva(Edata, mod, svobj, Edata)

Edb = fsvobj$db

# visualize the result
trainSV <- data.frame(Y, t(Edb), batch = train$batch)
ggplot(data = trainSV, aes(x = X2, y = X1, group = y , color = y)) +
  geom_point(aes(shape = factor(batch))) 

ggplot(data = trainSV, aes(x = as.factor(batch), y = X2)) + geom_boxplot()

ggplot(data = train, aes(x = as.factor(batch), y = X2)) + geom_boxplot()

Mdb = melt(Edb)
ggplot(Mdb, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",title = "1 Surrogate Variables") +
  geom_tile() + 
  scale_fill_gradient2()


Mdb.0 = melt(Edb[1:10,])
ggplot(Mdb.0, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",title = "1 Surrogate Variables") +
  geom_tile() + 
  scale_fill_gradient2()

# dig into the analysis ---------------------------------------------------

# estimate the coefficient of basis matrix
HatB = Edata %*% mod %*% solve(t(mod) %*% mod) 
R = Edata - HatB %*% t(mod)

R.pc = getPcaResult(R, varNames = colnames(R))
names(R.pc)

R.pc$npc
R.pv =  data.frame(Var1 = "Origin", Var2 = rownames(R.pc$varDf), value = R.pc$varDf[,2])

pcaScreePlot(pcaResult = R.pc, title="Scree Plot")

# use permutation test to find the unusual large eigenvalues --------------
B = 1000
pvMat = matrix(nrow = B, ncol = dim(R.pc$varDf)[1])
tempR = R
for(k in 1:B){
  # make permutation of each row independently
  newE = matrix(nrow = dim(R.pc$varDf)[1],ncol = dim(R)[2])
  for(i in 1:dim(R.pc$varDf)[1]){
    rank = sample(dim(R)[2])
    newE[i,1:dim(R)[2]] = tempR[i,rank]
  }
  # refit the model to get the new residual matrix
  newR = newE - newE %*% mod %*% solve(t(mod) %*% mod) %*% t(mod)
  
  # do pca on newR and take out the proprotion variance vector
  newR.pc = getPcaResult(newR, varNames = colnames(newR))
  pvMat[k, 1:dim(R.pc$varDf)[1]] = newR.pc$varDf[,2]
  temoR = newR
}

colnames(pvMat) = rownames(R.pc$varDf)
rownames(pvMat) = paste0("Run", 1:B)

EpvMat = melt(pvMat)
pv_stat = EpvMat %>% 
  group_by(Var2) %>%
  summarise(Median = median(value), Q975 = quantile(value, .975), Q025 = quantile(value, .025)) 

ggplot(data = pv_stat, aes(x = Var2, y = Median) )+ 
  geom_point(col = "blue") +
  geom_errorbar(aes(ymax = Q975, ymin = Q025)) + 
  geom_point(data = R.pv, aes(x = Var2, y = value), col = "red")

ggplot(data = EpvMat, aes(x = Var2, y = value, group = Var1)) + 
  geom_line(col = "blue") +
  geom_line(data = R.pv, aes(x = Var2, y = value), col = "red") +
  geom_point(data = R.pv, aes(x = Var2, y = value), col = "red", size = 2)

# analysis third data set -------------------------------------------------

M <- train3 %>% select(-y, -batch)
# Expression data
Edata = t(M)
colnames(Edata) = paste0("sample", c(1: dim(Edata)[2]))
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()

# primary variable of interest
Y <- train3 %>% select(y)

# make model matrix -------------------------------------------------------
## full model 
mod = model.matrix(~ y, data = train)
## null model
mod0 = model.matrix(~ 1, data = train)


# estimate the number of latent factors that need to be estimated ---------
n.sv = num.sv(Edata,mod,method="be")

# it fails to estimate the number of surogate variable

# estimate the surrogate variables
svobj = sva(t(M),mod,mod0,n.sv= 1, B = 5)


# Remove the batch effects
fsvobj = fsva(Edata, mod, svobj, Edata)

Edb = fsvobj$db

# visualize the result
trainSV <- data.frame(Y, t(Edb), batch = train$batch)
ggplot(data = trainSV, aes(x = X2, y = X1, group = y , color = y)) +
  geom_point(aes(shape = factor(batch))) 

ggplot(data = trainSV, aes(x = as.factor(batch), y = X2)) + geom_boxplot()

ggplot(data = train, aes(x = as.factor(batch), y = X2)) + geom_boxplot()

Mdb = melt(Edb)
ggplot(Mdb, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",title = "1 Surrogate Variables") +
  geom_tile() + 
  scale_fill_gradient2()


Mdb.0 = melt(Edb[1:10,])
ggplot(Mdb.0, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",title = "1 Surrogate Variables") +
  geom_tile() + 
  scale_fill_gradient2()

# dig into the analysis ---------------------------------------------------

# estimate the coefficient of basis matrix
HatB = Edata %*% mod %*% solve(t(mod) %*% mod) 
R = Edata - HatB %*% t(mod)

R.pc = getPcaResult(R, varNames = colnames(R))
names(R.pc)

R.pc$npc
R.pv =  data.frame(Var1 = "Origin", Var2 = rownames(R.pc$varDf), value = R.pc$varDf[,2])

pcaScreePlot(pcaResult = R.pc, title="Scree Plot")

# use permutation test to find the unusual large eigenvalues --------------
B = 1000
pvMat = matrix(nrow = B, ncol = dim(R.pc$varDf)[1])
tempR = R
for(k in 1:B){
  # make permutation of each row independently
  newE = matrix(nrow = dim(R.pc$varDf)[1],ncol = dim(R)[2])
  for(i in 1:dim(R.pc$varDf)[1]){
    rank = sample(dim(R)[2])
    newE[i,1:dim(R)[2]] = tempR[i,rank]
  }
  # refit the model to get the new residual matrix
  newR = newE - newE %*% mod %*% solve(t(mod) %*% mod) %*% t(mod)
  
  # do pca on newR and take out the proprotion variance vector
  newR.pc = getPcaResult(newR, varNames = colnames(newR))
  pvMat[k, 1:dim(R.pc$varDf)[1]] = newR.pc$varDf[,2]
  temoR = newR
}

colnames(pvMat) = rownames(R.pc$varDf)
rownames(pvMat) = paste0("Run", 1:B)

EpvMat = melt(pvMat)
pv_stat = EpvMat %>% 
  group_by(Var2) %>%
  summarise(Median = median(value), Q975 = quantile(value, .975), Q025 = quantile(value, .025)) 

ggplot(data = pv_stat, aes(x = Var2, y = Median) )+ 
  geom_point(col = "blue") +
  geom_errorbar(aes(ymax = Q975, ymin = Q025)) + 
  geom_point(data = R.pv, aes(x = Var2, y = value), col = "red")

ggplot(data = EpvMat, aes(x = Var2, y = value, group = Var1)) + 
  geom_line(col = "blue") +
  geom_line(data = R.pv, aes(x = Var2, y = value), col = "red") +
  geom_point(data = R.pv, aes(x = Var2, y = value), col = "red", size = 2)
