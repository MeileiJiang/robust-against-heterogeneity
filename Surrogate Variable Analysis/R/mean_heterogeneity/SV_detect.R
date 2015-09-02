######################################################################
## SV_detect.R
## Try to understand the details of permutation of Residual matrix.
##
## Meilei
######################################################################
library(dplyr)
library(ggplot2)
library(sva)

# load the data -----------------------------------------------------------

load("rdata/toyexample.RData")

M <- train100d %>% select(-y, -batch)
# Expression data
Edata = t(M)
colnames(Edata) = paste0("sample", c(1: dim(Edata)[2]))
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()

# primary variable of interest
Y <- train100d %>% select(y)

# generate the basis matrix -----------------------------------------------

## full model 
mod = model.matrix(~ y, data = train100d)
## null model
mod0 = model.matrix(~ 1, data = train100d)

# estimate the number of latent factors that need to be estimated ---------
n.sv = num.sv(Edata,mod,method="be")


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
pvMat = matrix(nrow = B, ncol = dim(R)[1])
tempR = R
for(k in 1:B){
  # make permutation of each row independently
  newE = matrix(nrow = dim(R)[1],ncol = dim(R)[2])
  for(i in 1:dim(R)[1]){
    rank = sample(dim(R)[2])
    newE[i,1:dim(R)[2]] = tempR[i,rank]
  }
  # refit the model to get the new residual matrix
  newR = newE - newE %*% mod %*% solve(t(mod) %*% mod) %*% t(mod)
  
  # do pca on newR and take out the proprotion variance vector
  newR.pc = getPcaResult(newR, varNames = colnames(newR))
  pvMat[k, 1:dim(R)[1]] = newR.pc$varDf[,2]
  temoR = newR
}

colnames(pvMat) = rownames(R)
rownames(pvMat) = paste0("Run", 1:B)

EpvMat = melt(pvMat)

ggplot(data = EpvMat, aes(x = Var2, y = value)) + 
  geom_boxplot(col = "blue") +
  geom_point(data = R.pv, aes(x = Var2, y = value), col = "red")

ggplot(data = EpvMat, aes(x = Var2, y = value, group = Var1)) + 
  geom_line(col = "blue") +
  geom_line(data = R.pv, aes(x = Var2, y = value), col = "red") +
  geom_point(data = R.pv, aes(x = Var2, y = value), col = "red", size = 2)

quant_pvMat = rep(0, dim(R)[1])
for(i in 1:dim(R)[1]){
  quant_pvMat[i] = quantile(pvMat[,i], 0.95)
}

which(R.pv$value > quant_pvMat)
# [1] 1