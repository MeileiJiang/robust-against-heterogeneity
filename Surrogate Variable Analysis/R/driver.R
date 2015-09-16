############################################################
## Driver.R
## Understand the basic ideas of surrogate variable analysis
## Repeat the analysis in the vignette of the sva package
## Meilei
#############################################################
library(dplyr)
library(ggplot2)
library(sva)
library(bladderbatch)
library(limma)
library(reshape2)
source('~/researchspace/robust-against-heterogeneity/Surrogate Variable Analysis/pcafuns/R/getPcaResult.R')
## data processing
data(bladderdata)

pheno = pData(bladderEset)
dim(pheno)
# [1] 57  4
edata = exprs(bladderEset)
dim(edata)
#[1] 22283    57
# This is an expression set consisting of 57 samples drawn from a study of bladder cancer. The samples 
# were collected on different dates which have been used to define 5 batches. These data are used as an 
# example in the sva package vignette.
Mdata = melt(edata)
g0 = ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value", 
       title = expression(atop("Expression from a study of bladder cancer", 
                               atop(italic("22283 Genes in rows and 57 samples in columns"), "")))) +
  geom_tile() + 
  scale_fill_gradient2(midpoint = 6) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank())

pdf(file = "figures/bladder.pdf", width = 6, height = 10)
print(g0)
dev.off()

pheno %>% group_by(batch,cancer) %>%
  summarize(n = n())

## make two model matrices
# full modelcontains both the adjustment varables and the variable of interest (cancer status) 
mod = model.matrix(~as.factor(cancer), data=pheno)
# null model contains only the adjustment variables
mod0 = model.matrix(~1,data=pheno)

# estimate the number of latent factors that need to be estimated
n.sv = num.sv(edata,mod,method="be", B = 1000)
print(paste0("Number of surrogate variable by permutation test: ", n.sv))

# Look at the details of permutation test

HatB = edata %*% mod %*% solve(t(mod) %*% mod) 
R = edata - HatB %*% t(mod)

MR = melt(R)
g1 = ggplot(MR, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",
       title = expression(atop("Expression after removing primary variable effects", 
                               atop(italic("22283 Genes in rows and 57 samples in columns"), "")))) +
  geom_tile() + 
  scale_fill_gradient2(midpoint = 0) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank())

pdf(file = "figures/bladder_remove_primary.pdf", width = 6, height = 8)
print(g1)
dev.off()


R.pc = getPcaResult(R, varNames = colnames(R), scale=F, center = F)

R.pv =  data.frame(Var1 = "Origin", Var2 = rownames(R.pc$varDf), value = R.pc$varDf[,2])

# use permutation test to find the unusual large eigenvalues --------------
# permutation by row as suggested in the literature
B = 1000
n1 = dim(R.pc$varDf)[1]; n2 = dim(R)[2]
pvMat1 = matrix(nrow = B, ncol = n1)
for(k in 1:B){
  # make permutation of each row independently
  newE = t(apply(R, 1, sample, replace = FALSE))
  # refit the model to get the new residual matrix
  newR = newE - newE %*% mod %*% solve(t(mod) %*% mod) %*% t(mod)
  colnames(newR) = c(1:n1)
  # do pca on newR and take out the proprotion variance vector
  newR.pc = getPcaResult(newR, varNames = colnames(newR), scale=F, center = F)
  pvMat1[k, 1:n1] = newR.pc$varDf[,2]
}

colnames(pvMat1) = rownames(R.pc$varDf)
rownames(pvMat1) = paste0("Run", 1:B)

EpvMat1 = melt(pvMat1)
pv_stat1 = EpvMat1 %>% 
  group_by(Var2) %>%
  summarise(Median = median(value), Q950 = quantile(value, .95), Q000 = quantile(value, .0)) 

g3 = ggplot(data = pv_stat1, aes(x = Var2, y = Median) )+ 
  geom_point(col = "blue") +
  geom_errorbar(aes(ymax = Q950, ymin = Q000)) + 
  geom_point(data = R.pv, aes(x = Var2, y = value), col = "red") +
  labs(x = "PC", y = "Proportion of Variance", title = "Eigenvalue Boxplot of Row Permuted Residual Matrix")


g4 = ggplot(data = EpvMat1, aes(x = Var2, y = value, group = Var1)) + 
  geom_line(col = "blue") +
  geom_line(data = R.pv, aes(x = Var2, y = value), col = "red") +
  geom_point(data = R.pv, aes(x = Var2, y = value), col = "red", size = 2)+
  labs(x = "PC", y = "Proportion of Variance", title = "Scree plot of Row Permuted Residual Matrix")

print(g3)
print(g4)
# permutation by column and row
B = 1000
n1 = dim(R.pc$varDf)[1]; n2 = dim(R)[2]
pvMat2 = matrix(nrow = B, ncol = n1)
for(k in 1:B){
  # make permutation of each row independently
  newE = apply(R, 2, sample, replace = FALSE)  
#  newE = apply(tempR, 2, sample, replace = FALSE)
  # refit the model to get the new residual matrix
  newR = newE - newE %*% mod %*% solve(t(mod) %*% mod) %*% t(mod)
  colnames(newR) = c(1:n1)
  # do pca on newR and take out the proprotion variance vector
  newR.pc = getPcaResult(newR, varNames = colnames(newR), scale=F, center = F)
  pvMat2[k, 1:n1] = newR.pc$varDf[,2]
}

colnames(pvMat2) = rownames(R.pc$varDf)
rownames(pvMat2) = paste0("Run", 1:B)

EpvMat2 = melt(pvMat2)
pv_stat2 = EpvMat2 %>% 
  group_by(Var2) %>%
  summarise(Median = median(value), Q950 = quantile(value, .95), Q000 = quantile(value, .0)) 

g5 = ggplot(data = pv_stat2, aes(x = Var2, y = Median) )+ 
  geom_point(col = "blue") +
  geom_errorbar(aes(ymax = Q950, ymin = Q000)) + 
  geom_point(data = R.pv, aes(x = Var2, y = value), col = "red") +
  labs(x = "PC", y = "Proportion of Variance", title = "Eigenvalue Boxplot of Fully Permuted Residual Matrix")

g6 = ggplot(data = EpvMat2, aes(x = Var2, y = value, group = Var1)) + 
  geom_line(col = "blue") +
  geom_line(data = R.pv, aes(x = Var2, y = value), col = "red") +
  geom_point(data = R.pv, aes(x = Var2, y = value), col = "red", size = 2) +
  labs(x = "PC", y = "Proportion of Variance", title = "Scree plot of Fully Permuted Residual Matrix")

print(g5)
print(g6)
## Leek's method to detect the surrogate variable
n.sv = num.sv(edata,mod,method="leek")
print(paste0("Number of surrogate variable by asymptotic approach: ", n.sv))



# estimate the surrogate variables
svobj = sva(edata,mod,mod0,n.sv=n.sv)
# adjusting for surrogate variables
pValues = f.pvalue(edata,mod,mod0)
qValues = p.adjust(pValues,method="BH")

modSv = cbind(mod,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)
pValuesSv = f.pvalue(edata,modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")

# applying the Combat function to adjust for known batches
batch = pheno$batch
modcombat = model.matrix(~1, data=pheno)

save(g0, g1, g3, g4, g5, g6, file = "rdata/bladder.RData")
