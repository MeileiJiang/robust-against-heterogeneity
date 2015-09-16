#################################################################################
## balanced_mean_hetero_SV_detect.R
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

M <- train0 %>% select(-y, -batch, -pi_1, -pi_2 )
pi_1 <- unique(train4$pi_1); pi_2 <- unique(train4$pi_2)
p_1 <- (pi_1 + pi_2)/2; p_2 <- 1 - p_1
# Expression data
Edata = t(M)
colnames(Edata) = paste0("sample", c(1: dim(Edata)[2]))
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()

# primary variable of interest
Y <- train0 %>% select(y)

# make model matrix -------------------------------------------------------
## full model 
mod = model.matrix(~ y, data = train0)
## null model
mod0 = model.matrix(~ 1, data = train0)


# estimate the number of latent factors that need to be estimated ---------
n.sv0 = num.sv(Edata,mod,method="leek")

n.sv = num.sv(Edata,mod,method="be", B = 1000)

print(paste0("The number of surrogate variable: ", n.sv))

if(n.sv == 0){
  print("Mannually setting the number of surrogate variable as 1 to apply the sva algorithm")
  svobj = sva(t(M),mod, mod0,n.sv= 1, B = 5)
} else{
  svobj = sva(t(M),mod, mod0,n.sv= n.sv, B = 5)
}



# Remove the batch effects
fsvobj = fsva(Edata, mod, svobj, Edata)

Edb = fsvobj$db

# visualize the result
trainSV <- data.frame(Y, t(Edb), batch = train0$batch)

ggplot(data = trainSV, aes(x = as.factor(batch), y = X31)) + geom_boxplot()

ggplot(data = train0, aes(x = as.factor(batch), y = X31)) + geom_boxplot()

ggplot(data = trainSV, aes(x = X31, group = batch, col = batch)) + 
  geom_density(linetype = "dashed") + 
  geom_density(data = train0, aes(x = X2, group = batch, col = batch))

Mdb = melt(Edb)
ggplot(Mdb, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",title = "1 Surrogate Variables") +
  geom_tile() + 
  scale_fill_gradient2()

# dig into the analysis ---------------------------------------------------

# estimate the coefficient of basis matrix
HatB = Edata %*% mod %*% solve(t(mod) %*% mod) 
R = Edata - HatB %*% t(mod)

MR = melt(R)
ggplot(MR, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",title = "1 Surrogate Variables") +
  geom_tile() + 
  scale_fill_gradient2()

R.pc = getPcaResult(R, varNames = colnames(R), scale=F, center = F)

R.pv =  data.frame(Var1 = "Origin", Var2 = rownames(R.pc$varDf), value = R.pc$varDf[,2])


# use permutation test to find the unusual large eigenvalues --------------
B = 1000
n1 = dim(R.pc$varDf)[1]; n2 = dim(R)[2]
pvMat = matrix(nrow = B, ncol = n1)
for(k in 1:B){
  # make permutation of each row independently
  # newE = apply(R, 2, sample, replace = FALSE)
  newE = t(apply(R, 1, sample, replace = FALSE))
  # refit the model to get the new residual matrix
  newR = newE - newE %*% mod %*% solve(t(mod) %*% mod) %*% t(mod)
  colnames(newR) = c(1:n1)
  # do pca on newR and take out the proprotion variance vector
  newR.pc = getPcaResult(newR, varNames = colnames(newR), scale=F, center = F)
  pvMat[k, 1:n1] = newR.pc$varDf[,2]
}

colnames(pvMat) = rownames(R.pc$varDf)
rownames(pvMat) = paste0("Run", 1:B)

EpvMat = melt(pvMat)
pv_stat = EpvMat %>% 
  group_by(Var2) %>%
  summarise(Median = median(value), Q950 = quantile(value, .95), Q000 = quantile(value, .0)) 

ggplot(data = pv_stat, aes(x = Var2, y = Median) )+ 
  geom_point(col = "blue") +
  geom_errorbar(aes(ymax = Q950, ymin = Q000)) + 
  geom_point(data = R.pv, aes(x = Var2, y = value), col = "red")

ggplot(data = EpvMat, aes(x = Var2, y = value, group = Var1)) + 
  geom_line(col = "blue") +
  geom_line(data = R.pv, aes(x = Var2, y = value), col = "red") +
  geom_point(data = R.pv, aes(x = Var2, y = value), col = "red", size = 2)

# Analysis of the angle between eigenvector and batch vector

bv = as.numeric(unlist(train0 %>% select(batch)))
y = as.numeric(unlist(Y))
cor(y, bv)
# [1] 0
cor(svobj$sv, bv)
# [1] -0.991649
R.dir = R.pc$dirDf
R.cor = data.frame(PC = colnames(R.dir), angle = acos(abs(cor(R.dir, bv)))/pi * 180 )


ggplot(data = R.cor, aes(x = PC, y = angle))+ 
  geom_point() + 
  geom_hline(yintercep = acos(abs(cor(svobj$sv, bv)))/pi * 180, col = "blue", linetype = "dashed") +
  labs(x = "PC", y = "Angle", 
       title = "Angle between PCs and Batch Effect \n blue dashed line: angle(SV, Batch Effect)") + 
  scale_y_continuous(limits = c(0, 90), breaks = seq(0, 90, by = 30))


