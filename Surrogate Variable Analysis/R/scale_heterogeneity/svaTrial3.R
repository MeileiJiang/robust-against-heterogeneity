##################################################################################
## svaTrial3.R
## Test sva on the training data which are gaussian mixture. 
## Author: Meilei Jiang
##################################################################################
library(dplyr)
library(ggplot2)
library(sva)
library(reshape2)

# data processing ---------------------------------------------------------

load("rdata/scaleToyexample.RData")

M <- train %>% select(-y, -batch)
# Expression data
Edata = t(M)
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


# adjusting for surrogate variables
pValues = f.pvalue(Edata,mod,mod0)
qValues = p.adjust(pValues,method="BH")
which(qValues < 0.05)
# 1


modSv = cbind(mod,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)
pValuesSv = f.pvalue(Edata,modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")
which(qValuesSv < 0.05)
# [1] 1 


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
