##################################################################################
## svaTrial.R
## Test sva on the training data which are gaussian mixture. 
## Author: Meilei Jiang
##################################################################################
library(dplyr)
library(ggplot2)
library(sva)


# data processing ---------------------------------------------------------

load("rdata/toyexample1.RData")

M <- train %>% select(-y, -label)
rownames(M) =  c("X1", "X2")

Y <- train %>% select(y)

glm.sol = glm(y ~ .,family = 'binomial', data = train)
summary(glm.sol)
# make model matrix -------------------------------------------------------
## full model 
mod = model.matrix(~ y, data = train)
## null model
mod0 = model.matrix(~ 1, data = train)


# estimate the number of latent factors that need to be estimated ---------
n.sv = num.sv(t(M),mod,method="be")



# estimate the surrogate variables
svobj = sva(t(M),mod,mod0,n.sv=1)
sv = svobj$sv
pprob.gam = svobj$pprob.gam


# adjusting for surrogate variables
pValues = f.pvalue(t(M),mod,mod0)
qValues = p.adjust(pValues,method="BH")

modSv = cbind(mod,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)
pValuesSv = f.pvalue(t(M),modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")

