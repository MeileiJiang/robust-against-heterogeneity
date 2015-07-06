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
## data processing
data(bladderdata)

pheno = pData(bladderEset)
edata = exprs(bladderEset)

pheno %>% group_by(batch,cancer) %>%
  summarize(n = n())

## make two model matrices
# full modelcontains both the adjustment varables and the variable of interest (cancer status) 
mod = model.matrix(~as.factor(cancer), data=pheno)
# null model contains only the adjustment variables
mod0 = model.matrix(~1,data=pheno)

# estimate the number of latent factors that need to be estimated
n.sv = num.sv(edata,mod,method="leek")
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
