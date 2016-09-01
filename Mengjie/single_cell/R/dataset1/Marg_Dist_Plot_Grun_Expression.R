##############################################################################################################
## Marg_Dist_Plot_Grun_Expression.R
## Make the marginal distribution for Expression Matrix of Grun reads cutoff dataset. 
## Author: Meilei Jiang
## Date: April 2016
##############################################################################################################

library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)
library(scales)
library(RUnit)
library(fpc)
library(moments)
library(gridExtra)
library(scales)

load(file = "rdata/Grun_reads_cutoff.RData")
load(file = "rdata/label.RData")

source("~/workspace/bdtbio_pcafuns/R/pcaSourceList.R")
source("R/utility/mytheme.R")


# visualize the dataset ---------------------------------------------------


dim(Expression)
# [1]  251 2795

Expression1 = data.frame(Expression)
dim(Expression1)
# [1]  251 2795
# Rows are cells and columns are genes.

mExpression =  melt(t(Expression1)) # long format of Expression matrix
colnames(mExpression) = c("gene","cell","expression")
mExpression2 <- mExpression %>% left_join(label.df)
mExpression2$cell = factor(mExpression$cell, levels = rownames(Expression1))


# Marginal Distribution Plot ----------------------------------------------

gene_mdplot = function(input, varname = NULL){
  inputgene = as.character(input$gene) 
  subExpression = mExpression2 %>% filter(gene == inputgene) 
  title = ""
  if(!is.na(varname)){
    value = round(input[,varname],5)
    title = paste(varname,"=",value)
  }
  
  
  g0 = ggplot(data = subExpression, aes(expression)) + 
    geom_density(aes(group = label, col  = label)) +
    geom_jitter(aes(x =expression, y = 1.8, col = label), height = 1, size = 1) +
    labs(x = inputgene, title  = title) 
  return(g0)
}

mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

pos = floor(quantile(c(1:dim(Expression1)[2]), seq(0,1,by = 1/14)))

geneMean = colMeans(Expression1)
geneStd = apply(Expression1, 2, sd)
geneSkew = apply(Expression1, 2, skewness)
geneKurtosis = apply(Expression1, 2, kurtosis)
geneMode = apply(Expression1, 2, mode)


gene.df = data.frame(gene = names(geneMean), gene_mean = geneMean, 
                     gene_sd = geneStd, gene_skewness = geneSkew, 
                     gene_kurtosis = geneKurtosis, gene_mode = geneMode)
plot(sort(gene.df$gene_sd))
plot(sort(gene.df$gene_skewness))
plot(sort(gene.df$gene_kurtosis))
plot(sort(gene.df$gene_mode))



# Standard deviation

gene.df.sd <- gene.df %>% arrange(gene_sd)
gene.df.sd$index <- c(1:length(geneMean))
pos.sd = c((dim(Expression1)[2]-14):dim(Expression1)[2])
sd.gene = gene.df.sd[pos.sd,]

gsd = ggplot(gene.df.sd, aes(x = index, y = gene_sd)) + 
  geom_point() +
  geom_vline(xintercept = pos.sd, lty = "dashed") +
  labs( x = "gene", y = "standard deviation" ) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())


pdf("figures/grun_sd_mdplot2.pdf", width = 12, height = 10)
grid.arrange(gsd,                     gene_mdplot(sd.gene[1,],"gene_sd"),  gene_mdplot(sd.gene[2,],"gene_sd"),  gene_mdplot(sd.gene[3,],"gene_sd"),  
             gene_mdplot(sd.gene[4,],"gene_sd"), gene_mdplot(sd.gene[5,],"gene_sd"),  gene_mdplot(sd.gene[6,],"gene_sd"),  gene_mdplot(sd.gene[7,],"gene_sd"),  
             gene_mdplot(sd.gene[8,],"gene_sd"), gene_mdplot(sd.gene[9,],"gene_sd"),  gene_mdplot(sd.gene[10,],"gene_sd"), gene_mdplot(sd.gene[11,],"gene_sd"), 
             gene_mdplot(sd.gene[12,],"gene_sd"),gene_mdplot(sd.gene[13,],"gene_sd"), gene_mdplot(sd.gene[14,],"gene_sd"), gene_mdplot(sd.gene[15,],"gene_sd"),
             ncol = 4)
dev.off()


# Skewness

gene.df.skewness <- gene.df %>% arrange(gene_skewness)
gene.df.skewness$index <- c(1:length(geneSkew))
pos.skew = c(1:15)
skewness.gene =gene.df.skewness[pos.skew,]

gskew = ggplot(gene.df.skewness, aes(x = index, y = gene_skewness)) + 
  geom_point() +
  geom_vline(xintercept = pos.skew, lty = "dashed") +
  labs( x = "gene", y = "skewness" ) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())


pdf("figures/grun_skew_mdplot2.pdf", width = 12, height = 10)
grid.arrange(gskew,                     gene_mdplot(skewness.gene[1,],"gene_skewness"),  gene_mdplot(skewness.gene[2,],"gene_skewness"),  gene_mdplot(skewness.gene[3,],"gene_skewness"),  
             gene_mdplot(skewness.gene[4,],"gene_skewness"), gene_mdplot(skewness.gene[5,],"gene_skewness"),  gene_mdplot(skewness.gene[6,],"gene_skewness"),  gene_mdplot(skewness.gene[7,],"gene_skewness"),  
             gene_mdplot(skewness.gene[8,],"gene_skewness"), gene_mdplot(skewness.gene[9,],"gene_skewness"),  gene_mdplot(skewness.gene[10,],"gene_skewness"), gene_mdplot(skewness.gene[11,],"gene_skewness"), 
             gene_mdplot(skewness.gene[12,],"gene_skewness"),gene_mdplot(skewness.gene[13,],"gene_skewness"), gene_mdplot(skewness.gene[14,],"gene_skewness"), gene_mdplot(skewness.gene[15,],"gene_skewness"),
             ncol = 4)
dev.off()


# Kurtosis

gene.df.kurtosis <- gene.df %>% arrange(gene_kurtosis)
gene.df.kurtosis$index <- c(1:length(geneKurtosis))
pos.kurt = c(1:15)
kurtosis.gene =gene.df.kurtosis[pos.kurt,]

gkurt = ggplot(gene.df.kurtosis, aes(x = index, y = gene_kurtosis)) + 
  geom_point() +
  geom_vline(xintercept = pos.kurt, lty = "dashed") +
  labs( x = "gene", y = "kurtosis" ) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

pdf("figures/grun_kurt_mdplot2.pdf", width = 12, height = 10)
grid.arrange(gkurt,                     gene_mdplot(kurtosis.gene[1,],"gene_kurtosis"),  gene_mdplot(kurtosis.gene[2,],"gene_kurtosis"),  gene_mdplot(kurtosis.gene[3,],"gene_kurtosis"),  
             gene_mdplot(kurtosis.gene[4,],"gene_kurtosis"), gene_mdplot(kurtosis.gene[5,],"gene_kurtosis"),  gene_mdplot(kurtosis.gene[6,],"gene_kurtosis"),  gene_mdplot(kurtosis.gene[7,],"gene_kurtosis"),  
             gene_mdplot(kurtosis.gene[8,],"gene_kurtosis"), gene_mdplot(kurtosis.gene[9,],"gene_kurtosis"),  gene_mdplot(kurtosis.gene[10,],"gene_kurtosis"), gene_mdplot(kurtosis.gene[11,],"gene_kurtosis"), 
             gene_mdplot(kurtosis.gene[12,],"gene_kurtosis"),gene_mdplot(kurtosis.gene[13,],"gene_kurtosis"), gene_mdplot(kurtosis.gene[14,],"gene_kurtosis"), gene_mdplot(kurtosis.gene[15,],"gene_kurtosis"),
             ncol = 4) 
dev.off()


# Mode

gene.df.mode <- gene.df %>% arrange(gene_mode)
gene.df.mode$index <- c(1:length(geneMode))
gmode = ggplot(gene.df.mode, aes(x = index, y = gene_mode)) + 
  geom_point() +
  geom_vline(xintercept = pos, lty = "dashed") +
  labs( x = "gene", y = "mode" ) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

mode.gene =gene.df.mode[pos,]
pdf("figures/grun_mode_mdplot.pdf", width = 12, height = 10)
grid.arrange(gmode,                     gene_mdplot(mode.gene[1,],"gene_mode"),  gene_mdplot(mode.gene[2,],"gene_mode"),  gene_mdplot(mode.gene[3,],"gene_mode"),  
             gene_mdplot(mode.gene[4,],"gene_mode"), gene_mdplot(mode.gene[5,],"gene_mode"),  gene_mdplot(mode.gene[6,],"gene_mode"),  gene_mdplot(mode.gene[7,],"gene_mode"),  
             gene_mdplot(mode.gene[8,],"gene_mode"), gene_mdplot(mode.gene[9,],"gene_mode"),  gene_mdplot(mode.gene[10,],"gene_mode"), gene_mdplot(mode.gene[11,],"gene_mode"), 
             gene_mdplot(mode.gene[12,],"gene_mode"),gene_mdplot(mode.gene[13,],"gene_mode"), gene_mdplot(mode.gene[14,],"gene_mode"), gene_mdplot(mode.gene[15,],"gene_mode"),
             ncol = 4) 
dev.off()


pdf("figures/grun_expression_mdplot.pdf", width = 12, height = 10)
grid.arrange(gsd,                     gene_mdplot(sd.gene[1,],"gene_sd"),  gene_mdplot(sd.gene[2,],"gene_sd"),  gene_mdplot(sd.gene[3,],"gene_sd"),  
             gene_mdplot(sd.gene[4,],"gene_sd"), gene_mdplot(sd.gene[5,],"gene_sd"),  gene_mdplot(sd.gene[6,],"gene_sd"),  gene_mdplot(sd.gene[7,],"gene_sd"),  
             gene_mdplot(sd.gene[8,],"gene_sd"), gene_mdplot(sd.gene[9,],"gene_sd"),  gene_mdplot(sd.gene[10,],"gene_sd"), gene_mdplot(sd.gene[11,],"gene_sd"), 
             gene_mdplot(sd.gene[12,],"gene_sd"),gene_mdplot(sd.gene[13,],"gene_sd"), gene_mdplot(sd.gene[14,],"gene_sd"), gene_mdplot(sd.gene[15,],"gene_sd"),
             ncol = 4)
grid.arrange(gskew,                     gene_mdplot(skewness.gene[1,],"gene_skewness"),  gene_mdplot(skewness.gene[2,],"gene_skewness"),  gene_mdplot(skewness.gene[3,],"gene_skewness"),  
             gene_mdplot(skewness.gene[4,],"gene_skewness"), gene_mdplot(skewness.gene[5,],"gene_skewness"),  gene_mdplot(skewness.gene[6,],"gene_skewness"),  gene_mdplot(skewness.gene[7,],"gene_skewness"),  
             gene_mdplot(skewness.gene[8,],"gene_skewness"), gene_mdplot(skewness.gene[9,],"gene_skewness"),  gene_mdplot(skewness.gene[10,],"gene_skewness"), gene_mdplot(skewness.gene[11,],"gene_skewness"), 
             gene_mdplot(skewness.gene[12,],"gene_skewness"),gene_mdplot(skewness.gene[13,],"gene_skewness"), gene_mdplot(skewness.gene[14,],"gene_skewness"), gene_mdplot(skewness.gene[15,],"gene_skewness"),
             ncol = 4)
grid.arrange(gkurt,                     gene_mdplot(kurtosis.gene[1,],"gene_kurtosis"),  gene_mdplot(kurtosis.gene[2,],"gene_kurtosis"),  gene_mdplot(kurtosis.gene[3,],"gene_kurtosis"),  
             gene_mdplot(kurtosis.gene[4,],"gene_kurtosis"), gene_mdplot(kurtosis.gene[5,],"gene_kurtosis"),  gene_mdplot(kurtosis.gene[6,],"gene_kurtosis"),  gene_mdplot(kurtosis.gene[7,],"gene_kurtosis"),  
             gene_mdplot(kurtosis.gene[8,],"gene_kurtosis"), gene_mdplot(kurtosis.gene[9,],"gene_kurtosis"),  gene_mdplot(kurtosis.gene[10,],"gene_kurtosis"), gene_mdplot(kurtosis.gene[11,],"gene_kurtosis"), 
             gene_mdplot(kurtosis.gene[12,],"gene_kurtosis"),gene_mdplot(kurtosis.gene[13,],"gene_kurtosis"), gene_mdplot(kurtosis.gene[14,],"gene_kurtosis"), gene_mdplot(kurtosis.gene[15,],"gene_kurtosis"),
             ncol = 4) 
grid.arrange(gmode,                     gene_mdplot(mode.gene[1,],"gene_mode"),  gene_mdplot(mode.gene[2,],"gene_mode"),  gene_mdplot(mode.gene[3,],"gene_mode"),  
             gene_mdplot(mode.gene[4,],"gene_mode"), gene_mdplot(mode.gene[5,],"gene_mode"),  gene_mdplot(mode.gene[6,],"gene_mode"),  gene_mdplot(mode.gene[7,],"gene_mode"),  
             gene_mdplot(mode.gene[8,],"gene_mode"), gene_mdplot(mode.gene[9,],"gene_mode"),  gene_mdplot(mode.gene[10,],"gene_mode"), gene_mdplot(mode.gene[11,],"gene_mode"), 
             gene_mdplot(mode.gene[12,],"gene_mode"),gene_mdplot(mode.gene[13,],"gene_mode"), gene_mdplot(mode.gene[14,],"gene_mode"), gene_mdplot(mode.gene[15,],"gene_mode"),
             ncol = 4) 
dev.off()
