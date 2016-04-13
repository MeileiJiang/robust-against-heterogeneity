##############################################################################################################
## Marg_Dist_Plot_Grun.R
## Make the marginal distribution for Grun reads cutoff dataset. Rank the columns (cell) by mean.
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

load(file = "rdata/Grun_reads_cutoff.RData")
source("~/workspace/bdtbio_pcafuns/R/pcaSourceList.R")
source("R/utility/mytheme.R")


# visualize the dataset ---------------------------------------------------

# Expression

dim(Expression)
# [1]  251 2795

Expression1 = data.frame(Expression)
Expression2 = t(Expression1)
dim(Expression2)

mExpression =  melt(Expression2)
colnames(mExpression) = c("gene","cell","expression")

ge = ggplot(data = mExpression, aes(x = cell, y = gene, fill = expression)) +
  labs(x = "Cell", y = "Gene", fill = "Value", title = "Heatmap of Expression value for Grun reads cutoff") +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

ge

cellMean = colMeans(Expression2)
cell.df = data.frame(cell = names(cellMean), cell_mean = cellMean)




cellRank = names(sort(cellMean)) 



sortExpression2 = Expression2[, cellRank]
msortExpression2 = melt(sortExpression2)
colnames(msortExpression2) = c("gene","cell","expression")

gse = ggplot(data = msortExpression2, aes(x = cell, y = gene, fill = expression)) +
  labs(x = "Cell", y = "Gene", fill = "Value", title = "Heatmap of Expression value for Grun reads cutoff \n 
       sorted by column gene expression mean") +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

gse




# Spikein

dim(Spikein)
# [1] 251  17

Spikein1 = data.frame(Spikein) 
checkEquals(rownames(Spikein1), rownames(Expression1)) 

Spikein2 = t(Spikein1)
mSpikein = melt(Spikein2)
colnames(mSpikein) = c("spikein","cell","expression")

gs = ggplot(data = mSpikein, aes(x = cell, y = spikein, fill = expression)) +
  labs(x = "Cell", y = "Spikein", fill = "Value", title = "Heatmap of Spikein for Grun reads cutoff") +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

gs

sortSpikein2 = Spikein2[, cellRank]
msortSpikein2 = melt(sortSpikein2)
colnames(msortSpikein2) = c("gene","cell","expression")

gss = ggplot(data = msortSpikein2, aes(x = cell, y = gene, fill = expression)) +
  labs(x = "Cell", y = "Gene", fill = "Value", title = "Heatmap of  Spikein for Grun reads cutoff \n sorted by column gene expression mean") +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

gss

pdf("figures/grun_sorted_heatmap.pdf", width = 12)
grid.arrange( ge, gs, ncol =2)
grid.arrange( gse, gss, ncol =2)
dev.off()

# principle components analysis -------------------------------------------




# Expression Matrix



##=========================================================================

exp.svd <- svd(Expression2)
exp.d = data.frame(eigen = exp.svd$d, index = c(1:dim(Expression2)[2]))
exp.pc <- exp.svd$v[,1:16]
pairs(exp.pc)
exp.pc.df = data.frame(exp.pc)
colnames(exp.pc.df) = paste0("PC", 1:16)

ge.scr = ggplot(data = exp.d, aes(x = index, y = eigen) ) +
  geom_point(col = "blue") +
  geom_line(col = "red") +
  labs(x = "Principal Components", y = "Eigenvalue", 
       title = "Scree plot of Expression Matrix")

# Cells seem to have three clusters on first two PCs. Identify the three PCs.

exp.2pc.df = data.frame(PC1 =exp.pc.df[,1], PC2 = 2 * exp.pc.df[,2])

exp.kmeans = kmeansruns(exp.2pc.df, krange = 3,critout=T, plot = F)

label = exp.kmeans$cluster

label.df = data.frame(label = as.factor(label), exp.2pc.df, cell.df)

ggplot(data = label.df, aes(x=PC1, y=PC2, col = label)) +
  geom_point()

label.df %>%
  filter(label == 1, PC2 < -0.14) %>%
  select(cell) 

# cell
# 1  SC_serum_2
# 2 SC_serum_44
# 3 SC_serum_58
# 4 SC_serum_72

change.cell = c("SC_serum_2","SC_serum_44","SC_serum_58","SC_serum_72")

label.df = label.df %>%
  mutate(label = replace(label, cell %in% change.cell, 3)) %>%
  as.data.frame() 

ggplot(data = label.df, aes(x=PC1, y=PC2, col = label)) +
  geom_point() +
  labs(title = "clusters under PC1 and PC2 of GRUN Expression data")

label.df = label.df %>%
  select(cell, label)

# Revisit the PC marginal distribution plot

exp.pc.df2 = exp.pc.df %>%
  mutate(cell = label.df$cell) %>%
  left_join(label.df, by = c("cell"))  
  

mexp.pc.df2 = melt(exp.pc.df2, id.vars = c("cell","label"))

mexp.pc.df2$variable = factor(mexp.pc.df2$variable, levels = paste0("PC",c(1:16)))

exp.pcscore = ggplot(data = mexp.pc.df2, aes(value)) + 
  facet_wrap(~variable, scales = "free_y") +
  geom_density(aes(col  = label, group = label)) +
  geom_jitter(aes(x = value, y = 6, col = label), height = 2) 

pdf("figures/grun_exp_pc_density.pdf")
print(exp.pcscore)
dev.off()

# Marginal Distribution Plot ----------------------------------------------


pos = floor(quantile(c(1:dim(Expression2)[1]), seq(0,1,by = 1/14)))

geneMean = rowMeans(Expression2)
geneStd = apply(Expression2, 1, sd)
geneSkew = apply(Expression2, 1, skewness)
geneKurtosis = apply(Expression2, 1, kurtosis)
gene.df = data.frame(gene = names(geneMean), gene_mean = geneMean, 
                     gene_sd = geneStd, gene_skewness = geneSkew, 
                     gene_kurtosis = geneKurtosis)
plot(sort(gene.df$gene_sd))
plot(sort(gene.df$gene_skewness))
plot(sort(gene.df$gene_kurtosis))



mExpression2 = mExpression %>% left_join(label.df, by = "cell") 
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
    geom_jitter(aes(x =expression, y = 1.8, col = label), height = 1) +
    labs(x = inputgene, title  = title) 
  return(g0)
}
# sum_expression = mExpression2 %>%
#   group_by(label, gene) %>%
#   summarise( mean_expression = mean(expression)) %>%
#   group_by()
#   
# gene_rank = sum_expression %>%
#   filter(label == 3) %>%
#   arrange(mean_expression) %>%
#   mutate(rank = row_number()) %>%
#   select(gene, rank)
# 
# sum_expression = sum_expression %>% left_join(gene_rank)
# 
# ggplot(data = sum_expression, aes(x = rank, y = mean_expression)) +
#   facet_wrap(~label) +
#   geom_line() +
#   theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

gene.df.sd <- gene.df %>% arrange(gene_sd)
gene.df.sd$index <- c(1:length(geneMean))
gsd = ggplot(gene.df.sd, aes(x = index, y = gene_sd)) + 
  geom_point() +
  geom_vline(xintercept = pos, lty = "dashed") +
  labs( x = "gene", y = "standard deviation" ) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

sd.gene = gene.df.sd[pos,]

pdf("figures/grun_sd_mdplot.pdf")
grid.arrange(gsd,                     gene_mdplot(sd.gene[1,],"gene_sd"),  gene_mdplot(sd.gene[2,],"gene_sd"),  gene_mdplot(sd.gene[3,],"gene_sd"),  
             gene_mdplot(sd.gene[4,],"gene_sd"), gene_mdplot(sd.gene[5,],"gene_sd"),  gene_mdplot(sd.gene[6,],"gene_sd"),  gene_mdplot(sd.gene[7,],"gene_sd"),  
             gene_mdplot(sd.gene[8,],"gene_sd"), gene_mdplot(sd.gene[9,],"gene_sd"),  gene_mdplot(sd.gene[10,],"gene_sd"), gene_mdplot(sd.gene[11,],"gene_sd"), 
             gene_mdplot(sd.gene[12,],"gene_sd"),gene_mdplot(sd.gene[13,],"gene_sd"), gene_mdplot(sd.gene[14,],"gene_sd"), gene_mdplot(sd.gene[15,],"gene_sd"),
             ncol = 4)
dev.off()

gene.df.skewness <- gene.df %>% arrange(gene_skewness)
gene.df.skewness$index <- c(1:length(geneSkew))
gskew = ggplot(gene.df.skewness, aes(x = index, y = gene_skewness)) + 
  geom_point() +
  geom_vline(xintercept = pos, lty = "dashed") +
  labs( x = "gene", y = "skewness" ) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

skewness.gene =gene.df.skewness[pos,]

pdf("figures/grun_skew_mdplot.pdf")
grid.arrange(gskew,                     gene_mdplot(skewness.gene[1,],"gene_skewness"),  gene_mdplot(skewness.gene[2,],"gene_skewness"),  gene_mdplot(skewness.gene[3,],"gene_skewness"),  
             gene_mdplot(skewness.gene[4,],"gene_skewness"), gene_mdplot(skewness.gene[5,],"gene_skewness"),  gene_mdplot(skewness.gene[6,],"gene_skewness"),  gene_mdplot(skewness.gene[7,],"gene_skewness"),  
             gene_mdplot(skewness.gene[8,],"gene_skewness"), gene_mdplot(skewness.gene[9,],"gene_skewness"),  gene_mdplot(skewness.gene[10,],"gene_skewness"), gene_mdplot(skewness.gene[11,],"gene_skewness"), 
             gene_mdplot(skewness.gene[12,],"gene_skewness"),gene_mdplot(skewness.gene[13,],"gene_skewness"), gene_mdplot(skewness.gene[14,],"gene_skewness"), gene_mdplot(skewness.gene[15,],"gene_skewness"),
             ncol = 4)
dev.off()


gene.df.kurtosis <- gene.df %>% arrange(gene_kurtosis)
gene.df.kurtosis$index <- c(1:length(geneKurtosis))
gkurt = ggplot(gene.df.kurtosis, aes(x = index, y = gene_kurtosis)) + 
  geom_point() +
  geom_vline(xintercept = pos, lty = "dashed") +
  labs( x = "gene", y = "kurtosis" ) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

kurtosis.gene =gene.df.kurtosis[pos,]
pdf("figures/grun_kurt_mdplot.pdf", width = 12, height = 12)
grid.arrange(gkurt,                     gene_mdplot(kurtosis.gene[1,],"gene_kurtosis"),  gene_mdplot(kurtosis.gene[2,],"gene_kurtosis"),  gene_mdplot(kurtosis.gene[3,],"gene_kurtosis"),  
             gene_mdplot(kurtosis.gene[4,],"gene_kurtosis"), gene_mdplot(kurtosis.gene[5,],"gene_kurtosis"),  gene_mdplot(kurtosis.gene[6,],"gene_kurtosis"),  gene_mdplot(kurtosis.gene[7,],"gene_kurtosis"),  
             gene_mdplot(kurtosis.gene[8,],"gene_kurtosis"), gene_mdplot(kurtosis.gene[9,],"gene_kurtosis"),  gene_mdplot(kurtosis.gene[10,],"gene_kurtosis"), gene_mdplot(kurtosis.gene[11,],"gene_kurtosis"), 
             gene_mdplot(kurtosis.gene[12,],"gene_kurtosis"),gene_mdplot(kurtosis.gene[13,],"gene_kurtosis"), gene_mdplot(kurtosis.gene[14,],"gene_kurtosis"), gene_mdplot(kurtosis.gene[15,],"gene_kurtosis"),
             ncol = 4) 
dev.off()




# PC Pairs scatter plot ---------------------------------------------------

exp.pca <- getPcaResult(t(Expression2), label.df$cell, scale = F, center = F)

## set up the colors and then make variables for PCA plots

nc = 3
color = c(hue_pal(h=c(0,360)+15,c=100,l=65,h.start=0,direction=1)(nc)) 
pnames="label"
colorsList = list(color)
names(colorsList) = pnames
colorsList
colVars = pnames

## scree plot

ge.scr = pcaScreePlot(exp.pca, size=1,title="Scree plot of Expression Matrix")
print(ge.scr)

pcaLoadingPlot(exp.pca, pc.choose=1:3)

## point size 3 works better on a monitor, point size 1 works better in a pdf

psize= 2           
lsize = 5
tsize = 4

checkTrue(all(as.character(exp.pca$caseNames)== label.df[,"cell"]))

exp.pca.plots = pcaPlots(exp.pca, label.df,colorsList,colVars,
                         pc.choices=list(1:6),sizes=psize,legendPoss=matrix(c(0.5,1),2,1),
                         legendSizes=lsize,
                         titles="PCA for Expression Matrix",
                         titlePoss=matrix(c(0,-3),2,1),titleSizes=tsize,
                         showPlots=FALSE)

pdf("figures/grun_pca_scatter.pdf", width = 12, height = 12)
pcaPlotter(exp.pca.plots)
dev.off()

