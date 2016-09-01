#######################################################################################
## PCA_Expression_Grun_reads.R
## Author: Meilei
## Principle Components Analysis on Expression matrix Grun_reads.
## Date: April 2016
#######################################################################################
library(dplyr)
library(ggplot2)
library(reshape2)


load(file = "rdata/Grun_reads_cutoff.RData")
load(file = "rdata/label.RData")

# load the r functions written by John
source("~/workspace/bdtbio_pcafuns/R/pcaSourceList.R")
source("R/utility/mytheme.R")



# data processing ---------------------------------------------------------

dim(Expression)
# [1]  251 2795

Expression1 = data.frame(Expression)
dim(Expression1)
# [1]  251 2795
# Rows are cells and columns are genes.



mExpression =  melt(t(Expression1)) # long format of Expression matrix
colnames(mExpression) = c("gene","cell","expression")

ge = ggplot(data = mExpression, aes(x = cell, y = gene, fill = expression)) +
  labs(x = "Cell", y = "Gene", fill = "Value", title = "Heatmap of Expression value for Grun reads cutoff") +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

cellMean = colMeans(Expression2)
cell.df = data.frame(cell = names(cellMean), cell_mean = cellMean)


# set up the colors and then make variables for PCA plots -----------------


# number of colors
nc = length(unique(label.df$samples))
#color = c(hue_pal(h=c(0,360)+15,c=100,l=65,h.start=0,direction=1)(nc)) 
colorsList = list(c("red","darkred","green","darkgreen"))

pnames="samples"
names(colorsList) = pnames
colorsList
# $samples
# [1] "#F8766D" "#7CAE00" "#00BFC4" "#C77CFF"

colVars = pnames


## point size 3 works better on a monitor, point size 1 works better in a pdf

psize=  2          
lsize = 12
tsize = 3


# Principle Component Analysis --------------------------------------------

## Not scaling the columns.

exp.pca <- getPcaResult(pcaDf = Expression1, varNames =  names(Expression1), scale = F, center = F)

# scree plot

ge.scr = pcaScreePlot(exp.pca, size=2, 
                      title="")
print(ge.scr)


# pca scatter plot

checkTrue(all(as.character(exp.pca$caseNames)== label.df[,"cell"]))

exp.pca.plots0 = pcaPlots(exp.pca, label.df[,c(1,3)], colorsList, colVars,
                         pc.choices=list(1:4),sizes=psize,legendPoss=matrix(c(0.5,2.5), nrow = 2, ncol = 1),
                         legendSizes=lsize,
                         titles="",
                         titlePoss=matrix(c(0, -2), nrow = 2, ncol = 1),titleSizes=tsize,
                         showPlots=FALSE)

exp.pca.plots = pcaPlots(exp.pca, label.df, colorsList, colVars,
                         pc.choices=list(1:6),sizes=psize,legendPoss=matrix(c(0.5,3.5), nrow = 2, ncol = 1),
                         legendSizes=lsize,
                         titles="Expression Matrix",
                         titlePoss=matrix(c(0, -2), nrow = 2, ncol = 1),titleSizes=tsize,
                         showPlots=FALSE)


exp.pca.plots2 = pcaPlots(exp.pca, label.df, colorsList, colVars,
                         pc.choices=list(7:12),sizes=psize,legendPoss=matrix(c(0.5,3.5), nrow = 2, ncol = 1),
                         legendSizes=lsize,
                         titles="Expression Matrix",
                         titlePoss=matrix(c(0, -2), nrow = 2, ncol = 1),titleSizes=tsize,
                         showPlots=FALSE)

png("figures/Grun_reads_cutoff/expression_pca_scatter.png", width = 720, height = 720)
pcaPlotter(exp.pca.plots0)
dev.off()

png("figures/Grun_reads_cutoff/expression_scree_plot.png", width = 720, height = 720)
print(ge.scr)
dev.off()

pdf("figures/grun_expression_pca_scatter.pdf", width = 12, height = 12)
pcaPlotter(exp.pca.plots)
pcaPlotter(exp.pca.plots2)
dev.off()

## Scaling the columns

exp.pca2 <- getPcaResult(pcaDf = Expression1, varNames =  names(Expression1), scale = T, center = F)

# scree plot

ge.scr2 = pcaScreePlot(exp.pca2, size=1, 
                      title="Scree plot of Scaled Expression Matrix for Grun Reads Cutoff")
print(ge.scr2)


# pca scatter plot

checkTrue(all(as.character(exp.pca2$caseNames)== label.df[,"cell"]))

exp.pca.plots3 = pcaPlots(exp.pca2, label.df, colorsList, colVars,
                         pc.choices=list(1:6),sizes=psize,legendPoss=matrix(c(0.5,3.5), nrow = 2, ncol = 1),
                         legendSizes=lsize,
                         titles="Scaled Expression Matrix",
                         titlePoss=matrix(c(0, -2), nrow = 2, ncol = 1),titleSizes=tsize,
                         showPlots=FALSE)

exp.pca.plots4 = pcaPlots(exp.pca, label.df, colorsList, colVars,
                          pc.choices=list(7:12),sizes=psize,legendPoss=matrix(c(0.5,3.5), nrow = 2, ncol = 1),
                          legendSizes=lsize,
                          titles="Scaled Expression Matrix",
                          titlePoss=matrix(c(0, -2), nrow = 2, ncol = 1),titleSizes=tsize,
                          showPlots=FALSE)

pdf("figures/grun_scaled_expression_pca_scatter.pdf", width = 12, height = 12)
pcaPlotter(exp.pca.plots3)
pcaPlotter(exp.pca.plots4)
dev.off()

pdf("figures/grun_expression_scree_plot.pdf")
print(ge.scr)
print(ge.scr2)
dev.off()
