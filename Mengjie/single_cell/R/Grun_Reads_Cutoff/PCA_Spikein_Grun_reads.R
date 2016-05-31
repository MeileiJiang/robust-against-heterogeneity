#######################################################################################
## PCA_Spikein_Grun_reads.R
## Author: Meilei
## Principle Components Analysis on Spiken matrix Grun_reads.
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

dim(Spikein)
# [1]  251 17
checkEquals(rownames(Spikein), rownames(Expression)) 

Spikein1 = data.frame(Spikein) 
dim(Spikein1)
# [1]  251 17
# Rows are cells and columns are spikeins.


mSpikein =  melt(t(Spikein1)) # long format of Expression matrix
colnames(mSpikein) = c("spikein","cell","expression")

gs = ggplot(data = mSpikein, aes(x = cell, y = spikein, fill = expression)) +
  labs(x = "Cell", y = "Spikein", fill = "Value", title = "Heatmap of Spikein for Grun reads cutoff") +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

gs


# set up the colors and then make variables for PCA plots -----------------


# number of colors
nc = length(unique(label.df$label))
color = c(hue_pal(h=c(0,360)+15,c=100,l=65,h.start=0,direction=1)(nc)) 
colorsList = list(color)

pnames="label"
names(colorsList) = pnames
colorsList
# $label
# [1] "#F8766D" "#00BA38" "#619CFF"
colVars = pnames


## point size 3 works better on a monitor, point size 1 works better in a pdf

psize= 1           
lsize = 6
tsize = 3


# Principle Component Analysis --------------------------------------------

## Not scaling the columns.

spi.pca <- getPcaResult(pcaDf = Spikein1, varNames =  names(Spikein1), scale = F, center = F)

# scree plot

gs.scr = pcaScreePlot(spi.pca, size=1, 
                      title="Scree plot of Spikein Matrix for Grun Reads Cutoff")
print(gs.scr)


# pca scatter plot

checkTrue(all(as.character(spi.pca$caseNames)== label.df[,"cell"]))

spi.pca.plots = pcaPlots(spi.pca, label.df, colorsList, colVars,
                         pc.choices=list(1:6),sizes=psize,legendPoss=matrix(c(0.5,3.5), nrow = 2, ncol = 1),
                         legendSizes=lsize,
                         titles="Spikein Matrix",
                         titlePoss=matrix(c(0, -2), nrow = 2, ncol = 1),titleSizes=tsize,
                         showPlots=FALSE)


spi.pca.plots2 = pcaPlots(spi.pca, label.df, colorsList, colVars,
                          pc.choices=list(7:12),sizes=psize,legendPoss=matrix(c(0.5,3.5), nrow = 2, ncol = 1),
                          legendSizes=lsize,
                          titles="Spikein Matrix",
                          titlePoss=matrix(c(0, -2), nrow = 2, ncol = 1),titleSizes=tsize,
                          showPlots=FALSE)


pdf("figures/grun_spikein_pca_scatter.pdf", width = 12, height = 12)
pcaPlotter(spi.pca.plots)
pcaPlotter(spi.pca.plots2)
dev.off()

## Scaling the columns

spi.pca2 <- getPcaResult(pcaDf = Spikein1, varNames =  names(Spikein1), scale = T, center = F)

# scree plot

gs.scr2 = pcaScreePlot(spi.pca2, size=1, 
                       title="Scree plot of Scaled Spikein Matrix for Grun Reads Cutoff")
print(gs.scr2)


# pca scatter plot

checkTrue(all(as.character(spi.pca2$caseNames)== label.df[,"cell"]))

spi.pca.plots3 = pcaPlots(spi.pca2, label.df, colorsList, colVars,
                          pc.choices=list(1:6),sizes=psize,legendPoss=matrix(c(0.5,3.5), nrow = 2, ncol = 1),
                          legendSizes=lsize,
                          titles="Scaled Spikein Matrix",
                          titlePoss=matrix(c(0, -2), nrow = 2, ncol = 1),titleSizes=tsize,
                          showPlots=FALSE)

spi.pca.plots4 = pcaPlots(spi.pca, label.df, colorsList, colVars,
                          pc.choices=list(7:12),sizes=psize,legendPoss=matrix(c(0.5,3.5), nrow = 2, ncol = 1),
                          legendSizes=lsize,
                          titles="Scaled Spikein Matrix",
                          titlePoss=matrix(c(0, -2), nrow = 2, ncol = 1),titleSizes=tsize,
                          showPlots=FALSE)

pdf("figures/grun_scaled_spikein_pca_scatter.pdf", width = 12, height = 12)
pcaPlotter(spi.pca.plots3)
pcaPlotter(spi.pca.plots4)
dev.off()

pdf("figures/grun_spikein_scree_plot.pdf")
print(gs.scr)
print(gs.scr2)
dev.off()
