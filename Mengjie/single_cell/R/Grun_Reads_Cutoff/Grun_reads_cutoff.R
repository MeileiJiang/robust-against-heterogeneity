#######################################################################################################
## Grun_reads_cutoff.R
## Author: Meilei
#######################################################################################################
library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)
library(scales)
library(RUnit)
library(moments)
library(gridExtra)
library(scales)

load(file = "rdata/Grun_reads_cutoff.RData")
load(file = "rdata/label.RData")

source("~/workspace/bdtbio_pcafuns/R/pcaSourceList.R")
source("R/utility/mytheme.R")



# data processing ---------------------------------------------------------

# Expression Matrix

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


# Spikein Matrix

dim(Spikein)
# [1] 251  17

Spikein1 = data.frame(Spikein) 
checkEquals(rownames(Spikein1), rownames(Expression1)) 

Spikein2 = t(Spikein1)
mSpikein = melt(Spikein2)
colnames(mSpikein) = c("spikein","cell","expression")
mSpikein2 <- mSpikein %>% left_join(label.df)
mSpikein2$cell = factor(mSpikein$cell, levels = rownames(Spikein1))


# visualize data set ------------------------------------------------------

# Expression Matrix

ge.1 = ggplot(data = mExpression2, aes(x = cell, y = gene, fill = expression)) +
  labs(x = "Cell", y = "Gene", fill = "Value", title = "Heatmap of Expression value for Grun reads cutoff") +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

ge.1

ge.2 = ggplot(data = mExpression2, aes(x = cell, y = gene, fill = expression)) +
  labs(x = "Cell", y = "Gene", fill = "Value", title = "Heatmap of Expression value for Grun reads cutoff") +
  facet_wrap(~label) +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

ge.2

# Spikenin Matrix

gs.1 = ggplot(data = mSpikein2, aes(x = cell, y = spikein, fill = expression)) +
  labs(x = "Cell", y = "Spikein", fill = "Value", title = "Heatmap of Spikein for Grun reads cutoff") +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

gs.1

gs.2 = ggplot(data = mSpikein2, aes(x = cell, y = spikein, fill = expression)) +
  labs(x = "Cell", y = "Spikein", fill = "Value", title = "Heatmap of Spikein for Grun reads cutoff") +
  facet_wrap(~label) +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

gs.2


# sort the column by cell expression mean ---------------------------------

cellMean = rowMeans(Expression1)
cell.df = data.frame(cell = names(cellMean), cell_mean = cellMean)
cellRank = names(sort(cellMean)) 


# sort Expression Matrix

msortExpression = mExpression2
msortExpression$cell = factor(mExpression2$cell, levels = cellRank)



gse.1 = ggplot(data = msortExpression, aes(x = cell, y = gene, fill = expression)) +
  labs(x = "Cell", y = "Gene", fill = "Value", title = "Heatmap of Expression value for Grun reads cutoff \n 
       sorted by column gene expression mean") +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

gse.1


gse.2 = ggplot(data = msortExpression, aes(x = cell, y = gene, fill = expression)) +
  labs(x = "Cell", y = "Gene", fill = "Value", title = "Heatmap of Expression value for Grun reads cutoff \n 
       sorted by column gene expression mean") +
  facet_wrap(~label) +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

gse.2



# sort Spikein matrix

msortSpikein = mSpikein2
msortSpikein$cell = factor(mSpikein2$cell, levels = cellRank)


gss.1 = ggplot(data = msortSpikein, aes(x = cell, y = spikein, fill = expression)) +
  labs(x = "Cell", y = "Spikein", fill = "Value", title = "Heatmap of  Spikein for Grun reads cutoff \n 
       sorted by column gene expression mean") +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

gss.1

gss.2 = ggplot(data = msortSpikein, aes(x = cell, y = spikein, fill = expression)) +
  labs(x = "Cell", y = "Spikein", fill = "Value", title = "Heatmap of  Spikein for Grun reads cutoff \n 
       sorted by column gene expression mean") +
  facet_wrap(~label) +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

gss.2




# summarize the result ----------------------------------------------------

pdf("figures/Grun_reads_cutoff2.pdf", width = 14)
grid.arrange(ge.1, gs.1, ncol =2)
print(ge.2)
print(gs.2)
dev.off()

pdf("figures/grun_sorted_heatmap.pdf", width = 14)
grid.arrange( gse.1, gss.1, ncol =2)
print(gse.2)
print(gss.2)
dev.off()

