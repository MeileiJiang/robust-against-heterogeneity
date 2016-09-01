############################################################################################
## grun_reads_cutoff_jive.R
## Author: Meilei
## Apply JIVE on grun_reads_cutoff.RData
############################################################################################

library(ggplot2)
library(dplyr)
library(reshape2)
library(r.jive)

load(file = "rdata/Grun_reads_cutoff.RData")
load(file = "rdata/label.RData")

# load the r functions written by John
source("~/workspace/bdtbio_pcafuns/R/pcaSourceList.R")
source("R/utility/mytheme.R")

# data processing ---------------------------------------------------------

dim(Expression)
# [1]  251 2795

Expression1 = t(Expression)
dim(Expression1)
# [1] 2795  251
# Rows are genes and columns are cells.

Spikein1 = t(Spikein)
dim(Spikein1)
# [1]  17 251

Data = list(Expression = Expression1, Spikein = Spikein1)

jive.result = jive(Data)

pdf(file = "figures/Grun_reads_cutoff/jive/VarExplained.pdf", height = 5, width = 8)
showVarExplained(jive.result)
dev.off()
png(file = "figures/Grun_reads_cutoff/jive/JIVE_Heatmap.png",height=465,width=705)
showHeatmaps(jive.result, order_by=-1)
dev.off()


# make PCA plot -----------------------------------------------------------

colorlist = c("red","darkred","green","skyblue")
colorlist = c("red", "green")
shapelist = c("*", "+")

names(colorlist) = unique(label.df$sample)
names(shapelist) = unique(label.df$medium)

Colors = colorlist[label.df$sample]
Shapes = shapelist[label.df$medium]
png("figures/Grun_reads_cutoff/jive/JointPCA_full0.png",height=800,width=800)
showPCA(jive.result, n_joint = 1,n_indiv =c(7,1),Colors=Colors, pch = Shapes)
dev.off()

png("figures/Grun_reads_cutoff/jive/JointPCA_joint_spikein.png",height=800,width=800)
showPCA(jive.result, n_joint = 1,n_indiv =c(0,1),Colors=Colors)
dev.off()

png("figures/Grun_reads_cutoff/jive/JointPCA_indiv_expression0.png",height=800,width=800)
showPCA(jive.result, n_joint = 0,n_indiv =c(7,0),Colors=Colors)
dev.off()

png("figures/Grun_reads_cutoff/jive/JointPCA_cluster.png",height=800,width=800)
showPCA(jive.result, n_joint = 1,n_indiv =c(2,0),Colors=Colors)
dev.off()


# save result -------------------------------------------------------------

save(jive.result, file = "rdata/grun_reads_cutoff_jive.RData")
