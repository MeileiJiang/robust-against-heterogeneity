#####################################################################################
## Non-iterative-JIVE_matlab.R
## Author: Meilei
## Visualized the non-iterative-JIVE from matlab.
#####################################################################################

library(R.matlab)
library(r.jive)
load(file = "rdata/label.RData")
load(file = "rdata/grun_reads_cutoff_jive.RData")
mat_result1 = readMat("rdata/grun_reads_cutoff_reasult1.mat")
# mat_result2 = readMat("rdata/grun_reads_cutoff_reasult2.mat")
# mat_result3 = readMat("rdata/grun_reads_cutoff_reasult3.mat")

jive_mat_result1 =  jive.result
jive_mat_result1$data = list(Expression = t(Expression), Spikein = t(Spikein))
jive_mat_result1$joint = list(mat_result1$joint.Expression, mat_result1$joint.Spikein)
jive_mat_result1$individual = list(mat_result1$indiv.Expression, mat_result1$indiv.Spikein)
jive_mat_result1$rankJ = mat_result1$rjoint
jive_mat_result1$rankA = c(mat_result1$rIndExpression, mat_result1$rIndSpikein)



pdf(file = "figures/Grun_reads_cutoff/jive/VarExplained1.pdf", height = 5, width = 8)
showVarExplained(jive_mat_result1)
dev.off()
png(file = "figures/Grun_reads_cutoff/jive/JIVE_Heatmap1.png",height=550,width=705)
showHeatmaps(jive_mat_result1, order_by=-1)
dev.off()


# make PCA plot -----------------------------------------------------------

colorlist = c("red", "green")
shapelist = c('O', "X")

names(colorlist) = unique(label.df$sample)
names(shapelist) = unique(label.df$medium)

Colors = colorlist[label.df$sample]
Shapes = shapelist[label.df$medium]

png("figures/Grun_reads_cutoff/jive/JointPCA_full1.png", height=1200, width=1200)
showPCA(jive_mat_result1, n_joint = 0,n_indiv =c(5,0),Colors=Colors, pch = Shapes)
dev.off()


