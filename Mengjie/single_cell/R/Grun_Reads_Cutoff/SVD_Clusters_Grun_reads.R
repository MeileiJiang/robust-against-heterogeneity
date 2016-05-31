#######################################################################################
## SVD_Clusters_Grun_reads.R
## SVD has been applied on Grun reads cutoff dataset. There seems to be three clusters. 
## Author: Meilei
## Date: April 2016
#######################################################################################
library(dplyr)
library(ggplot2)
library(reshape2)


load(file = "rdata/Grun_reads_cutoff.RData")


# Data Processing ---------------------------------------------------------

dim(Expression)
# [1]  251 2795
# Rows are cells and columns are genes.

Expression1 = data.frame(Expression)
Expression2 = t(Expression1)
dim(Expression2)

# [1] 2795  251
# Rows are genes and columns are cells.

mExpression =  melt(Expression2) # long format of Expression matrix
colnames(mExpression) = c("gene","cell","expression")

ge = ggplot(data = mExpression, aes(x = cell, y = gene, fill = expression)) +
  labs(x = "Cell", y = "Gene", fill = "Value", title = "Heatmap of Expression value for Grun reads cutoff") +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())

cellMean = colMeans(Expression2)
cell.df = data.frame(cell = names(cellMean), cell_mean = cellMean)

# Singular Value Decomposition and Clustering Labels ----------------------

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

ge.scr

# Cells seem to have three clusters on first two PCs. Identify the three PCs.

exp.2pc.df = data.frame(PC1 =exp.pc.df[,1], PC2 = 2 * exp.pc.df[,2])

exp.kmeans = kmeansruns(exp.2pc.df, krange = 3,critout=T, plot = F)
# 3  clusters  483.6274

label = exp.kmeans$cluster
label

label.df = data.frame(cell.df, label = as.factor(label), exp.2pc.df)

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
  mutate(label = replace(label, cell %in% change.cell, 2)) %>%
  as.data.frame() 

ggplot(data = label.df, aes(x=PC1, y=PC2, col = label)) +
  geom_point() +
  labs(title = "clusters under PC1 and PC2 of GRUN Expression data")

label.df = label.df %>%
  select(cell, label)

table(label.df$label)
#  1   2   3 
# 57  45 149 


save(label.df, file = "rdata/label.RData")


# revisit the heatmap with labels -----------------------------------------

mExpression2 <- mExpression %>% left_join(label.df, by = "cell")

ge = ggplot(data = mExpression2, aes(x = cell, y = gene, fill = expression)) +
  facet_wrap(~label)+
  labs(x = "Cell", y = "Gene", fill = "Value", title = "Heatmap of Expression value for Grun reads cutoff") +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())
ge

