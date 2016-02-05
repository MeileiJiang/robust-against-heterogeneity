######################################################################
## SimulationIllustration.R
## Illustrate the heatmap of each type of genes.
## Author: Meilei
######################################################################
library(ggplot2)
library(dplyr)
library(reshape2)
library(sva)

source('R/mean_heterogeneity/fixed_effects_mean_hetero_sva/Generate_mean_hetero_data.R')


n = 80; # Total number of samples
c1 = 40; # Total number of samples in the Class 1
c2 = 40; # Total number of samples in the Class 2


# Type A Genes ------------------------------------------------------------

t1 = c(40, 0, 0, 0) # Settings of features
Data1 = Generate_mean_hetero_data(pi_1 = 0.5, pi_2 = 0.5, t = t1)

M <- Data1 %>% select(-y, -batch, -pi_1, -pi_2, -c1, -c2)

# Expression data
Edata <- t(M)
colnames(Edata) <- paste0("sample", c(1: dim(Edata)[2]))

# visualize the dataset
Mdata = melt(Edata)

gg1 = ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value", title = 'Heatmap of Type A Genes') +
  geom_tile() + 
  scale_fill_gradient2(limits = c(-8, 8)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

pdf(file = 'figures/Modified_SVA/Type_A_Gene.pdf', height = 5)
print(gg1)
dev.off()

# Type B Genes ------------------------------------------------------------

t2 = c(0, 0, 40, 0) # Settings of features
Data2 = Generate_mean_hetero_data(pi_1 = 0.5, pi_2 = 0.5, t = t2)

M <- Data2 %>% select(-y, -batch, -pi_1, -pi_2, -c1, -c2)

# Expression data
Edata <- t(M)
colnames(Edata) <- paste0("sample", c(1: dim(Edata)[2]))

# visualize the dataset
Mdata = melt(Edata)

gg2 = ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value", title = 'Heatmap of Type B Genes') +
  geom_tile() + 
  scale_fill_gradient2(limits = c(-8, 8)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

pdf(file = 'figures/Modified_SVA/Type_B_Gene.pdf', height = 5)
print(gg2)
dev.off()

# Type C Genes ------------------------------------------------------------

t3 = c(0, 40, 0, 0) # Settings of features
Data3 = Generate_mean_hetero_data(pi_1 = 0.5, pi_2 = 0.5, t = t3)

M <- Data3 %>% select(-y, -batch, -pi_1, -pi_2, -c1, -c2)

# Expression data
Edata <- t(M)
colnames(Edata) <- paste0("sample", c(1: dim(Edata)[2]))

# visualize the dataset
Mdata = melt(Edata)

gg3 = ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value", title = 'Heatmap of Type C Genes') +
  geom_tile() + 
  scale_fill_gradient2(limits = c(-8, 8)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

pdf(file = 'figures/Modified_SVA/Type_C_Gene.pdf', height = 5)
print(gg3)
dev.off()

# Type D Genes ------------------------------------------------------------

t4 = c(0, 0, 0, 40) # Settings of features
Data4 = Generate_mean_hetero_data(pi_1 = 0.5, pi_2 = 0.5, t = t4)

M <- Data4 %>% select(-y, -batch, -pi_1, -pi_2, -c1, -c2)

# Expression data
Edata <- t(M)
colnames(Edata) <- paste0("sample", c(1: dim(Edata)[2]))

# visualize the dataset
Mdata = melt(Edata)

gg4 = ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value", title = 'Heatmap of Type D Genes') +
  geom_tile() + 
  scale_fill_gradient2(limits = c(-8, 8)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

pdf(file = 'figures/Modified_SVA/Type_D_Gene.pdf', height = 5)
print(gg4)
dev.off()

