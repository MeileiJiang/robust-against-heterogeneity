#################################################################################################
## angle_analysis_pc_x.R
## Try to understand the asymmetric angle between PC1 of R and x along the line pi_1 = pi_2.
##
## Author: Meilei
#################################################################################################
library(ggplot2)
library(dplyr)
library(sva)
library(reshape2)
source('R/mean_heterogeneity/fixed_effects_mean_hetero_sva/Generate_mean_hetero_data.R')


# case1: pi_1 = 0.1, pi_2 = 0.9 -------------------------------------------

Data1 = Generate_mean_hetero_data(pi_1 = 0.1, pi_2 = 0.9)
pi_1 <- unique(Data1$pi_1); pi_2 <- unique(Data1$pi_2)
p_1 <- (pi_1 + pi_2)/2; p_2 <- 1 - p_1

M <- Data1 %>% select(-y, -batch, -pi_1, -pi_2, -c1, -c2)

# Expression data
Edata <- t(M)
colnames(Edata) <- paste0("sample", c(1: dim(Edata)[2]))
Mdata <- melt(Edata)
gg1.1 <- ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value", 
       title = paste0("Simulated data (noise = 1)\n pi_1 = ", pi_1, 
                      ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)) +
  geom_tile() + 
  scale_fill_gradient2(limits=c(-10, 10)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

## full model 
mod <- model.matrix(~ y, data = Data1)
## null model
mod0 <- model.matrix(~ 1, data = Data1)

## estimate the number of surrogate variable
n.sv <- num.sv(Edata, mod, method="be", B = 100)
## estimate the surrogate varialbe
svobj <- sva(Edata, mod, mod0, n.sv= n.sv) 
## get the surrogate variable
sv <- svobj$sv

## consturct the hat matrix and get residual matrix R
HatB = Edata %*% mod %*% solve(t(mod) %*% mod) 
R = Edata - HatB %*% t(mod)

MR <- melt(R)
gg1.2 <- ggplot(MR, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",
       title = paste0("Residual Matrix After Removing Primary Effect (noise = 1)\n pi_1 = ", pi_1, 
                                ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)) +
  geom_tile() + 
  scale_fill_gradient2(limits=c(-10, 10)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

# get the PC1 from R
pc1 = svd(R)$v[,1]

## get the batch vector and class vector
bv = 2 - as.numeric(unlist(Data1 %>% select(batch)))
y = 2 - as.numeric(unlist(Data1 %>% select(y)))

cor(y, bv)
# [1] -0.8

cor(pc1, bv)
# [1] -0.5957799
acos(abs(cor(pc1, bv))) /pi * 180
# [1] 53.43175 
cor(sv,bv)
# [1] -0.2925123
acos(abs(cor(sv, bv))) /pi * 180
# [1] 72.99157
cor(pc1, sv)
# [1] -0.9311065
acos(abs(cor(pc1, sv))) /pi * 180
# [1] 21.39204

byx <- colMeans(Edata[21:40,]) # mean data row vector affected by both batch and class effects
bx <- colMeans(Edata[41:60,])  # mean data row vector affected by batch effects

rbyx <- colMeans(R[21:40,]) # mean residual data row vector afffected by both batch and class effects
rbx <- colMeans(R[41:60,]) # mean residual data row vector affected by batch effects.

mean(abs(cor(apply(R[21:40, ], 1, as.numeric), pc1)))
# [1] 0.7731743
cor(pc1, rbyx)
# [1] 0.9926291
acos(abs(cor(pc1, rbyx))) /pi * 180
# [1] 6.960916

mean(abs(cor(apply(R[41:60, ], 1, as.numeric), pc1)))
# [1] 0.7824425
cor(pc1, rbx)
# [1] 0.9926186
acos(abs(cor(pc1, rbx))) /pi * 180
# [1] 6.96585


mean(abs(cor(apply(Edata[21:40, ], 1, as.numeric), pc1)))
# [1] 0.7429105
cor(pc1, byx)
# [1] 0.9347363
acos(abs(cor(pc1, byx))) /pi * 180
# [1] 20.81441

mean(abs(cor(apply(Edata[41:60, ], 1, as.numeric), pc1)))
# [1] 0.5434639
cor(pc1, bx)
# [1] 0.6017414
acos(abs(cor(pc1, bx))) /pi * 180
# [1] 53.00528

rowvecDf <- rbind(y, bv, sqrt(40)*sv, sqrt(40)*pc1, byx, bx, rbyx, rbx)
rownames(rowvecDf) <- c('class_vector','batch_vector','surrogate_variable','pc1',
                        'gene_affected_by_batch_class','gene_affected_by_batch',
                        'residual_gene_affected_by_batch_class','residual_gene_affected_by_batch')

Mrowvec <- melt(rowvecDf)
gg1.3 <- ggplot(data = Mrowvec, aes(x = Var2, y = value)) +
  facet_wrap(~Var1, ncol = 2) + 
  geom_point() +
  scale_y_continuous(limits = c(-5, 5), breaks=c(-4,-2,-1,0,1,2,4)) + 
  labs(x = 'sample', y = 'value', 
       title = paste0("Vectors in the Row Space (noise = 1)\n pi_1 = ", pi_1, 
                      ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

pdf('figures/mean_heterogeneity/angle/angle_analysis_pc_x/gg1.1.pdf')
print(gg1.1)
dev.off()
pdf('figures/mean_heterogeneity/angle/angle_analysis_pc_x/gg1.2.pdf')
print(gg1.2)
dev.off()
pdf('figures/mean_heterogeneity/angle/angle_analysis_pc_x/gg1.3.pdf')
print(gg1.3)
dev.off()

# case2: pi_1 = 0.9, pi_2 = 0.1 -------------------------------------------

Data2 = Generate_mean_hetero_data(pi_1 = 0.9, pi_2 = 0.1)
pi_1 <- unique(Data2$pi_1); pi_2 <- unique(Data2$pi_2)
p_1 <- (pi_1 + pi_2)/2; p_2 <- 1 - p_1

M <- Data2 %>% select(-y, -batch, -pi_1, -pi_2, -c1, -c2)

# Expression data
Edata <- t(M)
colnames(Edata) <- paste0("sample", c(1: dim(Edata)[2]))
Mdata <- melt(Edata)
gg2.1 <- ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value", 
       title = paste0("Simulated data (noise = 1)\n pi_1 = ", pi_1, 
                      ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)) +
  geom_tile() + 
  scale_fill_gradient2(limits=c(-10, 10)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

## full model 
mod <- model.matrix(~ y, data = Data2)
## null model
mod0 <- model.matrix(~ 1, data = Data2)

## estimate the number of surrogate variable
n.sv <- num.sv(Edata, mod, method="be", B = 100)
## estimate the surrogate varialbe
svobj <- sva(Edata, mod, mod0, n.sv= n.sv) 
## get the surrogate variable
sv <- svobj$sv

## consturct the hat matrix and get residual matrix R
HatB <- Edata %*% mod %*% solve(t(mod) %*% mod) 
R <- Edata - HatB %*% t(mod)

MR <- melt(R)
gg2.2 <- ggplot(MR, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",
       title = paste0("Residual Matrix After Removing Primary Effect (noise = 1)\n pi_1 = ", pi_1, 
                      ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)) +
  geom_tile() + 
  scale_fill_gradient2(limits=c(-10, 10)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

# get the PC1 from R
pc1 = svd(R)$v[,1]

## get the batch vector and class vector
bv = 2 - as.numeric(unlist(Data2 %>% select(batch)))
y = 2 - as.numeric(unlist(Data2 %>% select(y)))

cor(y, bv)
# [1] 0.8

cor(pc1, bv)
# [1] 0.5944936
acos(abs(cor(pc1, bv))) / pi * 180
# [1] 53.52346

cor(sv, bv)
# [1] 0.08596965
acos(abs(cor(sv, bv))) / pi * 180
# [1] 85.06821

acos(abs(cor(sv, pc1))) / pi * 180
# [1] 72.36302

byx <- colMeans(Edata[21:40,]) # mean data row vector affected by both batch and class effects
bx <- colMeans(Edata[41:60,])  # mean data row vector affected by batch effects

rbyx <- colMeans(R[21:40,]) # mean residual data row vector afffected by both batch and class effects
rbx <- colMeans(R[41:60,]) # mean residual data row vector affected by batch effects.

mean(abs(cor(apply(R[21:40, ], 1, as.numeric), pc1)))
# [1] 0.7723611
cor(pc1, rbyx)
# [1] 0.9908265
acos(abs(cor(pc1, rbyx))) /pi * 180
# [1] 7.766705

mean(abs(cor(apply(R[41:60, ], 1, as.numeric), pc1)))
# [1] 0.7660364
cor(pc1, rbx)
# [1] 0.9905284
acos(abs(cor(pc1, rbx))) /pi * 180
# [1] 7.892098


mean(abs(cor(apply(Edata[21:40, ], 1, as.numeric), pc1)))
# [1] 0.2990941
cor(pc1, byx)
# [1] 0.3081424
acos(abs(cor(pc1, byx))) /pi * 180
# [1] 72.05268

mean(abs(cor(apply(Edata[41:60, ], 1, as.numeric), pc1)))
# [1] 0.5325497
cor(pc1, bx)
# [1] 0.5946876
acos(abs(cor(pc1, bx))) /pi * 180
# [1] 53.50964

rowvecDf <- rbind(y, bv, sqrt(40)*sv, sqrt(40)*pc1, byx, bx, rbyx, rbx)
rownames(rowvecDf) <- c('class_vector','batch_vector','surrogate_variable','pc1',
                        'gene_affected_by_batch_class','gene_affected_by_batch',
                        'residual_gene_affected_by_batch_class','residual_gene_affected_by_batch')

Mrowvec <- melt(rowvecDf)
gg2.3 <- ggplot(data = Mrowvec, aes(x = Var2, y = value)) +
  facet_wrap(~Var1, ncol = 2) + 
  geom_point() +
  labs(x = "Sample", y = "Value",
       title = paste0("Vectors in the Row Space (noise = 1)\n pi_1 = ", pi_1, 
                      ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)) +
  scale_y_continuous(limits = c(-5, 5), breaks=c(-4,-2,-1,0,1,2,4)) + 
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

pdf('figures/mean_heterogeneity/angle/angle_analysis_pc_x/gg2.1.pdf')
print(gg2.1)
dev.off()
pdf('figures/mean_heterogeneity/angle/angle_analysis_pc_x/gg2.2.pdf')
print(gg2.2)
dev.off()
pdf('figures/mean_heterogeneity/angle/angle_analysis_pc_x/gg2.3.pdf')
print(gg2.3)
dev.off()

# Case3: pi_1 = 0.5, pi_2 = 0.5 -------------------------------------------

Data3 = Generate_mean_hetero_data(pi_1 = 0.5, pi_2 = 0.5)
pi_1 <- unique(Data3$pi_1); pi_2 <- unique(Data3$pi_2)
p_1 <- (pi_1 + pi_2)/2; p_2 <- 1 - p_1

M <- Data3 %>% select(-y, -batch, -pi_1, -pi_2, -c1, -c2)

# Expression data
Edata <- t(M)
colnames(Edata) <- paste0("sample", c(1: dim(Edata)[2]))
Mdata <- melt(Edata)
gg3.1 <- ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value", 
       title = paste0("Simulated data (noise = 1)\n pi_1 = ", pi_1, 
                      ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)) +
  geom_tile() + 
  scale_fill_gradient2(limits=c(-10, 10)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

## full model 
mod <- model.matrix(~ y, data = Data3)
## null model
mod0 <- model.matrix(~ 1, data = Data3)

## estimate the number of surrogate variable
n.sv <- num.sv(Edata, mod, method="be", B = 100)
## estimate the surrogate varialbe
svobj <- sva(Edata, mod, mod0, n.sv= n.sv) 
## get the surrogate variable
sv <- svobj$sv

## consturct the hat matrix and get residual matrix R
HatB = Edata %*% mod %*% solve(t(mod) %*% mod) 
R = Edata - HatB %*% t(mod)

MR = melt(R)
gg3.2 <- ggplot(MR, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",
       title = paste0("Residual Matrix After Removing Primary Effect (noise = 1)\n pi_1 = ", pi_1, 
                      ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)) +
  geom_tile() + 
  scale_fill_gradient2(limits=c(-10, 10)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

# get the PC1 from R
pc1 = svd(R)$v[,1]

## get the batch vector and class vector
bv = 2 - as.numeric(unlist(Data3 %>% select(batch)))
y = 2 - as.numeric(unlist(Data3 %>% select(y)))

cor(bv, y)
# [1] 0

cor(pc1, bv)
# [1] 0.9961246
acos(abs(cor(pc1, bv))) / pi * 180
# [1] 5.045865
cor(sv, bv)
# [1] -0.991621
acos(abs(cor(sv, bv))) / pi * 180
# [1] 7.422269
acos(abs(cor(sv, pc1))) / pi * 180
# [1] 5.593343

byx <- colMeans(Edata[21:40,]) # mean data row vector affected by both batch and class effects
bx <- colMeans(Edata[41:60,])  # mean data row vector affected by batch effects

rbyx <- colMeans(R[21:40,]) # mean residual data row vector afffected by both batch and class effects
rbx <- colMeans(R[41:60,]) # mean residual data row vector affected by batch effects.

mean(abs(cor(apply(R[21:40, ], 1, as.numeric), pc1)))
# [1] 0.8939764
cor(pc1, rbyx)
# [1] 0.99607
acos(abs(cor(pc1, rbyx))) /pi * 180
# [1] 5.081339

mean(abs(cor(apply(R[41:60, ], 1, as.numeric), pc1)))
# [1] 0.8964634
cor(pc1, rbx)
# [1] 0.9959486
acos(abs(cor(pc1, rbx))) /pi * 180
# [1] 5.159263


mean(abs(cor(apply(Edata[21:40, ], 1, as.numeric), pc1)))
# [1] 0.6589212
cor(pc1, byx)
# [1] 0.6974893
acos(abs(cor(pc1, byx))) /pi * 180
# [1] 45.77408

mean(abs(cor(apply(Edata[41:60, ], 1, as.numeric), pc1)))
# [1] 0.8952507
cor(pc1, bx)
# [1] 0.99584
acos(abs(cor(pc1, bx))) /pi * 180
# [1] 5.227995

rowvecDf <- rbind(y, bv, sqrt(40)*sv, sqrt(40)*pc1, byx, bx, rbyx, rbx)
rownames(rowvecDf) <- c('class_vector','batch_vector','surrogate_variable','pc1',
                        'gene_affected_by_batch_class','gene_affected_by_batch',
                        'residual_gene_affected_by_batch_class','residual_gene_affected_by_batch')

Mrowvec <- melt(rowvecDf)
gg3.3 <- ggplot(data = Mrowvec, aes(x = Var2, y = value)) +
  facet_wrap(~Var1, ncol = 2) + 
  geom_point() +
  labs(x = "Sample", y = "Value",
       title = paste0("Vectors in the Row Space (noise = 1)\n pi_1 = ", pi_1, 
                      ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)) +
  scale_y_continuous(limits = c(-5, 5), breaks=c(-4,-2,-1,0,1,2,4)) + 
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

pdf('figures/mean_heterogeneity/angle/angle_analysis_pc_x/gg3.1.pdf')
print(gg3.1)
dev.off()
pdf('figures/mean_heterogeneity/angle/angle_analysis_pc_x/gg3.2.pdf')
print(gg3.2)
dev.off()
pdf('figures/mean_heterogeneity/angle/angle_analysis_pc_x/gg3.3.pdf')
print(gg3.3)
dev.off()
