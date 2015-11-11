#################################################################################################
## angle_analysis_pc_x2.R
## Try to understand the asymmetric angle between PC1 of R and x along the line pi_1 = pi_2.
## Look at the case where data has no random noise
##
## Author: Meilei
#################################################################################################
library(ggplot2)
library(dplyr)
library(sva)
library(reshape2)
source('R/mean_heterogeneity/fixed_effects_mean_hetero_sva/Generate_mean_hetero_data.R')


# case1: pi_1 = 0.1, pi_2 = 0.9 -------------------------------------------

Data4 = Generate_mean_hetero_data(pi_1 = 0.1, pi_2 = 0.9, noise = 0.01)
pi_1 <- unique(Data4$pi_1); pi_2 <- unique(Data4$pi_2)
p_1 <- (pi_1 + pi_2)/2; p_2 <- 1 - p_1

M <- Data4 %>% select(-y, -batch, -pi_1, -pi_2, -c1, -c2)

# Expression data
Edata <- t(M)
colnames(Edata) <- paste0("sample", c(1: dim(Edata)[2]))
Mdata <- melt(Edata)
gg4.1 <- ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value", 
       title = paste0("Simulated data (noise = 0.01)\n pi_1 = ", pi_1, 
                      ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)) +
  geom_tile() + 
  scale_fill_gradient2(limits=c(-10, 10)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

## full model 
mod <- model.matrix(~ y, data = Data4)
## null model
mod0 <- model.matrix(~ 1, data = Data4)

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
gg4.2 <- ggplot(MR, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",
       title = paste0("Residual Matrix After Removing Primary Effect (noise = 0.01) \n pi_1 = ", pi_1, 
                      ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)) +
  geom_tile() + 
  scale_fill_gradient2(limits=c(-10, 10)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

# get the PC1 from R
pc1 = svd(R)$v[,1]

## get the batch vector and class vector
bv = 2 - as.numeric(unlist(Data4 %>% select(batch)))
y = 2 - as.numeric(unlist(Data4 %>% select(y)))

cor(y, bv)
# [1] -0.8

cor(pc1, bv)
# [1] 0.5999995
acos(abs(cor(pc1, bv))) /pi * 180
# [1] 53.13014
cor(sv,bv)
# [1] -0.3166094
acos(abs(cor(sv, bv))) /pi * 180
# [1] 71.542
acos(abs(cor(sv, pc1))) /pi * 180
# [1] 18.41202

byx <- colMeans(Edata[21:40,]) # mean data row vector affected by both batch and class effects
bx <- colMeans(Edata[41:60,])  # mean data row vector affected by batch effects

rbyx <- colMeans(R[21:40,]) # mean residual data row vector afffected by both batch and class effects
rbx <- colMeans(R[41:60,]) # mean residual data row vector affected by batch effects.

acos(abs(mean(cor(apply(R[21:40, ], 1, as.numeric), pc1)))) /pi * 180
# [1] 0.4559973
cor(pc1, rbyx)
# [1] -0.9999992
acos(abs(cor(pc1, rbyx))) /pi * 180
# [1] 0.0745133

acos(abs(mean(cor(apply(R[41:60, ], 1, as.numeric), pc1)))) /pi * 180
# [1] 0.4802948
cor(pc1, rbx)
# [1] 0.9999992
acos(abs(cor(pc1, rbx))) /pi * 180
# [1] 0.07444972

acos(abs(mean(cor(apply(Edata[21:40, ], 1, as.numeric), pc1)))) /pi * 180
# [1] 18.42746
cor(pc1, byx)
# [1] 0.9487513
acos(abs(cor(pc1, byx))) /pi * 180
# [1] 18.42263

acos(abs(mean(cor(apply(Edata[41:60, ], 1, as.numeric), pc1)))) /pi * 180
# [1] 53.12177
cor(pc1, bx)
# [1] 0.6001239
acos(abs(cor(pc1, bx))) /pi * 180
# [1] 53.12123

rowvecDf <- rbind(y, bv, sqrt(40)*sv, sqrt(40)*pc1, byx, bx, rbyx, rbx)
rownames(rowvecDf) <- c('class_vector','batch_vector','surrogate_variable','pc1',
                        'gene_affected_by_batch_class','gene_affected_by_batch',
                        'residual_gene_affected_by_batch_class','residual_gene_affected_by_batch')

Mrowvec <- melt(rowvecDf)
gg4.3 <- ggplot(data = Mrowvec, aes(x = Var2, y = value)) +
  facet_wrap(~Var1, ncol = 2) + 
  geom_point() +
  scale_y_continuous(limits = c(-5, 5), breaks=c(-4,-2,-1,0,1,2,4)) + 
  labs(x = 'sample', y = 'value', 
       title = paste0("Vectors in the Row Space (noise = 0.01)\n pi_1 = ", pi_1, 
                      ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

pdf('figures/mean_heterogeneity/angle/angle_analysis_pc_x/gg4.1.pdf')
print(gg4.1)
dev.off()
pdf('figures/mean_heterogeneity/angle/angle_analysis_pc_x/gg4.2.pdf')
print(gg4.2)
dev.off()
pdf('figures/mean_heterogeneity/angle/angle_analysis_pc_x/gg4.3.pdf')
print(gg4.3)
dev.off()

# case2: pi_1 = 0.9, pi_2 = 0.1 -------------------------------------------

Data5 = Generate_mean_hetero_data(pi_1 = 0.9, pi_2 = 0.1, noise = 0.01)
pi_1 <- unique(Data5$pi_1); pi_2 <- unique(Data5$pi_2)
p_1 <- (pi_1 + pi_2)/2; p_2 <- 1 - p_1

M <- Data5 %>% select(-y, -batch, -pi_1, -pi_2, -c1, -c2)

# Expression data
Edata <- t(M)
colnames(Edata) <- paste0("sample", c(1: dim(Edata)[2]))
Mdata <- melt(Edata)
gg5.1 <- ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value", 
       title = paste0("Simulated data (noise = 0.01) \n pi_1 = ", pi_1, 
                      ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)) +
  geom_tile() + 
  scale_fill_gradient2(limits=c(-10, 10)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

## full model 
mod <- model.matrix(~ y, data = Data5)
## null model
mod0 <- model.matrix(~ 1, data = Data5)

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
gg5.2 <- ggplot(MR, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",
       title = paste0("Residual Matrix After Removing Primary Effect (noise = 0.01) \n pi_1 = ", pi_1, 
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
# [1] -0.5999995
acos(abs(cor(pc1, bv))) / pi * 180
# [1] 53.13014

cor(sv, bv)
# [1] -0.1949481
acos(abs(cor(sv, bv))) / pi * 180
# [1] 78.75831

acos(abs(cor(sv, pc1))) / pi * 180
# [1] 87.79553

byx <- colMeans(Edata[21:40,]) # mean data row vector affected by both batch and class effects
bx <- colMeans(Edata[41:60,])  # mean data row vector affected by batch effects

rbyx <- colMeans(R[21:40,]) # mean residual data row vector afffected by both batch and class effects
rbx <- colMeans(R[41:60,]) # mean residual data row vector affected by batch effects.

acos(abs(mean(cor(apply(R[21:40, ], 1, as.numeric), pc1)))) /pi * 180
# [1] 0.4540579
cor(pc1, rbyx)
# [1] -0.9999992
acos(abs(cor(pc1, rbyx))) / pi * 180
# [1] 0.07392095
acos(abs(mean(cor(apply(R[41:60, ], 1, as.numeric), pc1)))) /pi * 180
# [1] 0.456861

acos(abs(mean(cor(apply(Edata[21:40, ], 1, as.numeric), pc1)))) /pi * 180
# [1] 71.56871
cor(pc1, byx)
# [1] -0.3161682
acos(abs(cor(pc1, byx))) / pi * 180
# [1] 71.56865

acos(abs(mean(cor(apply(Edata[41:60, ], 1, as.numeric), pc1)))) /pi * 180
# [1] 53.12085
cor(pc1, bx)
# [1] -0.6001359
acos(abs(cor(pc1, bx))) / pi * 180
# [1] 53.12037

rowvecDf <- rbind(y, bv, sqrt(40)*sv, sqrt(40)*pc1, byx, bx, rbyx, rbx)
rownames(rowvecDf) <- c('class_vector','batch_vector','surrogate_variable','pc1',
                        'gene_affected_by_batch_class','gene_affected_by_batch',
                        'residual_gene_affected_by_batch_class','residual_gene_affected_by_batch')

Mrowvec <- melt(rowvecDf)
gg5.3 <- ggplot(data = Mrowvec, aes(x = Var2, y = value)) +
  facet_wrap(~Var1, ncol = 2) + 
  geom_point() +
  labs(x = "Sample", y = "Value",
       title = paste0("Vectors in the Row Space (noise = 0.01) \n pi_1 = ", pi_1, 
                      ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)) +
  scale_y_continuous(limits = c(-5, 5), breaks=c(-4,-2,-1,0,1,2,4)) + 
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

pdf('figures/mean_heterogeneity/angle/angle_analysis_pc_x/gg5.1.pdf')
print(gg5.1)
dev.off()
pdf('figures/mean_heterogeneity/angle/angle_analysis_pc_x/gg5.2.pdf')
print(gg5.2)
dev.off()
pdf('figures/mean_heterogeneity/angle/angle_analysis_pc_x/gg5.3.pdf')
print(gg5.3)
dev.off()

# Case3: pi_1 = 0.5, pi_2 = 0.5 -------------------------------------------

Data6 = Generate_mean_hetero_data(pi_1 = 0.5, pi_2 = 0.5, noise = 0.01)
pi_1 <- unique(Data6$pi_1); pi_2 <- unique(Data6$pi_2)
p_1 <- (pi_1 + pi_2)/2; p_2 <- 1 - p_1

M <- Data6 %>% select(-y, -batch, -pi_1, -pi_2, -c1, -c2)

# Expression data
Edata <- t(M)
colnames(Edata) <- paste0("sample", c(1: dim(Edata)[2]))
Mdata <- melt(Edata)
gg6.1 <- ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value", 
       title = paste0("Simulated data (noise = 0.01) \n pi_1 = ", pi_1, 
                      ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)) +
  geom_tile() + 
  scale_fill_gradient2(limits=c(-10, 10)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

## full model 
mod <- model.matrix(~ y, data = Data6)
## null model
mod0 <- model.matrix(~ 1, data = Data6)

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
gg6.2 <- ggplot(MR, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value",
       title = paste0("Residual Matrix After Removing Primary Effect (noise = 0.01) \n pi_1 = ", pi_1, 
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
# [1] 0.9999997
acos(abs(cor(pc1, bv))) / pi * 180
# [1] 0.04534867
cor(sv, bv)
# [1] -0.9999993
acos(abs(cor(sv, bv))) / pi * 180
# [1] 0.0681486
acos(abs(cor(sv, pc1))) / pi * 180
# [1] 0.04371337

byx <- colMeans(Edata[21:40,]) # mean data row vector affected by both batch and class effects
bx <- colMeans(Edata[41:60,])  # mean data row vector affected by batch effects

rbyx <- colMeans(R[21:40,]) # mean residual data row vector afffected by both batch and class effects
rbx <- colMeans(R[41:60,]) # mean residual data row vector affected by batch effects.

acos(abs(mean(cor(apply(R[21:40, ], 1, as.numeric), pc1)))) /pi * 180
# [1] 0.2829739
cor(pc1, rbyx)
# [1] 0.9999998
acos(abs(cor(pc1, rbyx))) / pi * 180
# [1] 0.03938041

acos(abs(mean(cor(apply(R[41:60, ], 1, as.numeric), pc1)))) /pi * 180
# [1] 0.286226
cor(pc1, rbx)
# [1] 0.9999998
acos(abs(cor(pc1, rbx))) / pi * 180
# [1] 0.03937629

acos(abs(mean(cor(apply(Edata[21:40, ], 1, as.numeric), pc1)))) /pi * 180
# [1] 44.99945
cor(pc1, byx)
# [1] 0.7071178
acos(abs(cor(pc1, byx))) / pi * 180
# [1] 44.99911

acos(abs(mean(cor(apply(Edata[41:60, ], 1, as.numeric), pc1)))) /pi * 180
# [1] 0.2880784
cor(pc1, bx)
# [1] 0.9999998
acos(abs(cor(pc1, bx))) / pi * 180
# [1] 0.03938698

rowvecDf <- rbind(y, bv, sqrt(40)*sv, sqrt(40)*pc1, byx, bx, rbyx, rbx)
rownames(rowvecDf) <- c('class_vector','batch_vector','surrogate_variable','pc1',
                        'gene_affected_by_batch_class','gene_affected_by_batch',
                        'residual_gene_affected_by_batch_class','residual_gene_affected_by_batch')

Mrowvec <- melt(rowvecDf)
gg6.3 <- ggplot(data = Mrowvec, aes(x = Var2, y = value)) +
  facet_wrap(~Var1, ncol = 2) + 
  geom_point() +
  labs(x = "Sample", y = "Value",
       title = paste0("Vectors in the Row Space (noise = 0.01) \n pi_1 = ", pi_1, 
                      ", pi_2 = ", pi_2, "\n p_1 = ", p_1, ", p_2 = ", p_2)) +
  scale_y_continuous(limits = c(-5, 5), breaks=c(-4,-2,-1,0,1,2,4)) + 
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

pdf('figures/mean_heterogeneity/angle/angle_analysis_pc_x/gg6.1.pdf')
print(gg6.1)
dev.off()
pdf('figures/mean_heterogeneity/angle/angle_analysis_pc_x/gg6.2.pdf')
print(gg6.2)
dev.off()
pdf('figures/mean_heterogeneity/angle/angle_analysis_pc_x/gg6.3.pdf')
print(gg6.3)
dev.off()
