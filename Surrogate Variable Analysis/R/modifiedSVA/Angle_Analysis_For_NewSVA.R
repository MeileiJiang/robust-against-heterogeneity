######################################################################################
## Angle_Analysis_For_NewSVA.R
## Analyze the angle between batch vector and primary variable under new gene settings
## 
## Author: Meilei Jiang
######################################################################################
library(ggplot2)
library(reshape2)
library(dplyr)
library(sva)
source('R/mean_heterogeneity/fixed_effects_mean_hetero_sva/Generate_mean_hetero_data.R')
source('R/modifiedSVA/helper.R')
source('R/modifiedSVA/newirwsva.R')

## parallel computing packages
library(doMC)
library(foreach)

n = 80; # Total number of samples
c1 = 40; # Total number of samples in the Class 1
c2 = 40; # Total number of samples in the Class 2
t = c(0, 50, 0, 50) # Settings of features

# Angle between batch vector and primary variable vector
Angle_bv_y = function(pi_1, pi_2){
  n11 = floor(c1 * pi_1); # Total number of samples in the Class 1 from Batch 1.
  n12 = c1 - n11; # Total number of samples in the Class 1 from Batch 2.
  n21 = floor(c2 * pi_2) # Total number of samples in the Class 2 from Batch 1.
  n22 = c2 - n21 # Total number of samples in the Class 2 from Batch 2.
  y = rep(c(1,1,0,0), c(n11, n12, n21, n22))
  bv = rep(c(1,0,1,0), c(n11, n12, n21, n22))
  acos(abs(cor(y,bv)))/pi * 180 
}

AngleMat = matrix(ncol = 41, nrow = 41)
for(i in 1:41){
  for(j in 1:41){
    AngleMat[i, j] = Angle_bv_y((i-1)/40, (j-1)/40)
  }
}

AngleMat[1, 1] = 90
AngleMat[41, 41] = 90


AngleDf = melt(AngleMat) %>%
  rename(pi_1 = Var1, pi_2 = Var2, angle = value)

g1 = ggplot(data = AngleDf, aes(x = pi_1, y = pi_2, fill = angle)) + 
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:9)*10, limits=c(0, 90), high = "white", mid = "blue") +
  labs(x= expression(pi[1]), y = expression(pi[2]),
       title = expression(atop("Angle Between Class Vector and Batch Vector", 
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    italic("Green letters represent the case batch effect is confounded with class effect"))))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red", size = 2) +
  geom_text(x = 1, y = 1, label ="O", col = "black", size = 2) +
  geom_text(x = 41, y = 41, label = "O", col = "black", size = 2) +
  geom_text(x = 1, y = 41, label = "N", col = "green", size = 2) +
  geom_text(x = 41, y = 1, label = "P", col = "green", size = 2)

# Angle between surrogate variable and batch vector
## To overcome the randomness in the simulation data, each case has been repeated 20 times


# parallel computing

## set up the index

n1 = 41 # number of pi_1
n2 = 41 # number of pi_2
indexdf = data.frame(outer=rep(1:n1,rep(n2,n1)),inner=rep(1:n2,n1))
index = 1:nrow(indexdf)

## set up parallel computing

registerDoMC(10)

Cor_sv_df = foreach(l = index, .combine = 'rbind') %:%
  foreach(k = 1:20, .combine = 'rbind') %dopar%{
    i <- indexdf[l,1]
    j <- indexdf[l,2]
    # Generate new data set 
    tempdata <- Generate_mean_hetero_data(pi_1 = (i-1)/40, pi_2 = (j-1)/40, t = t)
    M <- tempdata %>% select(-y, -batch, -pi_1, -pi_2, -c1, -c2)
    # Expression data
    Edata <- t(M)
    colnames(Edata) <- paste0("sample", c(1: dim(Edata)[2]))
    
    ## full model 
    mod <- model.matrix(~ y, data = tempdata)
    ## null model
    mod0 <- model.matrix(~ 1, data = tempdata)
    ## consturct the hat matrix and get residual matrix R
    HatB = Edata %*% mod %*% solve(t(mod) %*% mod) 
    R = Edata - HatB %*% t(mod)
    # get the PC1 from R
    pc1 = svd(R)$v[,1]
    
    ## get the batch vector and class vector
    bv = 2 - as.numeric(unlist(tempdata %>% select(batch)))
    y = 2 - as.numeric(unlist(tempdata %>% select(y)))
    if(sd(bv) == 0){
      cor_bv_pc1 <- 0
    }else{
      cor_bv_pc1 <- cor(bv, pc1)
    }
    cor_y_pc1 <- cor(y, pc1)
    cor_x_pc1 <- mean(abs(cor(apply(Edata[21:40, ], 1, as.numeric), pc1)))
    
    ## estimate the number of surrogate variable
    n.sv <- num.sv(Edata, mod, method="be", B = 100)
    ## estimate the surrogate varialbe
    svobj <- newirwsva(Edata, mod, mod0, n.sv= n.sv) 
    ## get the surrogate variable
    sv <- svobj$sv
    
    if(length(sv) == 1){
      cor_sv_bv <- NA
      cor_sv_y <- NA
      cor_sv_pc1 <- NA
    }else{
      if(sd(bv) == 0 | sd(sv) == 0){
        cor_sv_bv <- 0
      }else{
        cor_sv_bv <- max(abs(cor(sv,bv)))
      }
      cor_sv_y <- max(abs(cor(sv,y)))
      cor_sv_pc1 <- max(abs(cor(sv, pc1)))
    }
    data.frame(pi_1 = i, pi_2 = j, replicate = k, n.sv = n.sv, cor_sv_bv = cor_sv_bv, cor_sv_y = cor_sv_y, 
               cor_sv_pc1 = cor_sv_pc1, cor_bv_pc1 = cor_bv_pc1, cor_y_pc1 = cor_y_pc1, cor_x_pc1 = cor_x_pc1 )
  }

save(Cor_sv_df, file = "")

g2 = ggplot(data = Cor_sv_df, aes(x = pi_1, y = pi_2, fill = acos(cor_sv_bv) /pi * 180)) + 
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:9)*10, limits=c(0, 90), high = "white", mid = "blue") +
  labs(x = expression(pi[1]), y = expression(pi[2]), fill = 'angle',
       title = expression(atop("Angle Between Surrogate Variable and Batch Vector", 
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    italic("Green letters represent the case batch effect is confounded with class effect"))))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red", size = 2) +
  geom_text(x = 1, y = 1, label ="O", col = "black", size = 2) +
  geom_text(x = 41, y = 41, label = "O", col = "black", size = 2) +
  geom_text(x = 1, y = 41, label = "N", col = "green", size = 2) +
  geom_text(x = 41, y = 1, label = "P", col = "green", size = 2)

g3 <- ggplot(data = Cor_sv_df, aes(x = pi_1, y = pi_2, fill = acos(abs(cor_bv_pc1)) /pi * 180)) + 
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:9)*10, limits=c(0, 90), high = "white", mid = "blue") +
  labs(x = expression(pi[1]), y = expression(pi[2]), fill = 'angle',
       title = expression(atop("Angle Between Batch Vector and PC1", 
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    italic("Green letters represent the case batch effect is confounded with class effect"))))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red", size = 2) +
  geom_text(x = 1, y = 1, label ="O", col = "black", size = 2) +
  geom_text(x = 41, y = 41, label = "O", col = "black", size = 2) +
  geom_text(x = 1, y = 41, label = "N", col = "green", size = 2) +
  geom_text(x = 41, y = 1, label = "P", col = "green", size = 2)

g4 <- ggplot(data = Cor_sv_df, aes(x = pi_1, y = pi_2, fill = acos(cor_sv_pc1) /pi * 180)) + 
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:9)*10, limits=c(0, 90), high = "white", mid = "blue") +
  labs(x = expression(pi[1]), y = expression(pi[2]), fill = 'angle',
       title = expression(atop("Angle Between Surrogate Variable and PC1", 
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    italic("Green letters represent the case batch effect is confounded with class effect"))))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red", size = 2) +
  geom_text(x = 1, y = 1, label ="O", col = "black", size = 2) +
  geom_text(x = 41, y = 41, label = "O", col = "black", size = 2) +
  geom_text(x = 1, y = 41, label = "N", col = "green", size = 2) +
  geom_text(x = 41, y = 1, label = "P", col = "green", size = 2)

g5 <- ggplot(data = Cor_sv_df, aes(x = pi_1, y = pi_2, fill = acos(cor_x_pc1) /pi * 180)) + 
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:9)*10, limits=c(0, 90), high = "white", mid = "blue") +
  labs(x = expression(pi[1]), y = expression(pi[2]), fill = 'angle',
       title = expression(atop("Angle Between PC1 and Gene Affected by Both Batch & Class", 
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    italic("Green letters represent the case batch effect is confounded with class effect"))))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red", size = 2) +
  geom_text(x = 1, y = 1, label ="O", col = "black", size = 2) +
  geom_text(x = 41, y = 41, label = "O", col = "black", size = 2) +
  geom_text(x = 1, y = 41, label = "N", col = "green", size = 2) +
  geom_text(x = 41, y = 1, label = "P", col = "green", size = 2)

g6 <- ggplot(data = Cor_sv_df, aes(x = pi_1, y = pi_2, fill = acos(cor_y_pc1) /pi * 180)) + 
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:9)*10, limits=c(0, 90), high = "white", mid = "blue") +
  labs(x = expression(pi[1]), y = expression(pi[2]), fill = 'angle',
       title = expression(atop("Angle Between Class Vectors and PC1", 
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    italic("Green letters represent the case batch effect is confounded with class effect"))))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red", size = 2) +
  geom_text(x = 1, y = 1, label ="O", col = "black", size = 2) +
  geom_text(x = 41, y = 41, label = "O", col = "black", size = 2) +
  geom_text(x = 1, y = 41, label = "N", col = "green", size = 2) +
  geom_text(x = 41, y = 1, label = "P", col = "green", size = 2)

g7 <- ggplot(data = Cor_sv_df, aes(x = pi_1, y = pi_2, fill = acos(cor_sv_y) /pi * 180)) + 
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:9)*10, limits=c(0, 90), high = "white", mid = "blue") +
  labs(x = expression(pi[1]), y = expression(pi[2]), fill = 'angle',
       title = expression(atop("Angle Between Surrogate Variable and Class Vectors", 
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    italic("Green letters represent the case batch effect is confounded with class effect"))))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red", size = 2) +
  geom_text(x = 1, y = 1, label ="O", col = "black", size = 2) +
  geom_text(x = 41, y = 41, label = "O", col = "black", size = 2) +
  geom_text(x = 1, y = 41, label = "N", col = "green", size = 2) +
  geom_text(x = 41, y = 1, label = "P", col = "green", size = 2)


g8 <- ggplot(data = Cor_sv_df, 
             aes(x = pi_1, y = pi_2, 
                 fill = (acos(cor_sv_bv) - acos(abs(cor_bv_pc1)) )/pi * 180)) +
  geom_tile() + 
  scale_fill_gradient2(breaks =c(-9:9)*10, limits=c(-90, 90),
                       high = "red", mid = 'white', low = "blue") +
  labs(x = expression(pi[1]), y = expression(pi[2]), fill = 'angle',
       title = expression(atop("Angle Difference Between Surrogate Variable And PC1 With Batch Vector", 
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    italic("Green letters represent the case batch effect is confounded with class effect"))))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red", size = 2) +
  geom_text(x = 1, y = 1, label ="O", col = "black", size = 2) +
  geom_text(x = 41, y = 41, label = "O", col = "black", size = 2) +
  geom_text(x = 1, y = 41, label = "N", col = "green", size = 2) +
  geom_text(x = 41, y = 1, label = "P", col = "green", size = 2)




# Try to understand why the upper left corner is darker than the lower right corner.

mean_Cor_sv_df <- Cor_sv_df %>% 
  group_by(pi_1, pi_2) %>%
  summarise(mean_cor_sv_bv = mean(cor_sv_bv), mean_cor_sv_y = mean(cor_sv_y), mean_cor_sv_pc1 = mean(cor_sv_pc1),
            mean_cor_y_pc1 = mean(abs(cor_y_pc1)), mean_cor_x_pc1 = mean(abs(cor_x_pc1)), mean_diff_cor = mean(diff_cor))

ggplot(data = mean_Cor_sv_df, aes(x = pi_1, y = pi_2, fill = acos(mean_cor_sv_bv) /pi * 180)) + 
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:9)*10, limits=c(0, 90), high = "white", mid = "blue") +
  labs(x = expression(pi[1]), y = expression(pi[2]), fill = 'angle',
       title = expression(atop("Angle Between Surrogate Variable and Batch Vector", 
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    italic("Green letters represent the case batch effect is confounded with class effect"))))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red", size = 2) +
  geom_text(x = 1, y = 1, label ="O", col = "black", size = 2) +
  geom_text(x = 41, y = 41, label = "O", col = "black", size = 2) +
  geom_text(x = 1, y = 41, label = "N", col = "green", size = 2) +
  geom_text(x = 41, y = 1, label = "P", col = "green", size = 2)

ggplot(data = mean_Cor_sv_df, aes(x = pi_1, y = pi_2, fill = acos(abs(mean_diff_cor)) /pi * 180)) + 
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:9)*10, limits=c(0, 90), high = "white", mid = "blue") +
  labs(x = expression(pi[1]), y = expression(pi[2]), fill = 'angle',
       title = expression(atop("Angle Between Surrogate Variable and Class Vector", 
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    italic("Green letters represent the case batch effect is confounded with class effect"))))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red", size = 2) +
  geom_text(x = 1, y = 1, label ="O", col = "black", size = 2) +
  geom_text(x = 41, y = 41, label = "O", col = "black", size = 2) +
  geom_text(x = 1, y = 41, label = "N", col = "green", size = 2) +
  geom_text(x = 41, y = 1, label = "P", col = "green", size = 2)


ggplot(data = mean_Cor_sv_df, aes(x = pi_1, y = pi_2, fill = acos(mean_cor_sv_pc1) /pi * 180)) + 
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:9)*10, limits=c(0, 90), high = "white", mid = "blue") +
  labs(x = expression(pi[1]), y = expression(pi[2]), fill = 'angle',
       title = expression(atop("Angle Between Surrogate Variable and PC1", 
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    italic("Green letters represent the case batch effect is confounded with class effect"))))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red", size = 2) +
  geom_text(x = 1, y = 1, label ="O", col = "black", size = 2) +
  geom_text(x = 41, y = 41, label = "O", col = "black", size = 2) +
  geom_text(x = 1, y = 41, label = "N", col = "green", size = 2) +
  geom_text(x = 41, y = 1, label = "P", col = "green", size = 2)


ggplot(data = mean_Cor_sv_df, aes(x = pi_1, y = pi_2, fill = acos(mean_cor_x_pc1) /pi * 180)) + 
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:9)*10, limits=c(0, 90), high = "white", mid = "blue") +
  labs(x = expression(pi[1]), y = expression(pi[2]), fill = 'angle',
       title = expression(atop("Angle Between Gene Affected By Both Batch and Class and PC1", 
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    italic("Green letters represent the case batch effect is confounded with class effect"))))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red", size = 2) +
  geom_text(x = 1, y = 1, label ="O", col = "black", size = 2) +
  geom_text(x = 41, y = 41, label = "O", col = "black", size = 2) +
  geom_text(x = 1, y = 41, label = "N", col = "green", size = 2) +
  geom_text(x = 41, y = 1, label = "P", col = "green", size = 2)

ggplot(data = mean_Cor_sv_df, aes(x = pi_1, y = pi_2, fill = acos(mean_cor_y_pc1) /pi * 180)) + 
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:9)*10, limits=c(0, 90), high = "white", mid = "blue") +
  labs(x = expression(pi[1]), y = expression(pi[2]), fill = 'angle',
       title = expression(atop("Angle Between Class Vectors and PC1", 
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    italic("Green letters represent the case batch effect is confounded with class effect"))))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red", size = 2) +
  geom_text(x = 1, y = 1, label ="O", col = "black", size = 2) +
  geom_text(x = 41, y = 41, label = "O", col = "black", size = 2) +
  geom_text(x = 1, y = 41, label = "N", col = "green", size = 2) +
  geom_text(x = 41, y = 1, label = "P", col = "green", size = 2)



# save the computation result

save(Cor_sv_df, file = 'rdata/angle_analysis_newirwsva.RData')




# output the figures

pdf(file = 'figures/mean_heterogeneity/angle/angle_y_bv.pdf')
print(g1)
dev.off()

pdf(file = 'figures/mean_heterogeneity/angle/angle_sv_bv_3.pdf')
print(g2)
dev.off()

pdf(file = 'figures/mean_heterogeneity/angle/angle_bv_pc1_3.pdf')
print(g3)
dev.off()

pdf(file = 'figures/mean_heterogeneity/angle/angle_sv_pc1.pdf')
print(g4)
dev.off()

pdf(file = 'figures/mean_heterogeneity/angle/angle_x_pc1.pdf')
print(g5)
dev.off()

pdf(file = 'figures/mean_heterogeneity/angle/angle_y_pc1.pdf')
print(g6)
dev.off()

pdf(file = 'figures/mean_heterogeneity/angle/angle_sv_y.pdf')
print(g7)
dev.off()

pdf(file = 'figures/mean_heterogeneity/angle/comp_sv_pc1.pdf')
print(g8)
dev.off()

