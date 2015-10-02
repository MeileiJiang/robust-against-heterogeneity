################################################################################
## Angle_Analysis.R
## Analyze the angle between batch vector and primary variable vector against 
## unbalanced setting pi_1 and pi_2.
## Author: Meilei Jiang
################################################################################
library(ggplot2)
library(reshape2)
library(dplyr)
library(sva)
source('R/mean_heterogeneity/fixed_effects_mean_hetero_sva/Generate_mean_hetero_data.R')
## parallel computing packages
library(doMC)
library(foreach)

n = 80; # Total number of samples
c1 = 40; # Total number of samples in the Class 1
c2 = 40; # Total number of samples in the Class 2
t = c(20, 20, 20, 40) # Settings of features

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
                                    "")))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red") +
  geom_text(x = 1, y = 1, label ="O", col = "black") +
  geom_text(x = 41, y = 41, label = "O", col = "black")

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
    tempdata <- Generate_mean_hetero_data(pi_1 = (i-1)/40, pi_2 = (j-1)/40)
    M <- tempdata %>% select(-y, -batch, -pi_1, -pi_2, -c1, -c2)
    # Expression data
    Edata <- t(M)
    colnames(Edata) <- paste0("sample", c(1: dim(Edata)[2]))
    
    ## full model 
    mod <- model.matrix(~ y, data = tempdata)
    ## null model
    mod0 <- model.matrix(~ 1, data = tempdata)
    ## estimate the number of surrogate variable
    n.sv <- num.sv(Edata, mod, method="be", B = 100)
    ## estimate the surrogate varialbe
    svobj <- sva(Edata, mod, mod0, n.sv= n.sv) 
    ## get the surrogate variable
    sv <- svobj$sv
    ## get the batch vector
<<<<<<< HEAD
    bv = 2 - as.numeric(unlist(tempdata %>% select(batch)))
    y = 2 - as.numeric(unlist(tempdata %>% select(y)))
=======
    bv = as.numeric(unlist(tempdata %>% select(batch)))
    y = as.numeric(unlist(tempdata %>% select(y)))
>>>>>>> 12e3d54d0bc482013305553c0f66297710e7b53e
    
    if(length(sv) == 1){
      cor_sv_bv <- NA
      cor_sv_y <- NA
    }else{
      if(sd(bv) == 0 | sd(sv) == 0){
        cor_sv_bv <- 0
      }else{
        cor_sv_bv <- max(abs(cor(sv,bv)))
      }
      
      cor_sv_y <- max(abs(cor(sv,y)))
    }
    data.frame(pi_1 = i, pi_2 = j, replicate = k, n.sv = n.sv, cor_sv_bv = cor_sv_bv, cor_sv_y = cor_sv_y)
  }
<<<<<<< HEAD

Cor_sv_df <- Cor_sv_df %>%
  mutate(diff_cor = cor_sv_bv - cor_sv_y)

g2 = ggplot(data = Cor_sv_df, aes(x = pi_1, y = pi_2, fill = acos(cor_sv_bv) /pi * 180)) + 
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:9)*10, limits=c(0, 90), high = "white", mid = "blue") +
  labs(x = expression(pi[1]), y = expression(pi[2]), fill = 'angle',
       title = expression(atop("Angle Between Surrogate Variable and Batch Vector", 
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    "")))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red") +
  geom_text(x = 1, y = 1, label ="O", col = "black") +
  geom_text(x = 41, y = 41, label = "O", col = "black")

g3 = ggplot(data = Cor_sv_df, aes(x = pi_1, y = pi_2, fill = acos(abs(cor_sv_bv - cor_sv_y)) /pi * 180)) + 
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:9)*10, limits=c(0, 90), high = "white", mid = "blue") +
  labs(x = expression(pi[1]), y = expression(pi[2]), fill = 'angle',
       title = expression(atop("Angle Between Surrogate Variable and Class Vector", 
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    "")))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red") +
  geom_text(x = 1, y = 1, label ="O", col = "black") +
  geom_text(x = 41, y = 41, label = "O", col = "black")

# Try to understand why the upper left corner is darker than the lower right corner.

mean_Cor_sv_df <- Cor_sv_df %>% 
  group_by(pi_1, pi_2) %>%
  summarise(mean_cor_sv_bv = mean(cor_sv_bv), mean_cor_sv_y = mean(cor_sv_y))

ggplot(data = mean_Cor_sv_df, aes(x = pi_1, y = pi_2, fill = acos(mean_cor_sv_bv) /pi * 180)) + 
=======

g2 = ggplot(data = Cor_sv_df, aes(x = pi_1, y = pi_2, fill = acos(cor_sv_bv) /pi * 180)) + 
>>>>>>> 12e3d54d0bc482013305553c0f66297710e7b53e
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:9)*10, limits=c(0, 90), high = "white", mid = "blue") +
  labs(x = expression(pi[1]), y = expression(pi[2]), fill = 'angle',
       title = expression(atop("Angle Between Surrogate Variable and Batch Vector", 
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    "")))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red") +
  geom_text(x = 1, y = 1, label ="O", col = "black") +
  geom_text(x = 41, y = 41, label = "O", col = "black")

<<<<<<< HEAD
ggplot(data = mean_Cor_sv_df, aes(x = pi_1, y = pi_2, fill = acos(mean_cor_sv_bv - mean_cor_sv_y) /pi * 180)) + 
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:9)*10, limits=c(0, 90), high = "white", mid = "blue") +
  labs(x = expression(pi[1]), y = expression(pi[2]), fill = 'angle',
       title = expression(atop("Angle Between Surrogate Variable and Class Vector", 
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    "")))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red") +
  geom_text(x = 1, y = 1, label ="O", col = "black") +
  geom_text(x = 41, y = 41, label = "O", col = "black")

ggplot(data = mean_Cor_sv_df, aes(x = pi_1, y = pi_2, fill = mean_cor_sv_bv)) + 
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:10)/10, limits=c(0, 1), high = "blue", mid = "white") +
  labs(x = expression(pi[1]), y = expression(pi[2]), fill = 'angle',
       title = expression(atop("Correlation Coefficient Between Surrogate Variable and Class Vector", 
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    "")))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red") +
  geom_text(x = 1, y = 1, label ="O", col = "black") +
  geom_text(x = 41, y = 41, label = "O", col = "black")

leftupper_Cor_sv_df <- Cor_sv_df %>% 
  filter(pi_1 <= 11, pi_2 >= 31)

rightlower_Cor_sv_df <- Cor_sv_df %>%
  filter(pi_1 >= 31, pi_2 <= 11)

t.test(leftupper_Cor_sv_df$cor_sv_bv, rightlower_Cor_sv_df$cor_sv_y)

ggplot(data = leftupper_Cor_sv_df, aes(x = pi_1, y = pi_2, fill = acos(cor_sv_bv) /pi * 180)) + 
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:9)*10, limits=c(0, 90), high = "white", mid = "blue") +
  labs(x = expression(pi[1]), y = expression(pi[2]), fill = 'angle',
       title = expression(atop("Angle Between Surrogate Variable and Batch Vector: Left Upper Corner", 
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    "")))) +
  scale_x_discrete(breaks = seq(from = 1, to = 10, by = 1),labels = c(0:9)/40) +
  scale_y_discrete(breaks = seq(from = 32, to = 41, by = 1),labels = c(31:40)/40) 


ggplot(data = Cor_sv_df, aes(x = pi_1, y = pi_2, fill = acos(cor_sv_bv) /pi * 180)) + 
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:9)*10, limits=c(0, 90), high = "white", mid = "blue") +
  labs(x = expression(pi[1]), y = expression(pi[2]), fill = 'angle',
       title = expression(atop("Angle Between Surrogate Variable and Batch Vector", 
=======
g3 = ggplot(data = Cor_sv_df, aes(x = pi_1, y = pi_2, fill = acos(cor_sv_y) /pi * 180)) + 
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:9)*10, limits=c(0, 90), high = "white", mid = "blue") +
  labs(x = expression(pi[1]), y = expression(pi[2]), fill = 'angle',
       title = expression(atop("Angle Between Surrogate Variable and Class Vector", 
>>>>>>> 12e3d54d0bc482013305553c0f66297710e7b53e
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    "")))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red") +
  geom_text(x = 1, y = 1, label ="O", col = "black") +
  geom_text(x = 41, y = 41, label = "O", col = "black")
# save the computation result

save(Cor_sv_df, file = 'rdata/angle_analysis2.RData')




# save the computation result

save(Cor_sv_df, file = 'rdata/angle_analysis2.RData')




# output the figures

pdf(file = 'figures/mean_heterogeneity/angle/angle_y_bv.pdf')
print(g1)
dev.off()

pdf(file = 'figures/mean_heterogeneity/angle/angle_sv_bv.pdf')
print(g2)
dev.off()

pdf(file = 'figures/mean_heterogeneity/angle/angle_y_sv.pdf')
print(g3)
dev.off()



