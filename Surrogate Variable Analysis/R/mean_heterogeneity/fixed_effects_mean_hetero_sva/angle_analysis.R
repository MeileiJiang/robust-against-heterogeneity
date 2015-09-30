################################################################################
## Angle_Analysis.R
## Analyze the angle between batch vector and primary variable vector against 
## unbalanced setting pi_1 and pi_2.
## Author: Meilei Jiang
################################################################################
library(ggplot2)
library(dplyr)
library(sva)
source('R/mean_heterogeneity/fixed_effects_mean_hetero_sva/Generate_mean_hetero_data.R')
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

AngleMat_sv <- matrix(ncol = 41, nrow = 41)
AngleMat_y <- matrix(ncol = 41, nrow = 41)

Zero_count <- matrix(0, ncol = 41, nrow = 41)
for(i in 1:41){
  for(j in 1:41){
    tempcor_sv <- rep(NA, 20)
    tempcor_y <- rep(NA, 20)
    for(k in 1:20){
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
      bv = as.numeric(unlist(tempdata %>% select(batch)))
      y = as.numeric(unlist(tempdata %>% select(y)))
      
      if(length(sv) == 1){
        Zero_count[i,j] <- Zero_count[i,j] + 1
      } else{
        if(sd(bv) == 0 | sd(sv) == 0){
          tempcor_sv[k] <- 0
        }else{
          tempcor_sv[k] <- max(cor(sv,bv))
        }
        
        tempcor_y[k] <- max(cor(sv,y))
      }
    }
    AngleMat_sv[i,j] <- acos(abs(mean(tempcor_sv, na.rm = TRUE))) /pi * 180
    AngleMat_y[i,j] <- acos(abs(mean(tempcor_y, na.rm = TRUE))) /pi * 180
  }
}

# save the computation result

save(AngleMat_y, AngleMat_sv, Zero_count, file = 'rdata/angle_analysis.RData')

# visulize the result

Angle_sv_Df = melt(AngleMat_sv) %>%
  rename(pi_1 = Var1, pi_2 = Var2, angle = value)

g2 = ggplot(data = Angle_sv_Df, aes(x = pi_1, y = pi_2, fill = angle)) + 
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:9)*10, limits=c(0, 90), high = "white", mid = "blue") +
  labs(x = expression(pi[1]), y = expression(pi[2]), 
       title = expression(atop("Angle Between Surrogate Variable and Batch Vector", 
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    "")))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red") +
  geom_text(x = 1, y = 1, label ="O", col = "black") +
  geom_text(x = 41, y = 41, label = "O", col = "black")


Angle_y_Df = melt(AngleMat_y) %>%
  rename(pi_1 = Var1, pi_2 = Var2, angle = value)

g3 = ggplot(data = Angle_y_Df, aes(x = pi_1, y = pi_2, fill = angle)) + 
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:9)*10, limits=c(0, 90), high = "white", mid = "blue") +
  labs(x = expression(pi[1]), y = expression(pi[2]), 
       title = expression(atop("Angle Between Surrogate Variable and Class Vector", 
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    "")))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red") +
  geom_text(x = 1, y = 1, label ="O", col = "black") +
  geom_text(x = 41, y = 41, label = "O", col = "black")

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



