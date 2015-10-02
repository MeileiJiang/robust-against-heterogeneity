###################################################################################################
## angle_analysis2.R
## In order to understand the asymmetric of the angle matrix for SV and Batch Vector, I develop the
## angle between PC1 and batch vector and the angle between PC1 and data.
## Author: Meilei
###################################################################################################
library(ggplot2)
library(reshape2)
library(dplyr)


# f1 is the correlation coefficient between PC1 and Batch vector under different balanced settings
# This is derived theoretically.
f1 = function(pi_1, pi_2){
  if(pi_1 > 1 | pi_1 < 0) stop("pi_1 is out of range")
  if(pi_2 > 1 | pi_2 < 0) stop("pi_2 is out of range")

    if((pi_1 + pi_2 - pi_1^2 - pi_2^2) == 0){
      return(0)
  }
  return(sqrt(2 * (pi_1 + pi_2 - pi_1^2 - pi_2^2)/((pi_1 + pi_2)*(2 - pi_1 - pi_2))))
}

CorMat = matrix(ncol = 41, nrow = 41)
for(i in 1:41){
  for(j in 1:41){
    CorMat[i, j] = min(f1((i-1)/40, (j-1)/40),1)
  }
}

CorDf = melt(CorMat) %>%
  rename(pi_1 = Var1, pi_2 = Var2, Cor = value) %>%
  mutate(angle = acos(Cor)/pi * 180)

g1 = ggplot(data = CorDf, aes(x = pi_1, y = pi_2, fill = angle)) + 
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:9)*10, limits=c(0, 90), high = "white", mid = "blue") +
  labs(x= expression(pi[1]), y = expression(pi[2]), fill = 'angle',
       title = expression(atop("Angle Between PC1 and Batch Vector", 
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    italic("Green letters represent the case batch effect is confounded with class effect"))))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red", size = 2) +
  geom_text(x = 1, y = 1, label ="O", col = "black", size = 2) +
  geom_text(x = 41, y = 41, label = "O", col = "black", size = 2) +
  geom_text(x = 1, y = 41, label = "N", col = "green", size = 2) +
  geom_text(x = 41, y = 1, label = "P", col = "green", size = 2)

# f2 is the correlation coefficient between PC1 and Data under different balanced settings
# This is derived theoretically.
f2 = function(pi_1, pi_2){
  if(pi_1 > 1 | pi_1 < 0) stop("pi_1 is out of range")
  if(pi_2 > 1 | pi_2 < 0) stop("pi_2 is out of range")
  
  if((pi_1 + pi_2 - pi_1^2 - pi_2^2) == 0){
    return(0)
  }
  return(sqrt(2 * (pi_1 + pi_2 - pi_1^2 - pi_2^2)/(pi_1 + 1 - pi_2 -(1 - pi_1 - pi_2)^2/2)))
}

CorMat1 = matrix(ncol = 41, nrow = 41)
for(i in 1:41){
  for(j in 1:41){
    CorMat1[i, j] = min(f2((i-1)/40, (j-1)/40),1)
  }
}

CorDf1 = melt(CorMat1) %>%
  rename(pi_1 = Var1, pi_2 = Var2, Cor = value) %>%
  mutate(angle = acos(Cor)/pi * 180)

g2 = ggplot(data = CorDf1, aes(x = pi_1, y = pi_2, fill = angle)) + 
  geom_tile() +
  scale_fill_gradient2(breaks =c(0:9)*10, limits=c(0, 90), high = "white", mid = "blue") +
  labs(x= expression(pi[1]), y = expression(pi[2]), fill = 'angle',
       title = expression(atop("Angle Between PC1 and Gene with both batch effects and class effects", 
                               atop(italic("Red X represents the balanced case, Black O represents the none batch effect case"),
                                    italic("Green letters represent the case batch effect is confounded with class effect"))))) +
  scale_x_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  scale_y_discrete(breaks = seq(from = 1, to = 41, by = 10),labels = c(0:4)/4) +
  geom_text(x = 21, y = 21, label = "X", col = "red", size = 2) +
  geom_text(x = 1, y = 1, label ="O", col = "black", size = 2) +
  geom_text(x = 41, y = 41, label = "O", col = "black", size = 2) +
  geom_text(x = 1, y = 41, label = "N", col = "green", size = 2) +
  geom_text(x = 41, y = 1, label = "P", col = "green", size = 2)

pdf('figures/mean_heterogeneity/angle/angle_pc1_bv.pdf')
print(g1)
dev.off()

pdf('figures/mean_heterogeneity/angle/angle_pc1_x.pdf')
print(g2)
dev.off()

