##############################################################################################################
## Marg_Dist_Plot_Grun_Spikein.R
## Make the marginal distribution for Spikein Matrix of Grun reads cutoff dataset. 
## Author: Meilei Jiang
## Date: April 2016
##############################################################################################################

library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)
library(scales)
library(RUnit)
library(fpc)
library(moments)
library(gridExtra)
library(scales)

load(file = "rdata/Grun_reads_cutoff.RData")
load(file = "rdata/label.RData")

source("~/workspace/bdtbio_pcafuns/R/pcaSourceList.R")
source("R/utility/mytheme.R")


# data processing ---------------------------------------------------------

dim(Spikein)
# [1]  251 17
checkEquals(rownames(Spikein), rownames(Spikein)) 

Spikein1 = data.frame(Spikein) 
dim(Spikein1)
# [1]  251 17
# Rows are cells and columns are spikeins.


mSpikein =  melt(t(Spikein1)) # long format of Spikein matrix
colnames(mSpikein) = c("spikein","cell","expression")
mSpikein2 <- mSpikein%>% left_join(label.df)
mSpikein2$cell = factor(mSpikein$cell, levels = rownames(Spikein1))


# Marginal Distribution Plot ----------------------------------------------

spikein_mdplot = function(input, varname = NULL){
  inputspikein = as.character(input$spikein) 
  subSpikein = mSpikein2 %>% filter(spikein == inputspikein) 
  title = ""
  if(!is.na(varname)){
    value = round(input[,varname],5)
    title = paste(varname,"=",value)
  }
  
  
  g0 = ggplot(data = subSpikein, aes(expression)) + 
    geom_density(aes(group = label, col  = label)) +
    geom_jitter(aes(x =expression, y = 1.8, col = label), height = 1, size = 1) +
    labs(x = inputspikein, title  = title) 
  return(g0)
}

mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

pos = c(1:dim(Spikein1)[2])

spikeinMean = colMeans(Spikein1)
spikeinStd = apply(Spikein1, 2, sd)
spikeinSkew = apply(Spikein1, 2, skewness)
spikeinKurtosis = apply(Spikein1, 2, kurtosis)
spikeinMode = apply(Spikein1, 2, mode)


spikein.df = data.frame(spikein = names(spikeinMean), spikein_mean = spikeinMean, 
                     spikein_sd = spikeinStd, spikein_skewness = spikeinSkew, 
                     spikein_kurtosis = spikeinKurtosis, spikein_mode = spikeinMode)

plot(sort(spikein.df$spikein_mean))
plot(sort(spikein.df$spikein_sd))
plot(sort(spikein.df$spikein_skewness))
plot(sort(spikein.df$spikein_kurtosis))
plot(sort(spikein.df$spikein_mode))



# Standard deviation

spikein.df.sd <- spikein.df %>% arrange(spikein_sd)
spikein.df.sd$index <- c(1:length(spikeinMean))

sd.spikein = spikein.df.sd[pos,]

gsd = ggplot(spikein.df.sd, aes(x = index, y = spikein_sd)) + 
  geom_point() +
  geom_vline(xintercept = pos, lty = "dashed") +
  labs( x = "spikein", y = "standard deviation" ) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())


pdf("figures/grun_spikein_sd_mdplot.pdf", width = 10, height = 12)
grid.arrange(gsd,                                          spikein_mdplot(sd.spikein[1,],"spikein_sd"),  spikein_mdplot(sd.spikein[2,],"spikein_sd"),  
             spikein_mdplot(sd.spikein[3,],"spikein_sd"),  spikein_mdplot(sd.spikein[4,],"spikein_sd"),  spikein_mdplot(sd.spikein[5,],"spikein_sd"),  
             spikein_mdplot(sd.spikein[6,],"spikein_sd"),  spikein_mdplot(sd.spikein[7,],"spikein_sd"),  spikein_mdplot(sd.spikein[8,],"spikein_sd"), 
             spikein_mdplot(sd.spikein[9,],"spikein_sd"),  spikein_mdplot(sd.spikein[10,],"spikein_sd"), spikein_mdplot(sd.spikein[11,],"spikein_sd"), 
             spikein_mdplot(sd.spikein[12,],"spikein_sd"), spikein_mdplot(sd.spikein[13,],"spikein_sd"), spikein_mdplot(sd.spikein[14,],"spikein_sd"), 
             spikein_mdplot(sd.spikein[15,],"spikein_sd"), spikein_mdplot(sd.spikein[16,],"spikein_sd"), spikein_mdplot(sd.spikein[17,],"spikein_sd"),
             ncol = 3)
dev.off()


# Skewness

spikein.df.skewness <- spikein.df %>% arrange(spikein_skewness)
spikein.df.skewness$index <- c(1:length(spikeinSkew))
pos.skew = c(1:17)
skewness.spikein =spikein.df.skewness[pos.skew,]

gskew = ggplot(spikein.df.skewness, aes(x = index, y = spikein_skewness)) + 
  geom_point() +
  geom_vline(xintercept = pos.skew, lty = "dashed") +
  labs( x = "spikein", y = "skewness" ) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())


pdf("figures/grun_spikein_skew_mdplot.pdf", width = 10, height = 12)
grid.arrange(gskew,                                                    spikein_mdplot(skewness.spikein[1,],"spikein_skewness"),  spikein_mdplot(skewness.spikein[2,],"spikein_skewness"),  
             spikein_mdplot(skewness.spikein[3,],"spikein_skewness"),  spikein_mdplot(skewness.spikein[4,],"spikein_skewness"),  spikein_mdplot(skewness.spikein[5,],"spikein_skewness"),  
             spikein_mdplot(skewness.spikein[6,],"spikein_skewness"),  spikein_mdplot(skewness.spikein[7,],"spikein_skewness"),  spikein_mdplot(skewness.spikein[8,],"spikein_skewness"), 
             spikein_mdplot(skewness.spikein[9,],"spikein_skewness"),  spikein_mdplot(skewness.spikein[10,],"spikein_skewness"), spikein_mdplot(skewness.spikein[11,],"spikein_skewness"), 
             spikein_mdplot(skewness.spikein[12,],"spikein_skewness"), spikein_mdplot(skewness.spikein[13,],"spikein_skewness"), spikein_mdplot(skewness.spikein[14,],"spikein_skewness"), 
             spikein_mdplot(skewness.spikein[15,],"spikein_skewness"), spikein_mdplot(skewness.spikein[16,],"spikein_skewness"), spikein_mdplot(skewness.spikein[17,],"spikein_skewness"),
             ncol = 3)
dev.off()


# Kurtosis

spikein.df.kurtosis <- spikein.df %>% arrange(spikein_kurtosis)
spikein.df.kurtosis$index <- c(1:length(spikeinKurtosis))
pos.kurt = c(1:17)
kurtosis.spikein =spikein.df.kurtosis[pos.kurt,]

gkurt = ggplot(spikein.df.kurtosis, aes(x = index, y = spikein_kurtosis)) + 
  geom_point() +
  geom_vline(xintercept = pos.kurt, lty = "dashed") +
  labs( x = "spikein", y = "kurtosis" ) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

pdf("figures/grun_spikein_kurt_mdplot.pdf", width = 10, height = 12)
grid.arrange(gkurt,                                                    spikein_mdplot(kurtosis.spikein[1,],"spikein_kurtosis"),  spikein_mdplot(kurtosis.spikein[2,],"spikein_kurtosis"),  
             spikein_mdplot(kurtosis.spikein[3,],"spikein_kurtosis"),  spikein_mdplot(kurtosis.spikein[4,],"spikein_kurtosis"),  spikein_mdplot(kurtosis.spikein[5,],"spikein_kurtosis"),  
             spikein_mdplot(kurtosis.spikein[6,],"spikein_kurtosis"),  spikein_mdplot(kurtosis.spikein[7,],"spikein_kurtosis"),  spikein_mdplot(kurtosis.spikein[8,],"spikein_kurtosis"), 
             spikein_mdplot(kurtosis.spikein[9,],"spikein_kurtosis"),  spikein_mdplot(kurtosis.spikein[10,],"spikein_kurtosis"), spikein_mdplot(kurtosis.spikein[11,],"spikein_kurtosis"), 
             spikein_mdplot(kurtosis.spikein[12,],"spikein_kurtosis"), spikein_mdplot(kurtosis.spikein[13,],"spikein_kurtosis"), spikein_mdplot(kurtosis.spikein[14,],"spikein_kurtosis"), 
             spikein_mdplot(kurtosis.spikein[15,],"spikein_kurtosis"), spikein_mdplot(kurtosis.spikein[16,],"spikein_kurtosis"), spikein_mdplot(kurtosis.spikein[17,],"spikein_kurtosis"),
             ncol = 3)
dev.off()


# Mode

spikein.df.mode <- spikein.df %>% arrange(spikein_mode)
spikein.df.mode$index <- c(1:length(spikeinMode))
gmode = ggplot(spikein.df.mode, aes(x = index, y = spikein_mode)) + 
  geom_point() +
  geom_vline(xintercept = pos, lty = "dashed") +
  labs( x = "spikein", y = "mode" ) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

mode.spikein =spikein.df.mode[pos,]
pdf("figures/grun_spikein_mode_mdplot.pdf", width = 10, height = 12)
grid.arrange(gmode,                                                    spikein_mdplot(mode.spikein[1,],"spikein_mode"),  spikein_mdplot(mode.spikein[2,],"spikein_mode"),  
             spikein_mdplot(mode.spikein[3,],"spikein_mode"),  spikein_mdplot(mode.spikein[4,],"spikein_mode"),  spikein_mdplot(mode.spikein[5,],"spikein_mode"),  
             spikein_mdplot(mode.spikein[6,],"spikein_mode"),  spikein_mdplot(mode.spikein[7,],"spikein_mode"),  spikein_mdplot(mode.spikein[8,],"spikein_mode"), 
             spikein_mdplot(mode.spikein[9,],"spikein_mode"),  spikein_mdplot(mode.spikein[10,],"spikein_mode"), spikein_mdplot(mode.spikein[11,],"spikein_mode"), 
             spikein_mdplot(mode.spikein[12,],"spikein_mode"), spikein_mdplot(mode.spikein[13,],"spikein_mode"), spikein_mdplot(mode.spikein[14,],"spikein_mode"), 
             spikein_mdplot(mode.spikein[15,],"spikein_mode"), spikein_mdplot(mode.spikein[16,],"spikein_mode"), spikein_mdplot(mode.spikein[17,],"spikein_mode"),
             ncol = 3)
dev.off()


# combine together --------------------------------------------------------

pdf("figures/grun_spikein_mdplot.pdf", width = 10, height = 12)
grid.arrange(gsd,                                          spikein_mdplot(sd.spikein[1,],"spikein_sd"),  spikein_mdplot(sd.spikein[2,],"spikein_sd"),  
             spikein_mdplot(sd.spikein[3,],"spikein_sd"),  spikein_mdplot(sd.spikein[4,],"spikein_sd"),  spikein_mdplot(sd.spikein[5,],"spikein_sd"),  
             spikein_mdplot(sd.spikein[6,],"spikein_sd"),  spikein_mdplot(sd.spikein[7,],"spikein_sd"),  spikein_mdplot(sd.spikein[8,],"spikein_sd"), 
             spikein_mdplot(sd.spikein[9,],"spikein_sd"),  spikein_mdplot(sd.spikein[10,],"spikein_sd"), spikein_mdplot(sd.spikein[11,],"spikein_sd"), 
             spikein_mdplot(sd.spikein[12,],"spikein_sd"), spikein_mdplot(sd.spikein[13,],"spikein_sd"), spikein_mdplot(sd.spikein[14,],"spikein_sd"), 
             spikein_mdplot(sd.spikein[15,],"spikein_sd"), spikein_mdplot(sd.spikein[16,],"spikein_sd"), spikein_mdplot(sd.spikein[17,],"spikein_sd"),
             ncol = 3)
grid.arrange(gskew,                                                    spikein_mdplot(skewness.spikein[1,],"spikein_skewness"),  spikein_mdplot(skewness.spikein[2,],"spikein_skewness"),  
             spikein_mdplot(skewness.spikein[3,],"spikein_skewness"),  spikein_mdplot(skewness.spikein[4,],"spikein_skewness"),  spikein_mdplot(skewness.spikein[5,],"spikein_skewness"),  
             spikein_mdplot(skewness.spikein[6,],"spikein_skewness"),  spikein_mdplot(skewness.spikein[7,],"spikein_skewness"),  spikein_mdplot(skewness.spikein[8,],"spikein_skewness"), 
             spikein_mdplot(skewness.spikein[9,],"spikein_skewness"),  spikein_mdplot(skewness.spikein[10,],"spikein_skewness"), spikein_mdplot(skewness.spikein[11,],"spikein_skewness"), 
             spikein_mdplot(skewness.spikein[12,],"spikein_skewness"), spikein_mdplot(skewness.spikein[13,],"spikein_skewness"), spikein_mdplot(skewness.spikein[14,],"spikein_skewness"), 
             spikein_mdplot(skewness.spikein[15,],"spikein_skewness"), spikein_mdplot(skewness.spikein[16,],"spikein_skewness"), spikein_mdplot(skewness.spikein[17,],"spikein_skewness"),
             ncol = 3)
grid.arrange(gkurt,                                                    spikein_mdplot(kurtosis.spikein[1,],"spikein_kurtosis"),  spikein_mdplot(kurtosis.spikein[2,],"spikein_kurtosis"),  
             spikein_mdplot(kurtosis.spikein[3,],"spikein_kurtosis"),  spikein_mdplot(kurtosis.spikein[4,],"spikein_kurtosis"),  spikein_mdplot(kurtosis.spikein[5,],"spikein_kurtosis"),  
             spikein_mdplot(kurtosis.spikein[6,],"spikein_kurtosis"),  spikein_mdplot(kurtosis.spikein[7,],"spikein_kurtosis"),  spikein_mdplot(kurtosis.spikein[8,],"spikein_kurtosis"), 
             spikein_mdplot(kurtosis.spikein[9,],"spikein_kurtosis"),  spikein_mdplot(kurtosis.spikein[10,],"spikein_kurtosis"), spikein_mdplot(kurtosis.spikein[11,],"spikein_kurtosis"), 
             spikein_mdplot(kurtosis.spikein[12,],"spikein_kurtosis"), spikein_mdplot(kurtosis.spikein[13,],"spikein_kurtosis"), spikein_mdplot(kurtosis.spikein[14,],"spikein_kurtosis"), 
             spikein_mdplot(kurtosis.spikein[15,],"spikein_kurtosis"), spikein_mdplot(kurtosis.spikein[16,],"spikein_kurtosis"), spikein_mdplot(kurtosis.spikein[17,],"spikein_kurtosis"),
             ncol = 3)
grid.arrange(gmode,                                            spikein_mdplot(mode.spikein[1,],"spikein_mode"),  spikein_mdplot(mode.spikein[2,],"spikein_mode"),  
             spikein_mdplot(mode.spikein[3,],"spikein_mode"),  spikein_mdplot(mode.spikein[4,],"spikein_mode"),  spikein_mdplot(mode.spikein[5,],"spikein_mode"),  
             spikein_mdplot(mode.spikein[6,],"spikein_mode"),  spikein_mdplot(mode.spikein[7,],"spikein_mode"),  spikein_mdplot(mode.spikein[8,],"spikein_mode"), 
             spikein_mdplot(mode.spikein[9,],"spikein_mode"),  spikein_mdplot(mode.spikein[10,],"spikein_mode"), spikein_mdplot(mode.spikein[11,],"spikein_mode"), 
             spikein_mdplot(mode.spikein[12,],"spikein_mode"), spikein_mdplot(mode.spikein[13,],"spikein_mode"), spikein_mdplot(mode.spikein[14,],"spikein_mode"), 
             spikein_mdplot(mode.spikein[15,],"spikein_mode"), spikein_mdplot(mode.spikein[16,],"spikein_mode"), spikein_mdplot(mode.spikein[17,],"spikein_mode"),
             ncol = 3)
dev.off()

