#######################################################################################################
## Buettner_Log10.R
## Author: Meilei
#######################################################################################################
library(ggplot2)
library(dplyr)

load("~/researchspace/robust-against-heterogeneity/Mengjie/Single_cell_data/Buettner_Log10.RData")

dim(part1)
# [1]  182 8942
dim(part2)
# [1] 182 629
dim(Spikein)
# [1] 182  92

glimpse(part1)
