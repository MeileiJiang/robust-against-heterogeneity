###############################################################################
## myelement.R
##
## my element for plots
## 
## Author: Haaland
###############################################################################
require(grid)

my_theme = function (base_size = 12, base_family = "") 
{
	theme_grey(base_size = base_size, base_family = base_family) %+replace% 
			theme(axis.text = element_text(size = rel(0.8)), 
					axis.ticks = element_line(colour = "black"), 
					legend.key = element_rect(colour = "grey80"), 
					panel.background = element_rect(fill = "white", colour = NA), 
					panel.border = element_rect(fill = NA, colour = "grey10"), 
					panel.grid.major = element_line(colour = "grey70", size = 0.2), 
					panel.grid.minor = element_line(colour = "grey90", size = 0.5), 
					strip.background = element_rect(fill = "grey80", colour = "grey50"), 
					strip.background = element_rect(fill = "grey80", colour = "grey50"))
}

## color blind package from R graphics cookbook
##The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
## the standard ggplot2 palette but with 5 categories fixed and in the order that I like
fixedPalette = hue_pal(h=c(0,360)+15,c=100,l=65,h.start=0,direction=1)(5)[c(1,5,2,3,4)]

## blank background
theme_bare <- theme(
		axis.line = element_blank(), 
		axis.text.x = element_blank(), 
		axis.text.y = element_blank(),
		axis.ticks = element_blank(), 
		axis.title.x = element_blank(), 
		axis.title.y = element_blank(), 
		#axis.ticks.length = unit(0, "lines"), # Error 
		axis.ticks.margin = unit(c(0,0,0,0), "lines"), 
		legend.position = "none", 
		panel.background = element_rect(fill = "white",color=NA), 
		panel.border = element_blank(), 
		panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.margin = unit(c(0,0,0,0), "lines"), 
		plot.background = element_rect(fill = "white"),
		plot.margin = unit(c(0,0,0,0), "lines")
)