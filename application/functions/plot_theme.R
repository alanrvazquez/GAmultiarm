################################################################################
#
# THEME FOR PLOTS IN GGPLOT2
# 
#
# Author: Alan Vazquez
# Affiliation: Tecnologico de Monterrey
# Email: alanrvazquez@tec.mx
#
################################################################################

library(ggplot2)

ytext <- 25
xtext <- 20 
legtext <- 20
psize <- 8
xtitle <- 20
ytitle <- 20
legtitle <- 20
plot.theme <- theme(axis.text.y  = element_text(size=ytext, colour = 'black'), axis.text.x  = element_text(size=xtext, colour = 'black'),
                    axis.title.y  = element_text(size=ytitle, vjust=0.35), axis.title.x  = element_text(size=xtitle,vjust=0),
                    panel.background = element_rect(fill = "white", colour = 'black'),
                    legend.text = element_text(size = legtext), legend.title = element_text(size=legtitle),
                    plot.title = element_text(lineheight=.8, face="bold", size = 20,hjust = 0.5),
                    strip.text = element_text(size=20), strip.text.x = element_text(size=20, face="bold"), 
                    legend.position="none", 
                    panel.grid.major = element_line(linewidth = 0.5, linetype = 'dashed',
                                                    colour = "grey"), 
                    panel.grid.minor = element_line(linewidth = 0.25, linetype = 'dashed',
                                                    colour = "grey"))
