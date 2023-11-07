library(ggplot2)
library(data.table)
library(tidyverse)

#create separate variable for "score" (point sizes) & "Scaled Average Expression" (color gradient)
a=fread(file="Fig4.F.txt",header=T,sep="\t")

xdp <-
 ggplot( a,aes(reorder(Group,ID),reorder(Gene,-ID),color=score) )    +
 geom_point( aes(size=score) )                                   +
 scale_colour_gradient( low="green",high="red")                     +
 guides( color=guide_colorbar( title='Scaled Average Expression' ),size=guide_legend("score") ) +
 theme( axis.text.x=element_text(color="black", angle=90 ) )                       +
 theme_classic()+theme(axis.text.x= element_text(color="black"))
xdp+theme(axis.text.y= element_text(color="black"))
                                                        
ggsave("Fig4.F.pdf", width = 18, height = 30, units = "cm",dpi = 600)