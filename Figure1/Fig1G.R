library(forestplot)
library(dplyr)
library(data.table)
data=fread(file="Fig1G.txt",sep="\t")
labeltext = c("study", "p.value")
hrzl_lines = list("1" = gpar(lty = 2),"2" = gpar(lty = 2),"3" = gpar(lty = 2),"4" = gpar(lty = 2))
p=data %>%
  group_by(group) %>%
  forestplot(labeltext = c(study),clip = c(-.1, 3),
             shapes_gp = fpShapesGp(box = c("green","red","black") %>% lapply(function(x) gpar(fill = x, col = "#555555",cex=1.2)),
                                    default = gpar(vertices = TRUE)),
             ci.vertices = TRUE,
             ci.vertices.height = 0.05,txt_gp = fpTxtGp(ticks=gpar(cex=1),xlab=gpar(cex=1),label=gpar(cex=1),legend=gpar(cex=1)),
             boxsize = .15,zero = 1,line.margin = .1,hrzl_lines =hrzl_lines,
             col = fpColors(box = "black",line="black"),xlab  = "Hazard ratio")
plot(p)