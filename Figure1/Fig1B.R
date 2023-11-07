library(ComplexHeatmap)
library(circlize)
my_data=read.table(file="Fig1B.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)
mat = as.matrix(my_data)

mycols <- colorRamp2(breaks = c(0,1,2,3),colors = c("green","red","blue","gold"))

set.seed(123456)

circPlot = function() {

circos.par( gap.after=c(25) )
par(cex=2.5)
circos.heatmap( mat,split=NULL,col=mycols,track.height=0.4,cell_width =0.2,rownames.side="none",cluster=T )
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    if( CELL_META$sector.numeric.index==1 ) { # the last sector
        cn = colnames(mat)
        n  = length(cn)
        circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
            1:n - 0.5, cn, 
            cex = 0.5, adj = c(0, 0.5), facing = "inside")


        circos.segments( c(0,seq(nrow(mat))),0,                  c(0,seq(nrow(mat))),        ncol(mat)   )
        circos.segments( 0,                  c(0,seq(ncol(mat))),        nrow(mat),  c(0,seq(ncol(mat))) )

    }
}, bg.border = "black")

circos.clear()

}

lgd_exp = Legend( title="Expression",col_fun=mycols )

library(gridBase)
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region

pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
    just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circPlot()
upViewport()

h = dev.size()[2]
lgd_list = packLegend( lgd_exp, max_height = unit(0.9*h, "inch") )
draw(lgd_list, x = circle_size, just = "left")

