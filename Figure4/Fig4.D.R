library(Seurat)
library(dplyr)
library(MCPcounter)
library(ggplot2)
library(scales)
library(data.table)

MH_GENES_DT               <- fread( "ref.dataset.csv",header=T )
MH_GENES_DT[["CellType"]] <- make.names( MH_GENES_DT[["CellType"]] )

mh.genes <- unique( MH_GENES_DT[["Name"]] )
sc.genes <- unique( rownames( sample2_norm[["SCT"]]@data ) )

print( dim(MH_GENES_DT) )
MH_GENES_DT <- unique( MH_GENES_DT )
print( dim(MH_GENES_DT) )
MH_GENES_DT <- MH_GENES_DT[ Name %in% sc.genes ]
print( dim(MH_GENES_DT) )


SGS_types  <- unique( MH_GENES_DT[["CellType"]] )

for( CT in SGS_types ){

 SGS_scores   <- sample2_norm[["SCT"]]@data[ MH_GENES_DT[ CellType==CT , Name ] , ]
 if ( nrow(SGS_scores)!=length(MH_GENES_DT[ CellType==CT , Name ]) ) { stop("\n !!! not all genes from \"MH_GENE_LIST\" found in expression matrix \"sample2_norm\" !!!\n") }

 sample2_norm <- AddMetaData( object=sample2_norm,metadata=colMeans( SGS_scores ),col.name=CT )

}

hide_axis <- theme(axis.title.x=element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank(),
			panel.grid.minor = element_blank(),
panel.grid.major = element_blank())

plot_list_mcp <-
 lapply(
  X   = SGS_types,
  FUN = function(CT) {

   plot_list_mcp <-
    SpatialPlot(
     sample2_norm,
     features=CT,
     image.alpha=0,
     pt.size.factor=1.8
    )                                           +
    scale_fill_gradientn(colors = viridis_pal(option="C")(9), limits=c(0,1),na.value = "yellow")          + 
    DarkTheme()						+
hide_axis                                 +
    theme( text=element_text( size=14 ) )       + 
    theme( text=element_text( face="bold" ) )   +
    theme( legend.text=element_text( size=7 ) )

                    }
       )

library(gridExtra)

ggsave(
 filename="SGS_types.oneHeatmapPerPage.pdf",
 plot=marrangeGrob( plot_list_mcp,nrow=1,                                                   ncol=1 ),
 width=15,height=15
)
ggsave(
 filename="SGS_types.onePageAll.pdf",
 plot=marrangeGrob( plot_list_mcp,nrow=ceiling(sqrt( length( unique(MH_GENES_DT[["CellType"]]) ) )),ncol=ceiling(sqrt( length( unique(MH_GENES_DT[["CellType"]]) ) )) ),
 width=49,height=49
)

