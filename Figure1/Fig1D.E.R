library(Seurat)
library(dplyr)
library(MCPcounter)
library(ggplot2)

library(data.table)
MH_GENES_DT               <- fread( "Angelova_et_al.csv",header=T )
MH_GENES_DT[["CellType"]] <- make.names( MH_GENES_DT[["CellType"]] )

mh.genes <- unique( MH_GENES_DT[["Name"]] )
sc.genes <- unique( rownames( sample2_norm[["SCT"]]@data ) )

for ( ign in mh.genes[ !( mh.genes %in% sc.genes ) & !( mh.genes %in% c("HCP5","KIR2DL1","KIR2DL3","HLA-DR","IL12","COX2","TFGFB5") ) ] ) {

 tmp.dt <- gn2synonyms( GN=ign,DESC=T )
 if( nrow(tmp.dt)==1 ) {
  MH_GENES_DT[ Name==ign ][["Name"]] <- tmp.dt[["Symbol"]]
 } else                {
  cat( "\n" )
  print( ign )
  print( tmp.dt )
  cat( "\n" )
 }

}

print( dim(MH_GENES_DT) )
MH_GENES_DT <- unique( MH_GENES_DT )
print( dim(MH_GENES_DT) )
MH_GENES_DT <- MH_GENES_DT[ Name %in% sc.genes ]
print( dim(MH_GENES_DT) )

if ( any( table(MH_GENES_DT[["CellType"]])==1 ) ) { stop("\n !!! some \"CellTypes\" from \"MH_GENE_LIST\" contain only one single gene !!!\n") }



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
    scale_fill_viridis_c( option="C" )          + 
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

