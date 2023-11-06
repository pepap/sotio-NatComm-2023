library(data.table)
library(scales)
library(matrixStats)

#>> @pepap : create "cut2breaks" function
cut2breaks <- function(xVALS,xBREAKS=c(0,0.1,1,10,100,1000)) {

 xBREAKS <- sort(xBREAKS)
 xRANGE  <- range( x=xVALS,na.rm=T )

 if ( xBREAKS[1]               > min(xRANGE) ) { xBREAKS <- c( -Inf,xBREAKS     ) }
 if ( xBREAKS[length(xBREAKS)] < max(xRANGE) ) { xBREAKS <- c(      xBREAKS,Inf ) }

 OUT <- cut( x=xVALS,breaks=xBREAKS )

 return( table(OUT) )

}

#>> @pepap : load annotated variants from ILLUMINA
load( "ILL.mergedVars.annot.dt.rda",verbose=T )
ILL.mergedVars.annot.dt <- mergedVars.annot.dt
ILL.mergedVars.annot.dt <- ILL.mergedVars.annot.dt[ !is.na(cdna_start) ]
ILL.mergedVars.annot.dt <- ILL.mergedVars.annot.dt[ order( SAMPLE,VAR,ensembl_gene_stable_id,trx_len,decreasing=T ) ][ !duplicated(paste0(VAR,".",SAMPLE)) ]
rm( list=c("mergedVars.annot.dt") )

#>> @pepap : load annotated variants from AGILENT
load( "AGI.mergedVars.annot.dt.rda",verbose=T )
AGI.mergedVars.annot.dt <- mergedVars.annot.dt
AGI.mergedVars.annot.dt <- AGI.mergedVars.annot.dt[ !is.na(cdna_start) ]
AGI.mergedVars.annot.dt <- AGI.mergedVars.annot.dt[ order( SAMPLE,VAR,ensembl_gene_stable_id,trx_len,decreasing=T ) ][ !duplicated(paste0(VAR,".",SAMPLE)) ]
rm( list=c("mergedVars.annot.dt") )

#>> @pepap : formatting
ILL.mergedVars.annot.dt[["Chrom"]] <- ILL.mergedVars.annot.dt[["CHR"]]
ILL.mergedVars.annot.dt[["Start"]] <- ILL.mergedVars.annot.dt[["CRD"]]
ILL.mergedVars.annot.dt[["End"]]   <- ILL.mergedVars.annot.dt[["CRD"]]
AGI.mergedVars.annot.dt[["Chrom"]] <- AGI.mergedVars.annot.dt[["CHR"]]
AGI.mergedVars.annot.dt[["Start"]] <- AGI.mergedVars.annot.dt[["CRD"]]
AGI.mergedVars.annot.dt[["End"]]   <- AGI.mergedVars.annot.dt[["CRD"]]
br1.ann                            <- ILL.mergedVars.annot.dt[ SAMPLE=="BR5.1"   ]
br2.ann                            <- ILL.mergedVars.annot.dt[ SAMPLE=="BR5.2"   ]
br3.ann                            <- ILL.mergedVars.annot.dt[ SAMPLE=="BR5.3"   ]
id1.ann                            <- ILL.mergedVars.annot.dt[ SAMPLE=="ID8.1"   ]
id2.ann                            <- ILL.mergedVars.annot.dt[ SAMPLE=="ID8.2"   ]
id3.ann                            <- ILL.mergedVars.annot.dt[ SAMPLE=="ID8.3"   ]
sa1.ann                            <- AGI.mergedVars.annot.dt[ SAMPLE=="SAMPLE1" ]
sa2.ann                            <- AGI.mergedVars.annot.dt[ SAMPLE=="SAMPLE2" ]
sa3.ann                            <- AGI.mergedVars.annot.dt[ SAMPLE=="SAMPLE3" ]

#>> @pepap : define colours
COLS <- c( "black","limegreen" )
#>> @pepap : set number of points per gene to 100
i.brk <- c(0.01)

cat( " >> points : ",1/i.brk,"\n",sep="" )
BRKS=seq(0,1,i.brk)

#>> @pepap : use median of the three replicates as the ID8 sample
inp.tab <-
 rbind(
  cut2breaks( xVALS=id1.ann[ !is.na(cdna_start) , cdna_start/trx_len ],xBREAKS=BRKS ),
  cut2breaks( xVALS=id2.ann[ !is.na(cdna_start) , cdna_start/trx_len ],xBREAKS=BRKS ),
  cut2breaks( xVALS=id3.ann[ !is.na(cdna_start) , cdna_start/trx_len ],xBREAKS=BRKS )
 )
inp.tab <-
 rbind(
  apply( X=inp.tab,MARGIN=2,FUN=median ),
  cut2breaks( xVALS=sa2.ann[ !is.na(cdna_start) , cdna_start/trx_len ],xBREAKS=BRKS ),
  cut2breaks( xVALS=sa3.ann[ !is.na(cdna_start) , cdna_start/trx_len ],xBREAKS=BRKS )
 )
rownames(inp.tab)             <- c( "ID8","SAMPLE2","SAMPLE3" )
inp.tab                       <- inp.tab[ c("ID8","SAMPLE2") , ]

pdf( file=paste0("MH--mutAnal-ID8.S2-VarsPerGene.art.points",sprintf("%03d",1/i.brk),"-EMPTY-DIFFCOL.",pepDate(),".pdf"),pointsize=20,width=20,height=10 )
par( mfrow=c(1,1),bg="white" )
x.tmp <- colMedians( barplot( inp.tab,beside=T,width=5,plot=F ) )
barplot(
 inp.tab,width=5,col=COLS,legend.text=NULL,args.legend=NULL,beside=T,font.lab=2,main="",xaxt="n",xlab="",ylab="",plot=T
)
dev.off()
