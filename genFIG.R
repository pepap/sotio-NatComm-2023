library(data.table)
library(scales)

source("~/Scripts/R_func/cut2breaks.R")

load( "/storage/brno1-cerit/home/pepap/SOTIO/MH-mutAnal/BAM/00.SangerQual/GRCm38_68/STRELKA/mergedVars.annot.dt.rda",verbose=T )
ILL.mergedVars.annot.dt <- mergedVars.annot.dt
ILL.mergedVars.annot.dt <- ILL.mergedVars.annot.dt[ !is.na(cdna_start) ]
ILL.mergedVars.annot.dt <- ILL.mergedVars.annot.dt[ order( SAMPLE,VAR,ensembl_gene_stable_id,trx_len,decreasing=T ) ][ !duplicated(paste0(VAR,".",SAMPLE)) ]
rm( list=c("mergedVars.annot.dt") )

load( "/storage/brno1-cerit/home/pepap/SOTIO/MH--mutAnal-20230725/BAM/STRELKA/mergedVars.annot.dt.rda",              verbose=T )
AGI.mergedVars.annot.dt <- mergedVars.annot.dt
AGI.mergedVars.annot.dt <- AGI.mergedVars.annot.dt[ !is.na(cdna_start) ]
AGI.mergedVars.annot.dt <- AGI.mergedVars.annot.dt[ order( SAMPLE,VAR,ensembl_gene_stable_id,trx_len,decreasing=T ) ][ !duplicated(paste0(VAR,".",SAMPLE)) ]
rm( list=c("mergedVars.annot.dt") )

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

#!#COLS=c("cyan","red","limegreen")
COLS <- hue_pal()(9)

if (F) {
CLEN=30
mut.tab <- rbind( ILL.mergedVars.annot.dt,AGI.mergedVars.annot.dt )[ , table(consequence_allele_string) ]
br1.mat <- matrix( data=0,ncol=CLEN-1,nrow=length(mut.tab),dimnames=list( names(mut.tab),seq(CLEN-1) ) )
br2.mat <- matrix( data=0,ncol=CLEN-1,nrow=length(mut.tab),dimnames=list( names(mut.tab),seq(CLEN-1) ) )
br3.mat <- matrix( data=0,ncol=CLEN-1,nrow=length(mut.tab),dimnames=list( names(mut.tab),seq(CLEN-1) ) )
id1.mat <- matrix( data=0,ncol=CLEN-1,nrow=length(mut.tab),dimnames=list( names(mut.tab),seq(CLEN-1) ) )
id2.mat <- matrix( data=0,ncol=CLEN-1,nrow=length(mut.tab),dimnames=list( names(mut.tab),seq(CLEN-1) ) )
id3.mat <- matrix( data=0,ncol=CLEN-1,nrow=length(mut.tab),dimnames=list( names(mut.tab),seq(CLEN-1) ) )
sa1.mat <- matrix( data=0,ncol=CLEN-1,nrow=length(mut.tab),dimnames=list( names(mut.tab),seq(CLEN-1) ) )
sa2.mat <- matrix( data=0,ncol=CLEN-1,nrow=length(mut.tab),dimnames=list( names(mut.tab),seq(CLEN-1) ) )
sa3.mat <- matrix( data=0,ncol=CLEN-1,nrow=length(mut.tab),dimnames=list( names(mut.tab),seq(CLEN-1) ) )
b1MAX <- 0
b2MAX <- 0
b3MAX <- 0
i1MAX <- 0
i2MAX <- 0
i3MAX <- 0
s1MAX <- 0
s2MAX <- 0
s3MAX <- 0

II=1
LB=(-1)
for ( UB in seq( 0.0,1.0,length.out=CLEN )[-1] ) {

cat( " < ",LB," ; ",UB," >\n",sep="" )

br1.tmp                      <- br1.ann[ ( (cdna_start/trx_len)>LB ) & ( (cdna_start/trx_len)<=UB ) , table(consequence_allele_string) ]
br1.mat[ names(br1.tmp),II ] <- as.vector(br1.tmp)
br2.tmp                      <- br2.ann[ ( (cdna_start/trx_len)>LB ) & ( (cdna_start/trx_len)<=UB ) , table(consequence_allele_string) ]
br2.mat[ names(br2.tmp),II ] <- as.vector(br2.tmp)
br3.tmp                      <- br3.ann[ ( (cdna_start/trx_len)>LB ) & ( (cdna_start/trx_len)<=UB ) , table(consequence_allele_string) ]
br3.mat[ names(br3.tmp),II ] <- as.vector(br3.tmp)

id1.tmp                      <- id1.ann[ ( (cdna_start/trx_len)>LB ) & ( (cdna_start/trx_len)<=UB ) , table(consequence_allele_string) ]
id1.mat[ names(id1.tmp),II ] <- as.vector(id1.tmp)
id2.tmp                      <- id2.ann[ ( (cdna_start/trx_len)>LB ) & ( (cdna_start/trx_len)<=UB ) , table(consequence_allele_string) ]
id2.mat[ names(id2.tmp),II ] <- as.vector(id2.tmp)
id3.tmp                      <- id3.ann[ ( (cdna_start/trx_len)>LB ) & ( (cdna_start/trx_len)<=UB ) , table(consequence_allele_string) ]
id3.mat[ names(id3.tmp),II ] <- as.vector(id3.tmp)

sa1.tmp                      <- sa1.ann[ ( (cdna_start/trx_len)>LB ) & ( (cdna_start/trx_len)<=UB ) , table(consequence_allele_string) ]
sa1.mat[ names(sa1.tmp),II ] <- as.vector(sa1.tmp)
#sa1.mat[ names(sa1.tmp[ names(mut.tab) ])[ !is.na(names(sa1.tmp[ names(mut.tab) ])) ] , II ] <- as.vector(sa1.tmp[ names(mut.tab) ])[ !is.na(names(sa1.tmp[ names(mut.tab) ])) ]
sa2.tmp                      <- sa2.ann[ ( (cdna_start/trx_len)>LB ) & ( (cdna_start/trx_len)<=UB ) , table(consequence_allele_string) ]
sa2.mat[ names(sa2.tmp),II ] <- as.vector(sa2.tmp)
#sa2.mat[ names(sa2.tmp[ names(mut.tab) ])[ !is.na(names(sa2.tmp[ names(mut.tab) ])) ] , II ] <- as.vector(sa2.tmp[ names(mut.tab) ])[ !is.na(names(sa2.tmp[ names(mut.tab) ])) ]
sa3.tmp                      <- sa3.ann[ ( (cdna_start/trx_len)>LB ) & ( (cdna_start/trx_len)<=UB ) , table(consequence_allele_string) ]
sa3.mat[ names(sa3.tmp),II ] <- as.vector(sa3.tmp)
#sa3.mat[ names(sa3.tmp[ names(mut.tab) ])[ !is.na(names(sa3.tmp[ names(mut.tab) ])) ] , II ] <- as.vector(sa3.tmp[ names(mut.tab) ])[ !is.na(names(sa3.tmp[ names(mut.tab) ])) ]

b1MAX <- max( c(b1MAX,max(sum( as.vector(br1.tmp),na.rm=T ))),na.rm=T )
b2MAX <- max( c(b2MAX,max(sum( as.vector(br2.tmp),na.rm=T ))),na.rm=T )
b3MAX <- max( c(b3MAX,max(sum( as.vector(br3.tmp),na.rm=T ))),na.rm=T )
i1MAX <- max( c(i1MAX,max(sum( as.vector(id1.tmp),na.rm=T ))),na.rm=T )
i2MAX <- max( c(i2MAX,max(sum( as.vector(id2.tmp),na.rm=T ))),na.rm=T )
i3MAX <- max( c(i3MAX,max(sum( as.vector(id3.tmp),na.rm=T ))),na.rm=T )
s1MAX <- max( c(s1MAX,max(sum( as.vector(sa1.tmp),na.rm=T ))),na.rm=T )
s2MAX <- max( c(s2MAX,max(sum( as.vector(sa2.tmp),na.rm=T ))),na.rm=T )
s3MAX <- max( c(s3MAX,max(sum( as.vector(sa3.tmp),na.rm=T ))),na.rm=T )

II <- (II+1)
LB <- UB

}
allMAX <- max( c(b1MAX,b2MAX,b3MAX,i1MAX,i2MAX,i3MAX,s1MAX,s2MAX,s3MAX) )
}

if (F) {
COLS <- hue_pal()(9)
#!#for ( i.brk in c(0.01,0.02,0.05) ) {
for ( i.brk in c(0.01) ) {

cat( " >> points : ",1/i.brk,"\n",sep="" )
BRKS=seq(0,1,i.brk)

inp.tab <-
 rbind(
  cut2breaks( xVALS=br1.ann[ !is.na(cdna_start) , cdna_start/trx_len ],xBREAKS=BRKS ),
  cut2breaks( xVALS=br2.ann[ !is.na(cdna_start) , cdna_start/trx_len ],xBREAKS=BRKS ),
  cut2breaks( xVALS=br3.ann[ !is.na(cdna_start) , cdna_start/trx_len ],xBREAKS=BRKS ),
  cut2breaks( xVALS=id1.ann[ !is.na(cdna_start) , cdna_start/trx_len ],xBREAKS=BRKS ),
  cut2breaks( xVALS=id2.ann[ !is.na(cdna_start) , cdna_start/trx_len ],xBREAKS=BRKS ),
  cut2breaks( xVALS=id3.ann[ !is.na(cdna_start) , cdna_start/trx_len ],xBREAKS=BRKS ),
  cut2breaks( xVALS=sa1.ann[ !is.na(cdna_start) , cdna_start/trx_len ],xBREAKS=BRKS ),
  cut2breaks( xVALS=sa2.ann[ !is.na(cdna_start) , cdna_start/trx_len ],xBREAKS=BRKS ),
  cut2breaks( xVALS=sa3.ann[ !is.na(cdna_start) , cdna_start/trx_len ],xBREAKS=BRKS )
 )
rownames(inp.tab) <- c( "BR5.1","BR5.2","BR5.3","ID8.1","ID8.2","ID8.3","SAMPLE1","SAMPLE2","SAMPLE3" )

pdf( file=paste0("MH--mutAnal-all-VarsPerGene.art.points",sprintf("%03d",1/i.brk),".",pepDate(),".pdf"),pointsize=20,width=50,height=10 )
par( mfrow=c(1,1),bg="white" )
x.tmp <- colMedians( barplot( inp.tab,beside=T,plot=F ) )
barplot(
 inp.tab,col=COLS,legend.text=rownames(inp.tab),args.legend=list("topright",fill=COLS,bty="n"),beside=T,font.lab=2,
 main="Mutation variants : gene coverage",xaxt="n",xlab="",ylab="Number of SNVs",plot=T
)
axis( side=1,at=x.tmp[ c(1,ceiling(length(x.tmp)/2),length(x.tmp)) ],labels=c("5'","Relative position in gene","3'"),font=2 )
dev.off()

}
}

if (F) {
pdf( file=paste0("MH--mutAnal-all-VarsPerGene.art.BoxPlots.",pepDate(),".pdf"),pointsize=20,width=17,height=10 )
 par( mfrow=c(1,1),bg="white" )
 boxplot( as.data.table(t(inp.tab)),col=COLS,main="Mutation variants : gene coverage",ylab="Number of SNVs",font.lab=2 )
dev.off()
}

if (F) {
pdf( file=paste0("MH--mutAnal-all-VarsPerGene.art.BarPlots.",pepDate(),".pdf"),pointsize=20,width=17,height=10 )
 par( mfrow=c(1,1),bg="white" )
 barplot( rowSums(inp.tab),         col=COLS,main="Mutation variants : gene coverage",ylab="Number of SNVs",font.lab=2 )
dev.off()
}

#!#COLS <- c( "red","limegreen","blue" )
#!#COLS <- c( "red","blue" )
COLS <- c( "black","limegreen" )
for ( i.brk in c(0.01) ) {

cat( " >> points : ",1/i.brk,"\n",sep="" )
BRKS=seq(0,1,i.brk)

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
#!#inp.tab                       <- inp.tab[ c("ID8","SAMPLE3","SAMPLE2") , ]
mutPerRelGenePos.ID8.S2.dt    <- inp.tab
#!#mutPerRelGenePos.ID8.S3.S2.dt <- inp.tab
#!#save( mutPerRelGenePos.ID8.S3.S2.dt,file="mutPerRelGenePos.ID8.S3.S2.dt.rda" )
save( mutPerRelGenePos.ID8.S2.dt,file="mutPerRelGenePos.ID8.S2.dt.rda" )

#!#pdf( file=paste0("MH--mutAnal-ID8.S3.S2-VarsPerGene.art.points",sprintf("%03d",1/i.brk),".",pepDate(),".pdf"),pointsize=20,width=20,height=10 )
#!#pdf( file=paste0("MH--mutAnal-ID8.S2-VarsPerGene.art.points",sprintf("%03d",1/i.brk),".",pepDate(),".pdf"),pointsize=20,width=20,height=10 )
#!#pdf( file=paste0("MH--mutAnal-ID8.S2-VarsPerGene.art.points",sprintf("%03d",1/i.brk),"-EMPTY.",pepDate(),".pdf"),pointsize=20,width=20,height=10 )
pdf( file=paste0("MH--mutAnal-ID8.S2-VarsPerGene.art.points",sprintf("%03d",1/i.brk),"-EMPTY-DIFFCOL.",pepDate(),".pdf"),pointsize=20,width=20,height=10 )
par( mfrow=c(1,1),bg="white" )
x.tmp <- colMedians( barplot( inp.tab,beside=T,width=5,plot=F ) )
barplot(
 inp.tab,width=5,col=COLS,legend.text=NULL,args.legend=NULL,beside=T,font.lab=2,main="",xaxt="n",xlab="",ylab="",plot=T
#!# inp.tab,width=5,col=COLS,legend.text=c("ID8",expression( "SO1 "*italic("Brca"^"-/-")*" "*italic("Trp53"^"-/-")*', Hras, Myc')),args.legend=list("topright",fill=COLS,bty="n"),beside=T,font.lab=2,
#!# main="Mutation variants : gene coverage",xaxt="n",xlab="",ylab="Number of SNVs",plot=T
)
#!#axis( side=1,at=x.tmp[ c(1,ceiling(length(x.tmp)/2),length(x.tmp)) ],labels=c("5'","Relative position in gene","3'"),font=2 )
dev.off()

if (F) {
#!#pdf( file=paste0("MH--mutAnal-ID8.S3.S2-VarsPerGene.art.points",sprintf("%03d",1/i.brk),".",pepDate(),".2.pdf"),pointsize=20,width=15,height=10 )
pdf( file=paste0("MH--mutAnal-ID8.S2-VarsPerGene.art.points",sprintf("%03d",1/i.brk),".",pepDate(),".2.pdf"),pointsize=20,width=15,height=10 )
par( mfrow=c(1,1),bg="white" )
x.tmp <- colMedians( barplot( inp.tab,beside=T,width=5,plot=F ) )
barplot(
 inp.tab,width=5,col=rev(COLS),legend.text=rownames(inp.tab),args.legend=list("topright",fill=COLS,bty="n"),beside=T,font.lab=2,
 main="Mutation variants : gene coverage",xaxt="n",xlab="",ylab="Number of SNVs",plot=T
)
axis( side=1,at=x.tmp[ c(1,ceiling(length(x.tmp)/2),length(x.tmp)) ],labels=c("5'","Relative position in gene","3'"),font=2 )
dev.off()
}

}

