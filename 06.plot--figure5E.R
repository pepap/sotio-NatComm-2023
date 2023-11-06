library(data.table)
library(maftools)

ODIR="01.allSamples/"

#>> @pepap : unify annotations & clinical data
cli.tmp.dt                <- fread( "plot--figure5E-data/col_annot_allSamples.txt",header=T )
xSAMPLES                  <- fread( "plot--figure5E-data/samples.csv",header=T )
xSAMPLES                  <- xSAMPLES[ xSAMPLES[["Subject Identifier"]] %in% cli.tmp.dt[["Subject Identifier"]] ]
cli.dt                    <- fread( "plot--figure5E-data/oncoplot_annotation.txt",header=T )
cli.dt                    <- cli.dt[ ,c( "Subject Identifier","Treatment arm","OS_outcome","PFS_outcome" ),with=F ]
colnames(cli.dt)          <- c("Tumor_Sample_Barcode","Study_arm","OS","PFS")
cli.dt                    <- merge( cli.dt,cli.tmp.dt,by.x="Tumor_Sample_Barcode",by.y="Subject Identifier",all.y=T,sort=F )

outMAFsom="SomaticStatus.somatic"

#>> @pepap : unify treatments
treatmentArm      <- c("limegreen","cyan","red"); names(treatmentArm) <- c("Arm A","Arm B","Arm C")
OS                <- c("#FFCC00","#333333");      names(OS)           <- c("Alive","Dead")
PFS               <- c("#999999","#FF6600");      names(PFS)          <- c("Relapse","Relapse-free")
TMB_med           <- c("#FFCC00","#333333");      names(TMB_med)      <- c("High","Low")
TLS_med           <- c("#FFCC00","#333333");      names(TLS_med)      <- c("High","Low")
annotation.colors <- list( "Study_arm"=treatmentArm,TMB_med=TMB_med,TLS_med=TLS_med )

xcol <- c("forestgreen","forestgreen"); names(xcol) <- c("Missense_Mutation","Multi_Hit")

dir.create( path=ODIR,showWarnings=F )
for ( iii in seq_along(outMAFsom) ) {

cat( " => ",outMAFsom[iii],"\n",sep="" )

in.maf <- fread(paste0("plot--figure5E-data/",outMAFsom[iii],".maf"),header=T)[ Tumor_Sample_Barcode %in% xSAMPLES[["Subject Identifier"]] ]
write.table( x=in.maf,file=paste0(ODIR,outMAFsom[iii],".filtered.maf"),sep="\t",row.names=F,quote=F )

in.maf <- read.maf(    maf=paste0(ODIR,outMAFsom[iii],".filtered.maf"), clinicalData=cli.dt )
cat(" => ",ODIR,"\n",sep="" )

pdf(file=paste0(ODIR,outMAFsom[iii],"-top20-genes.pdf"),height=076,width=50,pointsize=40)
oncoplot(
 maf=in.maf,showTumorSampleBarcodes=F,barcode_mar=13,gene_mar=8,
 colors=xcol,
 annotationColor=annotation.colors,sortByAnnotation=T,groupAnnotationBySize=F,removeNonMutated=F,top=20,
 clinicalFeatures=c("Study_arm","TMB_med","TLS_med"),
 annotationOrder=c("Arm A","Arm B","Arm C"),
 fontSize=1.1,SampleNamefontSize=1.2,legendFontSize=2,annotationFontSize=1.7,bgCol="#CCCCCC",titleFontSize=1.8,
 borderCol="white",annoBorderCol="white"
)
dev.off()

}

