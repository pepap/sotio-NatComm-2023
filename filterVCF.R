library(data.table)
library(VariantAnnotation)
library(GenomicRanges)

cat("\n @ pepap-functions loaded : \"filterVCF\", \"purgeVCF\"\n\n",sep="")

#!#load( "/storage/brno1-cerit/home/pepap/SOTIO/ET-mutAnal/S21_032/Exome_regions/Twist_Mouse_Exome_Target_Rev1_7APR20.gr.rda",verbose=T )
load( "/storage/brno3-cerit/home/pepap/Exome_regions/agi.mm10.gr.rda",   verbose=T )
load( "/storage/brno3-cerit/home/pepap/Exome_regions/agi.ensembl.gr.rda",verbose=T )
load( "/storage/brno3-cerit/home/pepap/Exome_regions/Twist_Mouse_Exome_Target_Rev1_7APR20.gr.rda",verbose=T )

#= (1)
filterVCF <- function( VCFfile,FILTER.PASS=T,sele.gr=Twist_Mouse_Exome_Target_Rev1_7APR20.gr ) {

 cat( " > Reading VCF-file :\n    ",VCFfile,"\n",sep="" )
 inp.vcf <- readVcf( file=VCFfile )

 if ( FILTER.PASS ) {
  cat( " > Filter \"PASS\" variants only\n",sep="" )
  inp.vcf <- inp.vcf[ mcols(rowRanges(inp.vcf))[["FILTER"]]=="PASS" ]
 }

 cat( " > Variants from \"sele.gr\" regions only\n",sep="" )
 exo.ovr <- findOverlaps( query=rowRanges(inp.vcf),subject=sele.gr,ignore.strand=T )
 inp.vcf <- inp.vcf[ unique(queryHits(exo.ovr)) ]

 cat( " > Done\n",sep="" )
 return( inp.vcf )

}

#= (2)
purgeVCF  <- function( VCFfile,SNPrmv=NULL,INDrmv=NULL ) {

 cat( " > Reading VCF-file :\n    ",VCFfile,"\n",sep="" )
 inp.vcf <- readVcf( file=VCFfile )
 inp.dt  <- data.table( VAR=names(inp.vcf),REF=sub( "[/].*$","",sub( "^.*[_]","",names(inp.vcf) ) ),ALT=sub( "^.*[/]","",names(inp.vcf) ) )
 inp.dt[["TYPE"]]                                          <- "IND"
 inp.dt[ ( nchar(REF)==1 ) & ( nchar(ALT)==1 ) ][["TYPE"]] <- "SNP"

 snp.var <- c()
 if ( !is.null(SNPrmv) ) {

 cat( " > Reading SNP-file :\n    ",SNPrmv, "\n",sep="" )
 snp.vcf <- readVcf( file=SNPrmv  )

 cat( " > Variants from \"SNPrmv\" regions removed\n",sep="" )
 tmp.ovr <- findOverlaps( query=rowRanges(inp.vcf),subject=rowRanges(snp.vcf),ignore.strand=T )
 snp.var <- inp.dt[ unique(queryHits(tmp.ovr)) ][ TYPE=="SNP" ][["VAR"]]

 }

 ind.var <- c()
 if ( !is.null(INDrmv) ) {

 cat( " > Reading IND-file :\n    ",INDrmv, "\n",sep="" )
 ind.vcf <- readVcf( file=INDrmv  )

 cat( " > Variants from \"INDrmv\" regions removed\n",sep="" )
 tmp.ovr <- findOverlaps( query=rowRanges(inp.vcf),subject=rowRanges(ind.vcf),ignore.strand=T )
 ind.var <- inp.dt[ unique(queryHits(tmp.ovr)) ][ TYPE=="IND" ][["VAR"]]

 }

 rmv.var <- unique(c( snp.var,ind.var ))
 inp.vcf <- inp.vcf[ !( names(inp.vcf) %in% rmv.var ) ]

 cat( " > Done\n",sep="" )
 return( inp.vcf )

}

