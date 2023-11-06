#>> @pepap : load functions for variant filtering
source("04.mutation-analysis--filterVCF.R")

#>> @pepap : 
if ( !exists("NATIVE_VARS") ) {
 cat( " ** READING FVB_NJ & C57BL_6NJ native SNPs/INDELs **\n",sep="" )
 fvb.SNP.vcf <- readVcf("FVB_NJ.mgp.v5.snps.dbSNP142.vcf.gz")
 fvb.IND.vcf <- readVcf("FVB_NJ.mgp.v5.indels.dbSNP142.normed.vcf.gz")
 cbl.SNP.vcf <- readVcf("C57BL_6NJ.mgp.v5.snps.dbSNP142.vcf.gz")
 cbl.IND.vcf <- readVcf("C57BL_6NJ.mgp.v5.indels.dbSNP142.normed.vcf.gz")
 NATIVE_VARS <- "loaded"
 cat( "\n",sep="" )
}

VCFS="path/to/STRELKA/results/variants/variants.vcf.gz"
#>> @pepap : assign mouse-strain (C57BL_6NJ,FVB_NJ)
CONS <- c( "ID8","BR5" )

#>> @pepap : filter PASSED variants
if (T) {
cat( " ** READING STRELKA RESULTS **\n",sep="" )
mergedVars.dt <- data.table()
for ( i in seq_along(VCFS) ) {
 cat( "    -> ",CONS[i],"\n",sep="" )
 tmp.vcf <- filterVCF( VCFfile=paste0(STOR,VCFS[i]),FILTER.PASS=T,sele.gr=agi.ensembl.gr )
 cat( "       INPUT vs. FILTERED : ",length(tmp.vcf)," / ",sep="" )
 tmp.dt  <- data.table( VAR=names(tmp.vcf),REF=sub( "[/].*$","",sub( "^.*[_]","",names(tmp.vcf) ) ),ALT=sub( "^.*[/]","",names(tmp.vcf) ) )
 tmp.dt[["TYPE"]]                                          <- "IND"
 tmp.dt[ ( nchar(REF)==1 ) & ( nchar(ALT)==1 ) ][["TYPE"]] <- "SNP"

  if ( substr( x=CONS[i],start=1,stop=3 )=="ID8" ) {
   tmp.ovr <- findOverlaps( query=rowRanges(tmp.vcf),subject=rowRanges(cbl.SNP.vcf),ignore.strand=T )
   rmv.var <- tmp.dt[ unique(queryHits(tmp.ovr)) ][ TYPE=="SNP" ][["VAR"]]
   tmp.ovr <- findOverlaps( query=rowRanges(tmp.vcf),subject=rowRanges(cbl.IND.vcf),ignore.strand=T )
   rmv.var <- c( rmv.var,tmp.dt[ unique(queryHits(tmp.ovr)) ][ TYPE=="IND" ][["VAR"]] )
   tmp.vcf <- tmp.vcf[ !( names(tmp.vcf) %in% rmv.var ) ]
  }
  if ( substr( x=CONS[i],start=1,stop=3 )=="BR5" ) {
   tmp.ovr <- findOverlaps( query=rowRanges(tmp.vcf),subject=rowRanges(fvb.SNP.vcf),ignore.strand=T )
   rmv.var <- tmp.dt[ unique(queryHits(tmp.ovr)) ][ TYPE=="SNP" ][["VAR"]]
   tmp.ovr <- findOverlaps( query=rowRanges(tmp.vcf),subject=rowRanges(fvb.IND.vcf),ignore.strand=T )
   rmv.var <- c( rmv.var,tmp.dt[ unique(queryHits(tmp.ovr)) ][ TYPE=="IND" ][["VAR"]] )
   tmp.vcf <- tmp.vcf[ !( names(tmp.vcf) %in% rmv.var ) ]
  }

 cat( length(tmp.vcf)," (",substr( x=CONS[i],start=1,stop=3 ),")\n",sep="" )
 mrg.tmp.dt <-
  data.table(
   VAR    = names(tmp.vcf),
   CHR    = as.character(                  sub( "[:].*$","",names(tmp.vcf) )  ),
   CRD    = as.numeric(   sub( "^.*[:]","",sub( "[_].*$","",names(tmp.vcf) ) )),
   REF    = as.character( sub( "^.*[_]","",sub( "[/].*$","",names(tmp.vcf) ) )),
   ALT    = as.character( sub( "^.*[/]","",                 names(tmp.vcf) )  ),
   SAMPLE = CONS[i]
  )
 mergedVars.dt <- rbind( mergedVars.dt,mrg.tmp.dt )
 gc( verbose=T )
 cat( "    -> memory cleaning ...\n",sep="" )
 rm( list=c("tmp.vcf","tmp.dt","tmp.ovr","rmv.var","mrg.tmp.dt") )
 gc( verbose=T )
 cat( "\n",sep="" )
}
save( mergedVars.dt,file="mergedVars.dt.rda" )
} else {
load( "mergedVars.dt.rda",verbose=T )
}

#>> @pepap : filter variants from sequenced probes of ILLUMINA / AGILENT sequencing experiment
load( "targetReg.agi.mm10.gr.rda",verbose=T )
#load( "targetReg.Twist_Mouse_Exome_Target_Rev1_7APR20.mm10.gr.rda",verbose=T )
targetReg.gr            <- targetReg.agi.mm10.gr
TOTAL_TARR_WIDTH_MB     <- sum( width(targetReg.gr) )/1e06
PCDSS_TARR_WIDTH_MB     <- sum( width(targetReg.gr[ grepl( pattern="protein_coding",x=mcols(targetReg.gr)[["ANNOT"]] ) ]) )/1e06

tmp.bmrt.values      <- mergedVars.dt[ , unique(paste0(CHR,":",CRD,":",CRD)) ]
tmp.bmrt.attributes  <-
 c(
  "refsnp_id","refsnp_source","chr_name","chrom_start","allele","validated",
  "ensembl_transcript_stable_id","ensembl_gene_stable_id","ensembl_type",
  "consequence_type_tv","consequence_allele_string","ensembl_peptide_allele",
  "cdna_start","distance_to_transcript"
 )

#>> @pepap : annotate variants
library(biomaRt)

xFROM=1
rsid.dt <- data.table()
cat("\n")
for ( xTO in c(seq(from=1,to=length(tmp.bmrt.values)-1,by=10),length(tmp.bmrt.values))[-1] ) {
 cat( " -> from ",xFROM," to ",xTO," / ",length(tmp.bmrt.values),"\n",sep="" )
# print( seq(from=xFROM,to=xTO) )
 tmp.dt  <- as.data.table(getBM( attributes=tmp.bmrt.attributes,filters="chromosomal_region",values=tmp.bmrt.values[seq(from=xFROM,to=xTO)],mart=grcm38.mart ))
 rsid.dt <- rbind( rsid.dt,tmp.dt )
 xFROM <- ( xTO + 1 )
}

save( rsid.dt,file="rsid.dt.rda" )

#>> @pepap : load GENCODE gene annotation
load( "gencode.vM24.chr_patch_hapl_scaff.annotation.corrSeqLevels.trx2gene.PC.rda",verbose=T )
gencode.vM24.chr_patch_hapl_scaff.annotation.corrSeqLevels.trx2gene.PC[["tid_short"]] <-
 gencode.vM24.chr_patch_hapl_scaff.annotation.corrSeqLevels.trx2gene.PC[,sub( "[.].*$","",transcript_id )]

#load( "rsid.dt.rda",verbose=T )
rsid.dt[ is.na(ensembl_transcript_stable_id) ][["ensembl_transcript_stable_id"]] <- as.character("")
rsid.dt[ is.na(ensembl_gene_stable_id)       ][["ensembl_gene_stable_id"]]       <- as.character("")
rsid.dt[ is.na(ensembl_type)                 ][["ensembl_type"]]                 <- as.character("")
rsid.dt[ is.na(consequence_allele_string)    ][["consequence_allele_string"]]    <- as.character("")
rsid.dt[ is.na(consequence_type_tv)          ][["consequence_type_tv"]]          <- as.character("")
rsid.dt[ is.na(ensembl_peptide_allele)       ][["ensembl_peptide_allele"]]       <- as.character("")
rsid.dt[ ensembl_peptide_allele==T           ][["ensembl_peptide_allele"]]       <- as.character("T")
rsid.dt[["VAR"]]                                                                 <- rsid.dt[,paste0(chr_name,":",chrom_start,"_",consequence_allele_string)]

rsid.dt <-
 merge(
  rsid.dt,gencode.vM24.chr_patch_hapl_scaff.annotation.corrSeqLevels.trx2gene.PC[tag=="basic",c("tid_short","gene_name","trx_len"),with=F],
  by.x="ensembl_transcript_stable_id",by.y="tid_short",all.x=T,all.y=T,sort=F
 )
rsid.dt <- rsid.dt[ !is.na(trx_len) & !is.na(refsnp_id) ]

#>> @pepap : not interesting genes - lncRNAs, TEC transcripts ..
MANUALLY_REMOVED=c("ENSMUST00000177875","ENSMUST00000179982","ENSMUST00000179264")

mergedVars.annot.dt <- merge( mergedVars.dt,rsid.dt[ !( ensembl_transcript_stable_id %in% MANUALLY_REMOVED ) ],by="VAR",all.x=T,sort=F )
mergedVars.annot.dt[["TYPE"]]                                 <- "SNP"
mergedVars.annot.dt[ nchar(REF)>1 & nchar(ALT)<=1 ][["TYPE"]] <- "DEL"
mergedVars.annot.dt[ nchar(ALT)>1 & nchar(REF)<=1 ][["TYPE"]] <- "INS"
save( mergedVars.annot.dt,file="mergedVars.annot.dt.rda" )

#>> @pepap : report all samples, prepare TMB table
library(openxlsx)

xwb <- createWorkbook()
addWorksheet(   wb=xwb,sheetName="mutAnal--Appendix" )
writeDataTable( wb=xwb,sheet=    "mutAnal--Appendix", x=mergedVars.annot.dt[ !is.na(refsnp_id) ],   colNames=T,rowNames=F,headerStyle=createStyle(textDecoration="bold") )

TMB.dt <-
 data.table(
  condition         = c("SAMPLE1","SAMPLE2","SAMPLE3"),
  TMB.conding       = rep.int( as.numeric(""),3 ),
  TMB.nonsynonymous = rep.int( as.numeric(""),3 ),
  VAR.coding        = c(
   sum(!duplicated(mergedVars.annot.dt[ !is.na(refsnp_id) & SAMPLE=="SAMPLE1" ][["VAR"]])),
   sum(!duplicated(mergedVars.annot.dt[ !is.na(refsnp_id) & SAMPLE=="SAMPLE2" ][["VAR"]])),
   sum(!duplicated(mergedVars.annot.dt[ !is.na(refsnp_id) & SAMPLE=="SAMPLE3" ][["VAR"]]))
                       ),
  VAR.nonsynonymous = c(
   sum(!duplicated(mergedVars.annot.dt[ !is.na(refsnp_id) & consequence_type_tv!="synonymous_variant" & SAMPLE=="SAMPLE1" ][["VAR"]])),
   sum(!duplicated(mergedVars.annot.dt[ !is.na(refsnp_id) & consequence_type_tv!="synonymous_variant" & SAMPLE=="SAMPLE2" ][["VAR"]])),
   sum(!duplicated(mergedVars.annot.dt[ !is.na(refsnp_id) & consequence_type_tv!="synonymous_variant" & SAMPLE=="SAMPLE3" ][["VAR"]]))
                       ),
  TotalRegionInMB   = rep.int( c(TOTAL_TARR_WIDTH_MB),c(3) ),
  CodingRegionInMB  = rep.int( c(PCDSS_TARR_WIDTH_MB),c(3) )
 )
TMB.dt[["TMB.conding"]]       <- TMB.dt[ , VAR.coding/CodingRegionInMB        ]
TMB.dt[["TMB.nonsynonymous"]] <- TMB.dt[ , VAR.nonsynonymous/CodingRegionInMB ]

addWorksheet(   wb=xwb,sheetName="TMB" )
writeDataTable( wb=xwb,sheet=    "TMB", x=TMB.dt,   colNames=T,rowNames=F,headerStyle=createStyle(textDecoration="bold") )
saveWorkbook(   wb=xwb,file="MH--mutAnal--Appendix.20230728.xlsx",overwrite=T )

}

save( TMB.dt,file="TMB.dt.rda" )
