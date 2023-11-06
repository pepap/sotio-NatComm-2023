#!/bin/bash

#>> @pepap : download reference genome
wget --show-progress ftp://ftp-mouse.sanger.ac.uk/ref/GRCm38_68.fa.gz
wget --show-progress ftp://ftp-mouse.sanger.ac.uk/ref/GRCm38_68.fa.fai

FTP_STOR="ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/strain_specific_vcfs"
SUFF0001=".mgp.v5.indels.dbSNP142.normed.vcf.gz"
SUFF0002=".mgp.v5.indels.dbSNP142.normed.vcf.gz.md5"
SUFF0003=".mgp.v5.indels.dbSNP142.normed.vcf.gz.tbi"
SUFF0004=".mgp.v5.indels.dbSNP142.normed.vcf.gz.tbi.md5"
SUFF0005=".mgp.v5.snps.dbSNP142.vcf.gz"
SUFF0006=".mgp.v5.snps.dbSNP142.vcf.gz.md5"
SUFF0007=".mgp.v5.snps.dbSNP142.vcf.gz.tbi"
SUFF0008=".mgp.v5.snps.dbSNP142.vcf.gz.tbi.md5"

#>> @pepap : download SNPs & INDELs for selected mice strains
MSTRAINS=( BALB_cJ C3H_HeJ C57BL_6NJ FVB_NJ )

for mstr in `echo ${MSTRAINS[@]}`

do echo ${mstr}

wget --show-progress ${FTP_STOR}/${mstr}${SUFF0001}
wget --show-progress ${FTP_STOR}/${mstr}${SUFF0002}
wget --show-progress ${FTP_STOR}/${mstr}${SUFF0003}
wget --show-progress ${FTP_STOR}/${mstr}${SUFF0004}

wget --show-progress ${FTP_STOR}/${mstr}${SUFF0005}
wget --show-progress ${FTP_STOR}/${mstr}${SUFF0006}
wget --show-progress ${FTP_STOR}/${mstr}${SUFF0007}
wget --show-progress ${FTP_STOR}/${mstr}${SUFF0008}

echo ""

done
