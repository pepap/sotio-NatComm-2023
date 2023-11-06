#!/bin/bash

exeCMD="/storage/brno1-cerit/home/pepap/pLIBRARY/STAR-2.7.0c/bin/Linux_x86_64_static/STAR"
inpFA="/storage/brno3-cerit/home/pepap/Annotations/SNP-InDels-mice-strains/GRCm38_68.fa"

${exeCMD} --runMode          genomeGenerate \
          --runThreadN       16             \
          --genomeDir        .              \
          --genomeFastaFiles ${inpFA}
