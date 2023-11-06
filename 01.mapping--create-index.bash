#!/bin/bash

exeCMD="STAR"
inpFA="GRCm38_68.fa"

#>> @pepap : generate STAR-index for the reference genome
${exeCMD} --runMode          genomeGenerate \
          --runThreadN       16             \
          --genomeDir        .              \
          --genomeFastaFiles ${inpFA}
