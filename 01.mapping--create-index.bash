#!/bin/bash

exeCMD="STAR"
inpFA="GRCm38_68.fa"

${exeCMD} --runMode          genomeGenerate \
          --runThreadN       16             \
          --genomeDir        .              \
          --genomeFastaFiles ${inpFA}
