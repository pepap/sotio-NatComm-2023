#!/bin/bash

module add strelka-2.9.10

REFA="GRCm38_68.fa"
BAM="path/to/aligned/BAM/file"

#>> @pepap : prepare & configure germline STRELKA run
configureStrelkaGermlineWorkflow.py --bam ${BAM} \
                                    --referenceFasta ${REFA} \
                                    --exome                  \
                                    --runDir .
#>> @pepap : run STRELKA mutation analysis
./runWorkflow.py -m local -j 4
