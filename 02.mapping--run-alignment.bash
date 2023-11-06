#!/bin/bash

exeCMD="STAR"
IND="path/to/previously/created/STAR_index"
FASTQS1="path/to/FASTQ_read_sequences"
FASTQS2="path/to/FASTQ_mate_sequences"
OUTBAM="output_BAM_name"

#>> @pepap : run an alignment
${exeCMD} --runMode alignReads \
          --runThreadN  4 \
          --genomeDir ${IND} \
          --genomeLoad LoadAndRemove \
          --readFilesIn ${FASTQS1} ${FASTQS2} \
          --readFilesCommand zcat \
          --readStrand Unstranded \
          --limitBAMsortRAM  20000000000 \
          --outFileNamePrefix ${OUTBAM}.pe. \
          --outReadsUnmapped Fastx \
          --outQSconversionAdd +32 \
          --outSAMtype BAM SortedByCoordinate \
          --outFilterMultimapNmax 99999 \
          --outFilterMismatchNoverLmax 0.1 \
          --outFilterMatchNminOverLread 0.66 \
          --alignSJoverhangMin   10 \
          --alignSJDBoverhangMin 10 \
          --chimOutType Junctions \
          --chimSegmentMin 20

#>> @pepap : create an index file
samtools index ${OUTBAM}.pe.Aligned.sortedByCoord.out.bam
