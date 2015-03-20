#!/bin/bash

samfile=${1}
bamfile=${1%.sam}

samtools view -bS ${samfile} | samtools sort - ${bamfile}
samtools index ${bamfile}.bam