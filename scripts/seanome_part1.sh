#!/bin/bash

# TODO: Run from temp dir
# TODO: organize the directory as you'd expect it to be if run from the web.
# TODO : Delete files that are not necessary

: ${1?"Usage: $0 Forward Reverse Prefix threads"}
: ${2?"Usage: $0 Forward Reverse Prefix threads"}
: ${3?"Usage: $0 Forward Reverse Prefix threads"}
: ${4?"Usage: $0 Forward Reverse Prefix threads"}

input_for=${1}
input_rev=${2}
short_name=${3} 
THREADS=${4}

minClustSize=3
maxClustSize=200


cd $( dirname "${input_for}" )
echo "Running from: " && pwd 

echo "#1: Assembling pairs"
echo "pear -f $input_for -r $input_rev -o $short_name -j ${THREADS}"
pear -f $input_for -r $input_rev -o $short_name -j ${THREADS}
echo -e "\n\n"

pidArr=()
echo "#2: Quality filtering reads"
echo "/home/mahdi/bin/fastx/fastq_quality_filter -q 20 -p 75 -i ${short_name}.assembled.fastq -o ${short_name}.assembled.cleaned.fastq -Q 33"
echo "/home/mahdi/bin/fastx/fastq_quality_filter -q 20 -p 75 \
      -i ${short_name}.unassembled.forward.fastq -o ${short_name}.unassembled.cleaned.forward.fastq -Q 33"
echo "/home/mahdi/bin/fastx/fastq_quality_filter -q 20 -p 75 \
      -i ${short_name}.unassembled.reverse.fastq -o ${short_name}.unassembled.cleaned.reverse.fastq -Q 33"
fastq_quality_filter -q 20 -p 75 -i ${short_name}.assembled.fastq \
      -o ${short_name}.assembled.cleaned.fastq -Q 33 &
pidArr+=($!)
fastq_quality_filter -q 20 -p 75 \
      -i ${short_name}.unassembled.forward.fastq \
      -o ${short_name}.unassembled.cleaned.forward.fastq -Q 33 & 
pidArr+=($!)
fastq_quality_filter -q 20 -p 75 \
      -i ${short_name}.unassembled.reverse.fastq \
      -o ${short_name}.unassembled.cleaned.reverse.fastq -Q 33 &
pidArr+=($!)
wait  ${pidArr[@]}
echo -e "\n\n"

echo "#3: Generating input fastq file (renaming)"
echo "cat ${short_name}.assembled.cleaned.fastq ${short_name}.unassembled.cleaned.forward.fastq \
     ${short_name}.unassembled.cleaned.reverse.fastq > ${short_name}.fastq"
echo "fastq_rename ${short_name}.fastq > -o ${short_name}.fastq_"
cat ${short_name}.assembled.cleaned.fastq ${short_name}.unassembled.cleaned.forward.fastq \
     ${short_name}.unassembled.cleaned.reverse.fastq > ${short_name}.fastq
fastq_rename ${short_name}.fastq > ${short_name}.fastq_
mv ${short_name}.fastq_ ${short_name}.fastq
echo -e "\n\n"


echo "#4: Converting fastq sequences into fasta"
echo "seqtk seq -a ${short_name}.fastq > ${short_name}.fasta"
seqtk seq -A ${short_name}.fastq > ${short_name}.fasta
echo -e "\n\n"


echo "#5: Masking Fasta file"
echo "vsearch -maskfasta ${short_name}.fasta --hardmask --output ${short_name}.masked.fasta --threads ${THREADS}"
vsearch -maskfasta ${short_name}.fasta --hardmask \
 --output ${short_name}.masked.fasta --threads ${THREADS}
echo -e "\n\n"


echo "#6: Sorting Fasta file"
echo "vsearch -sortbylength ${short_name}.masked.fasta --output ${short_name}.masked.sorted.fasta --threads ${THREADS}"
vsearch -sortbylength ${short_name}.masked.fasta --output ${short_name}.masked.sorted.fasta --threads ${THREADS}
echo -e "\n\n"


echo "#7: Clustering-- First iteration (Plus strand)"
echo "vsearch --cluster_smallmem ${short_name}.masked.sorted.fasta \
      --strand plus --id 0.95  --consout ${short_name}_1.cons  --msaout ${short_name}_1.msa --userout ${short_name}_1.out \
      --userfields query+target+caln+qstrand+tstrand \
       --mincols 80 --maxdiffs 10 --threads ${THREADS}"
vsearch --cluster_smallmem ${short_name}.masked.sorted.fasta \
      --strand plus --id 0.95  --consout ${short_name}_1.cons \
      --msaout ${short_name}_1.msa --userout ${short_name}_1.out \
      --userfields query+target+caln+qstrand+tstrand \
      --mincols 80 --maxdiffs 10 --threads ${THREADS}
echo -e "\n\n"


echo "#8: Extending the consensus sequences and fixing the outfiles for 1st iteation Clustering"
echo "muso.py  -t ${THREADS} -i1 ${short_name}_1.msa -o1 ${short_name}_mod_1.cons -i2 ${short_name}_1.out -o2 ${short_name}_mod_1.out -g 1 -m 3 -n 2000 "
muso.py  -t ${THREADS} -i1 ${short_name}_1.msa -o1 ${short_name}_mod_1.cons \
 -i2 ${short_name}_1.out -o2 ${short_name}_mod_1.out -g 1 -m 3 -n 2000 
echo -e "\n\n"


echo "#9: sorting consensus sequence file"
echo "vsearch -sortbylength ${short_name}_mod_1.cons --output ${short_name}_mod_1.sorted.fasta -threads ${THREADS}"
vsearch -sortbylength ${short_name}_mod_1.cons \
 --output ${short_name}_mod_1.sorted.fasta -threads ${THREADS}
echo -e "\n\n"


echo "#10: Clustering-- First iteration (both strands)"
echo "vsearch --cluster_smallmem ${short_name}_mod_1.sorted.fasta \
  --strand both --id 0.95  --consout ${short_name}_2.cons  --msaout ${short_name}_2.msa --userout ${short_name}_2.out \
  --userfields query+target+caln+qstrand+tstrand \
  --mincols 80 --maxdiffs 10 --threads ${THREADS}"
vsearch --cluster_smallmem ${short_name}_mod_1.sorted.fasta \
  --strand both --id 0.95  --consout ${short_name}_2.cons  \
  --msaout ${short_name}_2.msa --userout ${short_name}_2.out \
  --userfields query+target+caln+qstrand+tstrand \
  --mincols 80 --maxdiffs 10 --threads ${THREADS}
echo -e "\n\n"


echo "#11: Extending the consensus sequences and fixing the outfiles for 2nd iteation Clustering"
echo "muso.py  -t ${THREADS} -i1 ${short_name}_2.msa -o1 ${short_name}.final.contigs -i2 ${short_name}_2.out -o2 ${short_name}.final.out -g 2 -m 3 -n 2000 "
muso.py  -t ${THREADS} -i1 ${short_name}_2.msa -o1 ${short_name}.final.contigs \
 -i2 ${short_name}_2.out -o2 ${short_name}.final.out -g 2 -m 3 -n 2000 
echo -e "\n\n"


echo "#12: Tracking overlaps by following indirection"
echo "trackOverlaps.py -i1 ${short_name}_mod_1.out  -i2 ${short_name}.final.out  -o ${short_name}.mapping_to_cons"
trackOverlaps.py -i1 ${short_name}_mod_1.out  -i2 ${short_name}.final.out  -o ${short_name}.mapping_to_cons

echo -e "\n\n"


echo "#13: Dusting the contigs generated from the assembly"
echo "vsearch -maskfasta ${short_name}.final.contigs \
      -hardmask  -output ${short_name}.final.contigs.masked -qmask dust -threads ${THREADS}"
echo "sed -i -E '/>/!s/[acgt]/N/g' ${short_name}.final.contigs.masked"
vsearch -maskfasta ${short_name}.final.contigs \
      -hardmask  -output ${short_name}.final.contigs.masked -qmask dust -threads ${THREADS}
sed -i -E '/>/!s/[acgt]/N/g' ${short_name}.final.contigs.masked
echo -e "\n\n"


echo "#14: Concatenating the contigs into a temporary pseudo-genome"
echo "grep -v '>'  ${short_name}.final.contigs.masked |  sed -e '1i>${short_name}' | seqret stdin ${short_name}.temp.pseudo.fasta"
grep -v '>'  ${short_name}.final.contigs.masked |  sed -e '1i>${short_name}' | seqret stdin ${short_name}.temp.pseudo.fasta
echo -e "\n\n"


echo "#15: Finding repeats in the pseudogenome"
echo "build_lmer_table -l 12 -sequence ${short_name}.temp.pseudo.fasta -freq freqs_${short_name}"
echo "RepeatScout -sequence ${short_name}.temp.pseudo.fasta -output repeats_${short_name}.fa -freq freqs_${short_name} -l 12"
echo "usearch -usearch_local repeats_${short_name}.fa  -db ${short_name}.final.contigs.masked -id 0.80 -strand both \
 -userfields query+target+id+qcov+tcov+qlo+qhi+tlo+thi -userout ${short_name}_bad.ids \
 -query_cov 0.75  -threads ${THREADS} --maxaccepts 0 --maxrejects 0"
echo "maskSeqs.py -i ${short_name}_bad.ids  -f ${short_name}.final.contigs.masked -o ${short_name}.final.contigs.masked_"
echo "mv ${short_name}.final.contigs.masked_ ${short_name}.final.contigs.masked"
echo "grep "\>" ${short_name}.final.contigs.masked | sed 's/>//' |  sed 's/<unknown description>//' > ${short_name}_clean.ids"
build_lmer_table -l 12 -sequence ${short_name}.temp.pseudo.fasta -freq freqs_${short_name}
RepeatScout -sequence ${short_name}.temp.pseudo.fasta -output repeats_${short_name}.fa -freq freqs_${short_name} -l 12
# find matches with repeats in the contigs. Hits need to be of at least 75% similarity
usearch -usearch_local repeats_${short_name}.fa  -db ${short_name}.final.contigs.masked -id 0.80 -strand both \
 -userfields query+target+id+qcov+tcov+qlo+qhi+tlo+thi -userout ${short_name}_bad.ids \
 -query_cov 0.75  -threads ${THREADS} --maxaccepts 0 --maxrejects 0
maskSeqs.py -i ${short_name}_bad.ids  -f ${short_name}.final.contigs.masked -o ${short_name}.final.contigs.masked_
mv ${short_name}.final.contigs.masked_ ${short_name}.final.contigs.masked
grep ">" ${short_name}.final.contigs.masked | sed 's/>//' |  sed 's/<unknown description>//' > ${short_name}_clean.ids
echo -e "\n\n"


#this enerates the _pseudo_ref.fasta files, as well as the corresponding bam file
echo "#16: Building consensus sequences with min $minClustSize and max $maxClustSize"
echo "make_sam_with_cons.py \
      -u  ${short_name}.mapping_to_cons -q ${short_name}.fastq -c ${short_name}_clean.ids \
      -f ${short_name}.final.contigs.masked -l ${minClustSize} -m ${maxClustSize} multiple -o ${short_name}"
make_sam_with_cons.py \
  -u  ${short_name}.mapping_to_cons -q ${short_name}.fastq -c ${short_name}_clean.ids \
  -f ${short_name}.final.contigs.masked -l ${minClustSize} -m ${maxClustSize} multiple -o ${short_name}
echo -e "\n\n"
