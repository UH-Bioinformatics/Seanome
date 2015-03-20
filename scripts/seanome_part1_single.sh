#!/bin/bash
#SBATCH -J TestPL
#SBATCH --ntasks=10
#SBATCH -t 2-16:00
#SBATCH -p lm.q


# TODO: Run from temp dir
# TODO: organize the directory as you'd expect it to be if run from the web.
# TODO : Delete files that are not necessary

: ${1?"Usage: $0 Forward Reverse Prefix threads"}
: ${2?"Usage: $0 Forward Reverse Prefix threads"}
: ${3?"Usage: $0 Forward Reverse Prefix threads"}
: ${4?"Usage: $0 Forward Reverse Prefix threads"}


#get scripts location and add it to the path... 
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ) 
export PATH=${DIR}:${PATH}


input_for=${1}
input_rev=${2}
short_name=${3} 
THREADS=${4}

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
fastq_quality_filter -q 20 -p 75 -i ${short_name}.assembled.fastq -o ${short_name}.assembled.cleaned.fastq -Q 33 &
pidArr+=($!)
fastq_quality_filter -q 20 -p 75 \
      -i ${short_name}.unassembled.forward.fastq -o ${short_name}.unassembled.cleaned.forward.fastq -Q 33 & 
pidArr+=($!)
fastq_quality_filter -q 20 -p 75 \
      -i ${short_name}.unassembled.reverse.fastq -o ${short_name}.unassembled.cleaned.reverse.fastq -Q 33 &
pidArr+=($!)
wait  ${pidArr[@]}
echo -e "\n\n"


echo "#5: Generating input fastq file (renaming)"
echo "cat ${short_name}.assembled.cleaned.fastq ${short_name}.unassembled.cleaned.forward.fastq \
     ${short_name}.unassembled.cleaned.reverse.fastq > ${short_name}.fastq"
echo "/home/mahdi/bin/fastx/fastx_renamer -n COUNT -i IH.fastq -o IH.fastq_ -Q 33"
cat ${short_name}.assembled.cleaned.fastq ${short_name}.unassembled.cleaned.forward.fastq \
     ${short_name}.unassembled.cleaned.reverse.fastq > ${short_name}.fastq
fastx_renamer -n COUNT -i ${short_name}.fastq -o ${short_name}.fastq_ -Q 33
mv ${short_name}.fastq_ ${short_name}.fastq
echo -e "\n\n"


# Stop here and to the whole combine file step.
# prefix all the sequence ids with the short_name
# TODO: DLS -- Need to write my own in C
renamer.py ${short_name}.fastq ${short_name} ${short_name}.rename.fastq
