#!/bin/bash

: ${1?"Usage: $0 short_nameA short_namebB databasename threads"}
: ${2?"Usage: $0 short_nameA short_namebB databasename threads"}
: ${3?"Usage: $0 short_nameA short_namebB databasename threads"}
: ${4?"Usage: $0 short_nameA short_namebB databasename threads"}

short_nameA=${1}
short_nameB=${2}
dbname=${3}.db3
THREADS=${4}

echo "#1 creating directory csr to do work in and change to that directory"
echo "mkdir csr"
echo "cd csr"
mkdir csr 
cd csr


### Move the bam files genearted in the last step of the run.sh
echo "#2 copy bam files from each library"
echo "mkdir cluster_bams"
echo "cp ../${short_nameA}/*.bam* cluster_bams/"
echo "cp ../${short_nameB}/*.bam* cluster_bams/"
mkdir cluster_bams
cp ../${short_nameA}/*.bam* cluster_bams/
cp ../${short_nameB}/*.bam* cluster_bams/
echo -e "\n\n"


### move the pseudo_ref files (concat and parts) taht were generated in the run.sh
echo "#3 copy pseudo refrences from each library"
echo "mkdir pseudo_refs"
echo "cp ../${short_nameA}/*_pseudo_ref*.fasta pseudo_refs/"
echo "cp ../${short_nameB}/*_pseudo_ref*.fasta pseudo_refs/"
mkdir pseudo_refs
cp ../${short_nameA}/*_pseudo_ref*.fasta pseudo_refs/
cp ../${short_nameB}/*_pseudo_ref*.fasta pseudo_refs/
echo -e "\n\n"


echo "#4 find seed CSRs (lastz replacement)"
echo "usearch -usearch_local pseudo_refs/${short_nameA}_pseudo_ref_parts.fasta  \
 -db pseudo_refs/${short_nameB}_pseudo_ref_parts.fasta -id 0.80 -strand both \
 -userfields query+target+id+alnlen+bits+qstrand+tstrand+qlo+qhi+tlo+thi+qrow+trow \
 -userout  ${short_nameA}_${short_nameB}.out -threads ${THREADS} -evalue 1e-5 "
echo "sort -k 1 -k4nr ${short_nameA}_${short_nameB}.out > sorted_hits"
echo "getSeed.py sorted_hits pseudo_refs/${short_nameA}_pseudo_ref_parts.fasta \
 pseudo_refs/${short_nameB}_pseudo_ref_parts.fasta > ${short_nameA}_${short_nameB}.lastz"
echo "awk '{$1=A; $2=B; print $0}' A="${short_nameA}" B="${short_nameB}" \
 ${short_nameA}_${short_nameB}.lastz > ${short_nameA}_${short_nameB}.lastz_"
echo "mv ${short_nameA}_${short_nameB}.lastz_ ${short_nameA}_${short_nameB}.lastz "
usearch -usearch_local pseudo_refs/${short_nameA}_pseudo_ref_parts.fasta  \
 -db pseudo_refs/${short_nameB}_pseudo_ref_parts.fasta -id 0.80 -strand both \
 -userfields query+target+id+alnlen+bits+qstrand+tstrand+qlo+qhi+tlo+thi+qrow+trow \
 -userout  ${short_nameA}_${short_nameB}.out -threads ${THREADS} -evalue 1e-5 
sort -k 1 -k4nr ${short_nameA}_${short_nameB}.out > sorted_hits
getSeed.py sorted_hits pseudo_refs/${short_nameA}_pseudo_ref_parts.fasta \
 pseudo_refs/${short_nameB}_pseudo_ref_parts.fasta > ${short_nameA}_${short_nameB}.lastz
awk '{$1=A; $2=B; print $0}' A="${short_nameA}" B="${short_nameB}" \
 ${short_nameA}_${short_nameB}.lastz > ${short_nameA}_${short_nameB}.lastz_
mv ${short_nameA}_${short_nameB}.lastz_ ${short_nameA}_${short_nameB}.lastz 
echo -e "\n\n"

# This step also will generate the database file.
echo "#5 Filter and insert seed CSRs into sqlite database"
echo "Seanome.py -t ${THREADS} seed_csr -i ${short_nameA}_${short_nameB}.lastz -d ${dbname} -l 150 -s 0.94"
Seanome.py -t ${THREADS} seed_csr -i ${short_nameA}_${short_nameB}.lastz -d ${dbname} -l 150 -s 0.94
echo -e "\n\n"

#### ####
## find_csr step would go here in a run with more than 2 species/groups.
#### ####

echo "#6 Create initial consensus sequences"
echo "makeCons.py -t ${THREADS} -d  ${dbname}"
makeCons.py -t ${THREADS} -d  ${dbname}
echo -e "\n\n"

echo "#7 based on the CSRs, generate a sam/bam file for each.  Also, update the consensus"
echo "Seanome.py -t ${THREADS} inferSAM  -s cluster_bams/ -d  ${dbname}"
Seanome.py -t ${THREADS} inferSAM  -s cluster_bams/ -d  ${dbname}
echo -e "\n\n"

echo "#8 With the CSRs, run trimal and store the cleanned sequences and the colnumber log"
echo "trimall.py -t ${THREADS} -d  ${dbname}"
trimall.py -t ${THREADS} -d  ${dbname}
echo -e "\n\n"

echo "#9 store modified sam files based on what the trimall.py step provided"
echo "cleanSAM.py -t ${THREADS} -d  ${dbname}"
cleanSAM.py -t ${THREADS} -d  ${dbname}
echo -e "\n\n"

echo "#10 generate a vcf file for each CSR "
echo "vcf_generator.py -t ${THREADS} -d  ${dbname}"
vcf_generator.py -t ${THREADS} -d  ${dbname}
echo -e "\n\n"

echo "#11 Manually go through the pile at each location the vcf specified, and store a modified vcf file"
echo "vcfmod.py -t ${THREADS} -d  ${dbname}"
vcfmod.py -t ${THREADS} -d  ${dbname}
echo -e "\n\n"

cd ..