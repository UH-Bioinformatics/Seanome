#!/bin/bash
#SBATCH -J TestPL
#SBATCH --ntasks=10
#SBATCH -t 2-16:00
#SBATCH -p lm.q

: ${1?"Usage: $0 path_to_.fastq-1 path_to_.fastq-2 Prefix threads"}
: ${2?"Usage: $0 path_to_.fastq-1 path_to_.fastq-2 Prefix threads"}
: ${3?"Usage: $0 path_to_.fastq-1 path_to_.fastq-2 Prefix threads"}
: ${4?"Usage: $0 path_to_.fastq-1 path_to_.fastq-2 Prefix threasd"}

#get scripts location and add it to the path... 
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ) 
export PATH=${DIR}:${PATH}

combined_name=${3} 
dbname=${combined_name}.db3
THREADS=${4}

minClustSize=3
maxClustSize=200


cd $( dirname "${input_for}" )
echo "Running from: " && pwd 


# combine all the renamed fastq files from the part 1 scrip
echo "#2: Merge the two libraries into one"
echo "cat ${1} ${2}   > ${combined_name}.fastq"
cat ${1} ${2} > ${combined_name}.fastq
echo -e "\n\n"


echo "#2: Converting fastq sequences into fasta"
echo "seqtk seq -a ${combined_name}.fastq > ${combined_name}.fasta"
seqtk seq -A ${combined_name}.fastq > ${combined_name}.fasta
echo -e "\n\n"


echo "#3: Converting fastq sequences into fasta"
echo "vsearch --derep_fulllength  ${combined_name}.fasta -minseqlength 1 -output ${combined_name}.dedup.fasta -uc ${combined_name}.uc  --threads ${THREADS}"
vsearch --derep_fulllength  ${combined_name}.fasta -minseqlength 1 -output ${combined_name}.dedup.fasta -uc ${combined_name}.uc --threads ${THREADS}
echo -e "\n\n"


echo "#4: Masking Fasta file"
echo "vsearch -maskfasta ${combined_name}.fasta --hardmask --output ${combined_name}.masked.fasta --threads ${THREADS}"
vsearch -maskfasta ${combined_name}.fasta --hardmask \
 --output ${combined_name}.masked.fasta --threads ${THREADS}
echo -e "\n\n"


echo "#5: Sorting Fasta file"
echo "vsearch -sortbylength ${combined_name}.masked.fasta --output ${combined_name}.masked.sorted.fasta --threads ${THREADS}"
vsearch -sortbylength ${combined_name}.masked.fasta --output ${combined_name}.masked.sorted.fasta --threads ${THREADS}
echo -e "\n\n"


echo "#6: Clustering-- First iteration (Plus strand)"
echo "vsearch --cluster_smallmem ${combined_name}.masked.sorted.fasta \                                                                                                            
      --strand plus --id 0.95  --consout ${combined_name}_1.cons  --msaout ${combined_name}_1.msa --userout ${combined_name}_1.out \                                                     
      --userfields query+target+caln+qstrand+tstrand \                                                                                                                          
       --mincols 80 --maxdiffs 10 --threads ${THREADS}"
vsearch --cluster_smallmem ${combined_name}.masked.sorted.fasta \
      --strand plus --id 0.95  --consout ${combined_name}_1.cons \
      --msaout ${combined_name}_1.msa --userout ${combined_name}_1.out \
      --userfields query+target+caln+qstrand+tstrand \
      --mincols 80 --maxdiffs 10 --threads ${THREADS}
echo -e "\n\n"


echo "#7: Extending the consensus sequences and fixing the outfiles for 1st iteation Clustering"
echo "muso.py  -t ${THREADS} -i1 ${combined_name}_1.msa -o1 ${combined_name}_mod_1.cons -i2 ${combined_name}_1.out -o2 ${combined_name}_mod_1.out -g 1 -m 3 -n 2000 "
muso.py  -t ${THREADS} -i1 ${combined_name}_1.msa -o1 ${combined_name}_mod_1.cons \
 -i2 ${combined_name}_1.out -o2 ${combined_name}_mod_1.out -g 1 -m 3 -n 2000
echo -e "\n\n"


echo "#8: sorting consensus sequence file"
echo "vsearch -sortbylength ${combined_name}_mod_1.cons --output ${combined_name}_mod_1.sorted.fasta -threads ${THREADS}"
vsearch -sortbylength ${combined_name}_mod_1.cons \
 --output ${combined_name}_mod_1.sorted.fasta -threads ${THREADS}
echo -e "\n\n"


echo "#9: Clustering-- First iteration (both strands)"
echo "vsearch --cluster_smallmem ${combined_name}_mod_1.sorted.fasta \                                                                                                             
  --strand both --id 0.95  --consout ${combined_name}_2.cons  --msaout ${combined_name}_2.msa --userout ${combined_name}_2.out \                                                         
  --userfields query+target+caln+qstrand+tstrand \                                                                                                                              
  --mincols 80 --maxdiffs 10 --threads ${THREADS}"
vsearch --cluster_smallmem ${combined_name}_mod_1.sorted.fasta \
  --strand both --id 0.95  --consout ${combined_name}_2.cons  \
  --msaout ${combined_name}_2.msa --userout ${combined_name}_2.out \
  --userfields query+target+caln+qstrand+tstrand \
  --mincols 80 --maxdiffs 10 --threads ${THREADS}
echo -e "\n\n"


echo "#10: Extending the consensus sequences and fixing the outfiles for 2nd iteation Clustering"
echo "muso.py  -t ${THREADS} -i1 ${combined_name}_2.msa -o1 ${combined_name}.final.contigs -i2 ${combined_name}_2.out -o2 ${combined_name}.final.out -g 2 -m 3 -n 2000 "
muso.py -t ${THREADS} -i1 ${combined_name}_2.msa -o1 ${combined_name}.final.contigs \
 -i2 ${combined_name}_2.out -o2 ${combined_name}.final.out -g 2 -m 3 -n 2000
echo -e "\n\n"


echo "#11: Tracking overlaps by following indirection"
echo "trackOverlaps.py -i1 ${combined_name}_mod_1.out  -i2 ${combined_name}.final.out  -o ${combined_name}.mapping_to_cons"
trackOverlaps.py -i1 ${combined_name}_mod_1.out  -i2 ${combined_name}.final.out  -o ${combined_name}.mapping_to_cons
echo -e "\n\n"


# Insert step to update the mapping_to_cons to contain all the dereplicated hits.
echo "#12: Update the mapping with deduplicated sequences"
echo "update_mapping.py -i ${combined_name}.uc -m ${combined_name}.mapping_to_cons"
update_mapping.py -i ${combined_name}.uc -m ${combined_name}.mapping_to_cons


echo "#13: Dusting the contigs generated from the assembly"
echo "usearch -maskfasta ${combined_name}.final.contigs -hardmask  -output ${combined_name}.final.contigs.masked -qmask dust"
echo "sed -i 's/[acgt]/N/g' ${combined_name}.final.contigs.masked"
vsearch -maskfasta ${combined_name}.final.contigs \
      -hardmask  -output ${combined_name}.final.contigs.masked -qmask dust
sed -i -E '/>/!s/[acgt]/N/g' ${combined_name}.final.contigs.masked
echo -e "\n\n"


echo "#14: Concatenating the contigs into a temporary pseudo-genome"
echo "grep -v '>'  ${combined_name}.final.contigs.masked |  sed -e '1i>${combined_name}' | seqret stdin ${combined_name}.temp.pseudo.fasta"
grep -v '>'  ${combined_name}.final.contigs.masked |  sed -e '1i>${combined_name}' | seqret stdin ${combined_name}.temp.pseudo.fasta
echo -e "\n\n"


echo "#15: Finding repeats in the pseudogenome"
echo "build_lmer_table -l 12 -sequence ${combined_name}.temp.pseudo.fasta -freq freqs_${combined_name}"
echo "RepeatScout -sequence ${combined_name}.temp.pseudo.fasta -output repeats_${combined_name}.fa -freq freqs_${combined_name} -l 12"
echo "usearch -usearch_local repeats_${combined_name}.fa  -db ${combined_name}.final.contigs.masked -id 0.80 -strand both \                                                           
 -userfields query+target+id+qcov+tcov+qlo+qhi+tlo+thi -userout ${combined_name}_bad.ids \                                                                                         
 -query_cov 0.75  -threads ${THREADS} --maxaccepts 0 --maxrejects 0"
echo "maskSeqs.py -i ${combined_name}_bad.ids  -f ${combined_name}.final.contigs.masked -o ${combined_name}.final.contigs.masked_"
echo "mv ${combined_name}.final.contigs.masked_ ${combined_name}.final.contigs.masked"
echo "grep '>' ${combined_name}.final.contigs.masked | sed 's/>//' |  sed 's/<unknown description>//' > ${combined_name}_clean.ids"
build_lmer_table -l 12 -sequence ${combined_name}.temp.pseudo.fasta -freq freqs_${combined_name}
RepeatScout -sequence ${combined_name}.temp.pseudo.fasta -output repeats_${combined_name}.fa -freq freqs_${combined_name} -l 12
# find matches with repeats in the contigs. Hits need to be of at least 75% similarity                                                                                          
usearch -usearch_local repeats_${combined_name}.fa  -db ${combined_name}.final.contigs.masked -id 0.80 -strand both \
 -userfields query+target+id+qcov+tcov+qlo+qhi+tlo+thi -userout ${combined_name}_bad.ids \
 -query_cov 0.75  -threads ${THREADS} --maxaccepts 0 --maxrejects 0 
maskSeqs.py -i ${combined_name}_bad.ids  -f ${combined_name}.final.contigs.masked -o ${combined_name}.final.contigs.masked_
mv ${combined_name}.final.contigs.masked_ ${combined_name}.final.contigs.masked
grep '>' ${combined_name}.final.contigs.masked | sed 's/>//' |  sed 's/<unknown description>//' > ${combined_name}_clean.ids
echo -e "\n\n"


#Takes fastq which are the assembled and non-assembled sequences as input
#we need the Fastq because we need the qualities of the reads to generate the sams
#this also generates the _pseudo_ref.fasta
echo "#16: Building consensus sequences with min $minClustSize and max $maxClustSize"
echo "make_sam_with_cons.py -u  ${combined_name}.mapping_to_cons -q ${combined_name}.fastq \
 -c ${combined_name}_clean.ids -f ${combined_name}.final.contigs.masked -l ${minClustSize} \ 
 -m ${maxClustSize}  single -d ${combined_name}.db3"
make_sam_with_cons.py -u  ${combined_name}.mapping_to_cons -q ${combined_name}.fastq \
-c ${combined_name}_clean.ids -f ${combined_name}.final.contigs.masked \
-l ${minClustSize} -m ${maxClustSize}  single -d ${combined_name}.db3
echo -e "\n\n"


echo "#17 generate a vcf file for each CSR "
echo "vcf_generator.py -t ${THREADS} -d  ${dbname}"
vcf_generator.py -t ${THREADS} -d  ${dbname}
echo -e "\n\n"


echo "#18 Manually go through the pile at each location the vcf specified, and store a modified vcf file"
echo "vcfmod.py -t ${THREADS} -d  ${dbname}"
vcfmod.py -t ${THREADS} -d  ${dbname}
echo -e "\n\n"

