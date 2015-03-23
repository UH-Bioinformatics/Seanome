import os
import math
import yaml
import argparse
#YAML FORMAT
"""
libraries:
  IH: 
    - /full/path/to/IH_forward.fastq
    - /full/path/to/IH_reverse.fastq
  IGA: 
    - /full/path/to/IGA_forward.fastq
    - /full/path/to/IGA_reverse.fastq
    
"""

MAIN_SCRIPT="primary_script.sh"
CHILD_SCRIPT="job_%s_-_%s_-_child.sh"

def printPearLine(oscript):
   genericBlock(oscript, """pear -f '${input_forward_path}' -r '${input_reverse_path}' -o ${short_name} -j ${THREADS} """)


def printQualFilterAndMerge(oscript, threads):
   command = """fastq_quality_filter -q 20 -p 75 -i ${short_name}.%(type)s.fastq -o ${short_name}.%(type)s.cleaned.fastq -Q 33 """
   if threads >= 3:
      print >> oscript, """pidArr=()"""
      print >> oscript, """date; echo "START: Quality filter"""

      print >> oscript, """date; echo "%s" """%(command%dict(type="assembled"))
      print >> oscript, """%s &"""%(command%dict(type="assembled"))
      print >> oscript, """pidArr+=($!)"""

      print >> oscript, """date; echo "%s" """%(command%dict(type="unassembled.forward"))
      print >> oscript, """%s &"""%(command%dict(type="unassembled.forward"))
      print >> oscript, """pidArr+=($!)"""

      print >> oscript, """date; echo "%s" """%(command%dict(type="unassembled.reverse"))
      print >> oscript, """%s &"""%(command%dict(type="unassembled.reverse"))
      print >> oscript, """pidArr+=($!)"""
      print >> oscript, """wait ${pidArr[@]}"""
      print >> oscript, """date; echo "END: Quality filter"\n"""

   elif threads == 2:
      print >> oscript, """pidArr=()"""
      print >> oscript, """date; echo "START: Quality filter"""
      print >> oscript, """date; echo "%s" """%(command%dict(type="assembled"))
      print >> oscript, """%s &"""%(command%dict(type="assembled"))
      print >> oscript, """pidArr+=($!)"""

      print >> oscript, """date; echo "%s" """%(command%dict(type="unassembled.forward"))
      print >> oscript, """%s &"""%(command%dict(type="unassembled.forward"))
      print >> oscript, """pidArr+=($!)"""
      print >> oscript, """wait ${pidArr[@]}"""

      print >> oscript, """date; echo "%s" """%(command%dict(type="unassembled.reverse"))
      print >> oscript, """%s"""%(command%dict(type="unassembled.reverse"))
      print >> oscript, """date; echo "END: Quality filter"\n"""

   else:
      genericBlock(oscript, """%s"""%(command%dict(type="assembled")))
      genericBlock(oscript, """%s"""%(command%dict(type="unassembled.forward")))
      genericBlock(oscript, """%s"""%(command%dict(type="unassembled.reverse")))

   genericBlock(oscript, """cat ${short_name}.assembled.cleaned.fastq ${short_name}.unassembled.forward.cleaned.fastq  ${short_name}.unassembled.reverse.cleaned.fastq > ${short_name}.fastq""")
   print >> oscript, """merged_fastq='${short_name}.fastq'"""
   print >> oscript, """merged_fasta='${short_name}.fasta'"""


def renameFastq(oscript, addPrefix):
   if not addPrefix:
      print >> oscript, """date; echo "START fastq_rename ${merged_fastq} > ${merged_fastq}_ && mv ${merged_fastq}_ ${merged_fastq}" """
      print >> oscript, """fastq_rename ${merged_fastq} > ${merged_fastq}_ && mv ${merged_fastq}_ ${merged_fastq}"""
      print >> oscript, """date; echo "END fastq_rename ${merged_fastq} > ${merged_fastq}_ && mv ${merged_fastq}_ ${merged_fastq}" """
   else:
      print >> oscript, """date; echo "START fastq_rename ${merged_fastq} ${short_name} > ${merged_fastq}_ && mv ${merged_fastq}_ ${merged_fastq}" """
      print >> oscript, """fastq_rename ${merged_fastq} ${short_name} > ${merged_fastq}_ && mv ${merged_fastq}_ ${merged_fastq}"""
      print >> oscript, """date; echo "END fastq_rename ${merged_fastq} ${short_name} > ${merged_fastq}_ && mv ${merged_fastq}_ ${merged_fastq}" """


def genericBlock(oscript, line):
   print >> oscript, """date; echo "START %s" """%(line)
   print >> oscript, line
   print >> oscript, """date; echo "END %s"\n"""%(line)

def generateMulti(args):
   threads = max(1, int(math.floor( float(args.threads) / float(args.jobs))))
   parameters = yaml.load(open(args.config))
   children_scripts = []
   cid = 0
   passalongs = []		
   for ident, pairs in parameters['libraries'].iteritems():
      ident = ident.replace("_","-").replace(" ","-")
      scriptname = os.path.join(args.workdir, CHILD_SCRIPT%(cid, ident))
      cid += 1
      children_scripts.append(scriptname)
      with open(scriptname, "w") as oscript:
	 print >> oscript, "#!/bin/bash"
	 print >> oscript, "input_forward_path='%s'"%(pairs[0]) 
	 print >> oscript, "input_reverse_path='%s'"%(pairs[1]) 
	 print >> oscript, "short_name='%s'"%(ident) 
	 print >> oscript, "THREADS=%s\n"%(threads) 
	 print >> oscript, "mkdir ${short_name}"
	 print >> oscript, "cd ${short_name}"
	 printPearLine(oscript)
	 printQualFilterAndMerge(oscript, threads)
	 renameFastq(oscript, False)
	 genericBlock(oscript, """seqtk seq -A ${merged_fastq} > ${merged_fasta}""")
	 genericBlock(oscript, """vsearch -maskfasta ${merged_fasta} --hardmask --output ${short_name}.masked.fasta --threads ${THREADS}""")
	 genericBlock(oscript, """vsearch -sortbylength ${short_name}.masked.fasta --output ${short_name}.masked.sorted.fasta --threads ${THREADS}""")
	 genericBlock(oscript, """vsearch --cluster_smallmem ${short_name}.masked.sorted.fasta --strand plus --id 0.95  --consout ${short_name}_1.cons --msaout ${short_name}_1.msa --userout ${short_name}_1.out --userfields query+target+caln+qstrand+tstrand --mincols 80 --maxdiffs 10 --threads ${THREADS}""")
	 genericBlock(oscript, """muso.py  -t ${THREADS} -i1 ${short_name}_1.msa -o1 ${short_name}_mod_1.cons -i2 ${short_name}_1.out -o2 ${short_name}_mod_1.out -g 1 -m 3 -n 2000""")
         genericBlock(oscript, """vsearch -sortbylength ${short_name}_mod_1.cons --output ${short_name}_mod_1.sorted.fasta -threads ${THREADS}""")
         genericBlock(oscript, """vsearch --cluster_smallmem ${short_name}_mod_1.sorted.fasta --strand both --id 0.95  --consout ${short_name}_2.cons --msaout ${short_name}_2.msa --userout ${short_name}_2.out --userfields query+target+caln+qstrand+tstrand --mincols 80 --maxdiffs 10 --threads ${THREADS}""")
         genericBlock(oscript, """muso.py  -t ${THREADS} -i1 ${short_name}_2.msa -o1 ${short_name}.final.contigs -i2 ${short_name}_2.out -o2 ${short_name}.final.out -g 2 -m 3 -n 2000""")
         genericBlock(oscript, """trackOverlaps.py -i1 ${short_name}_mod_1.out  -i2 ${short_name}.final.out  -o ${short_name}.mapping_to_cons""")
         genericBlock(oscript, """vsearch -maskfasta ${short_name}.final.contigs -hardmask  -output ${short_name}.final.contigs.masked -qmask dust -threads ${THREADS}""")
         genericBlock(oscript, """sed -i -E '/>/!s/[acgt]/N/g' ${short_name}.final.contigs.masked""")
         genericBlock(oscript, """grep -v '>'  ${short_name}.final.contigs.masked |  sed -e '1i>${short_name}' | seqret stdin ${short_name}.temp.pseudo.fasta""")
         genericBlock(oscript, """build_lmer_table -l 12 -sequence ${short_name}.temp.pseudo.fasta -freq freqs_${short_name}""")
         genericBlock(oscript, """RepeatScout -sequence ${short_name}.temp.pseudo.fasta -output repeats_${short_name}.fa -freq freqs_${short_name} -l 12""")
         genericBlock(oscript, """usearch -usearch_local repeats_${short_name}.fa  -db ${short_name}.final.contigs.masked -id 0.80 -strand both -userfields query+target+id+qcov+tcov+qlo+qhi+tlo+thi -userout ${short_name}_bad.ids -query_cov 0.75  -threads ${THREADS} --maxaccepts 0 --maxrejects 0""")
         genericBlock(oscript, """maskSeqs.py -i ${short_name}_bad.ids  -f ${short_name}.final.contigs.masked -o ${short_name}.final.contigs.masked_ && mv ${short_name}.final.contigs.masked_ ${short_name}.final.contigs.masked""")
         genericBlock(oscript, """grep ">" ${short_name}.final.contigs.masked | sed 's/>//' |  sed 's/<unknown description>//' > ${short_name}_clean.ids""")
         genericBlock(oscript, """make_sam_with_cons.py -u  ${short_name}.mapping_to_cons -q ${short_name}.fastq -c ${short_name}_clean.ids -f ${short_name}.final.contigs.masked -l ${minClustSize} -m ${maxClustSize} multiple -o ${short_name}""")
	 print >> oscript, """cd .."""
      passalongs.append( (ident, "%s.bam"%(ident), "%s.bam,bai"%(ident), "%s_pseudo_ref_parts.fasta"%(ident), "%s_pseudo_ref.fasta"%(ident), )  )

   threads = args.threads
   with open(MAIN_SCRIPT, "w") as o:
      print >> o, "#!/bin/bash"
      print >> o , "OLDDIR=`pwd`"
      print >> o, "THREADS=%s"%(threads) 
      print >> o, "DB_NAME='%s'\n"%(args.database) 
      print >> o, "mkdir -p %s"%(args.workdir)	
      print >> o, "cd %s\n"%(args.workdir)	      
      if args.jobs != 1:
	 for c, cname in enumerate(children_scripts):
	    if c % args.jobs == 0:
	       if c != 0:
		  print >> o, "wait ${pidArr[@]}"
	       print >> o, """pidArr=()"""	
	    genericBlock(o, """bash %s &"""%(cname))
	 print >> o, "wait ${pidArr[@]}"
      else:	
	 for cname in children_scripts:
	    genericBlock(o, """bash %s """%(cname))
      print >> o, "mkidr csr"
      print >> o, "cd csr"
      print >> o, "mkdir cluster_bams"
      for d in passalongs:
	 genericBlock(o, """mv ../%(parent)s/%(bam)s ../%(parent)s/%(bamidx)s cluster_bams/"""%dict(parent = d[0], bam = d[1], bamidx = d[2]) )
      seedA, seedB= passalongs[:2]	
      genericBlock(o, """Seanome.py -t ${THREADS}  -d ${DB_NAME} seed_csr -i1 ../%(parentA)s/%(refA)s -n1 %(parentA)s -i2 ../%(parentB)s/%(refB)s -n2 %(parentB)s -l 150 -s 0.94"""%dict(parentA = seedA[0], refA=seedA[3], parentB = seedB[0], refB=seedB[3]))	
      for d in passalongs[2:]:
	 genericBlock(o, """Seanome.py -t ${THREADS}  -d ${DB_NAME} find_csr -g ../%(parent)s/%(ref)s -l 150 -s 0.94"""%dict(parent = d[0], ref = d[4]))			
      genericBlock(o, """Seanome.py -t ${THREADS} -d ${DB_NAME} consensus""")
      genericBlock(o, """Seanome.py -t ${THREADS} -d ${DB_NAME} inferSAM  -s cluster_bams""")
      genericBlock(o, """Seanome.py -t ${THREADS} -d ${DB_NAME} trimAL""")
      genericBlock(o, """Seanome.py -t ${THREADS} -d ${DB_NAME} cleanSAM""")
      genericBlock(o, """Seanome.py -t ${THREADS} -d ${DB_NAME} generateVCF""")
      genericBlock(o, """Seanome.py -t ${THREADS} -d ${DB_NAME} updateVCF""")
      print >> o, "cd ${OLDDIR}"
		


def generateSingle(args):
   threads = max(1, int(math.floor( float(args.threads) / float(args.jobs))))
   parameters = yaml.load(open(args.config))
   children_scripts = []
   cid = 0
   passalongs = []		
   for ident, pairs in parameters['libraries'].iteritems():
      ident = ident.replace("_","-").replace(" ","-")	
      scriptname = os.path.join(args.workdir, CHILD_SCRIPT%(cid, ident))
      cid += 1
      children_scripts.append(scriptname)
      with open(scriptname, "w") as oscript:
	 print >> oscript, "#!/bin/bash"
	 print >> oscript, "input_forward_path='%s'"%(pairs[0]) 
	 print >> oscript, "input_reverse_path='%s'"%(pairs[1]) 
	 print >> oscript, "short_name='%s'"%(ident) 
	 print >> oscript, "THREADS=%s\n"%(threads) 
	 print >> oscript, "mkdir ${short_name}"	
	 print >> oscript, "cd ${short_name}"	
	 printPearLine(oscript)
	 printQualFilterAndMerge(oscript, threads)
	 renameFastq(oscript, True)
	 print >> oscript, """cd .."""
      passalongs.append( ident  )
   threads = args.threads
   with open(MAIN_SCRIPT, "w") as o:
      print >> o, "#!/bin/bash"
      print >> o , "OLDDIR=`pwd`"
      print >> o, "THREADS=%s"%(threads) 
      print >> o, "DB_NAME='%s'\n"%(args.database) 
      print >> o, "mkdir -p %s"%(args.workdir)	
      print >> o, "cd %s"%(args.workdir)	
      if args.jobs != 1:
	 for c, cname in enumerate(children_scripts):
	    if c % args.jobs == 0:
	       if c != 0:
		  print >> o, "wait ${pidArr[@]}"
	       print >> o, """pidArr=()"""	
	    genericBlock(o, """bash %s &"""%(cname))
	 print >> o, "wait ${pidArr[@]}"
      else:	
	 for cname in children_scripts:
	    genericBlock(o, """bash %s """%(cname))
      print >> o, "mkidr csr"
      print >> o, "cd csr"
      print >> o, """combined_name="%s" """%("_".join(passalongs))
      genericBlock(o, """cat %s > ${combined_name}.fastq """%( " " .join( [ """ "../%(ident)s/%(ident)s.fastq" """%dict(ident=d) for d in passalongs] ) ) )
      genericBlock(o, """seqtk seq -A ${combined_name}.fastq > ${combined_name}.fasta""")
      genericBlock(o, """vsearch --derep_fulllength  ${combined_name}.fasta -minseqlength 1 -output ${combined_name}.dedup.fasta -uc ${combined_name}.uc --threads ${THREADS}""")
      genericBlock(o, """vsearch -maskfasta ${combined_name}.fasta --hardmask --output ${combined_name}.masked.fasta --threads ${THREADS}""")
      genericBlock(o, """vsearch -sortbylength ${combined_name}.masked.fasta --output ${combined_name}.masked.sorted.fasta --threads ${THREADS}""")
      genericBlock(o, """vsearch --cluster_smallmem ${combined_name}.masked.sorted.fasta --strand plus --id 0.95  --consout ${combined_name}_1.cons --msaout ${combined_name}_1.msa --userout ${combined_name}_1.out --userfields query+target+caln+qstrand+tstrand --mincols 80 --maxdiffs 10 --threads ${THREADS}""")
      genericBlock(o, """muso.py  -t ${THREADS} -i1 ${combined_name}_1.msa -o1 ${combined_name}_mod_1.cons -i2 ${combined_name}_1.out -o2 ${combined_name}_mod_1.out -g 1 -m 3 -n 2000""")
      genericBlock(o, """vsearch -sortbylength ${combined_name}_mod_1.cons --output ${combined_name}_mod_1.sorted.fasta -threads ${THREADS}""")
      genericBlock(o, """vsearch --cluster_smallmem ${combined_name}_mod_1.sorted.fasta --strand both --id 0.95  --consout ${combined_name}_2.cons --msaout ${combined_name}_2.msa --userout ${combined_name}_2.out --userfields query+target+caln+qstrand+tstrand --mincols 80 --maxdiffs 10 --threads ${THREADS}""")
      genericBlock(o, """muso.py -t ${THREADS} -i1 ${combined_name}_2.msa -o1 ${combined_name}.final.contigs -i2 ${combined_name}_2.out -o2 ${combined_name}.final.out -g 2 -m 3 -n 2000""")
      genericBlock(o, """trackOverlaps.py -i1 ${combined_name}_mod_1.out  -i2 ${combined_name}.final.out  -o ${combined_name}.mapping_to_cons""")
      genericBlock(o, """update_mapping.py -i ${combined_name}.uc -m ${combined_name}.mapping_to_cons""")
      genericBlock(o, """vsearch -maskfasta ${combined_name}.final.contigs -hardmask  -output ${combined_name}.final.contigs.masked -qmask dust""")
      genericBlock(o, """sed -i -E '/>/!s/[acgt]/N/g' ${combined_name}.final.contigs.masked""")
      genericBlock(o, """grep -v '>'  ${combined_name}.final.contigs.masked |  sed -e '1i>${combined_name}' | seqret stdin ${combined_name}.temp.pseudo.fasta""")
      genericBlock(o, """build_lmer_table -l 12 -sequence ${combined_name}.temp.pseudo.fasta -freq freqs_${combined_name}""")
      genericBlock(o, """RepeatScout -sequence ${combined_name}.temp.pseudo.fasta -output repeats_${combined_name}.fa -freq freqs_${combined_name} -l 12""")
      genericBlock(o, """usearch -usearch_local repeats_${combined_name}.fa  -db ${combined_name}.final.contigs.masked -id 0.80 -strand both -userfields query+target+id+qcov+tcov+qlo+qhi+tlo+thi -userout ${combined_name}_bad.ids -query_cov 0.75  -threads ${THREADS} --maxaccepts 0 --maxrejects 0 """)
      genericBlock(o, """maskSeqs.py -i ${combined_name}_bad.ids  -f ${combined_name}.final.contigs.masked -o ${combined_name}.final.contigs.masked_""")
      genericBlock(o, """mv ${combined_name}.final.contigs.masked_ ${combined_name}.final.contigs.masked""")
      genericBlock(o, """grep '>' ${combined_name}.final.contigs.masked | sed 's/>//' |  sed 's/<unknown description>//' > ${combined_name}_clean.ids""")
      genericBlock(o, """make_sam_with_cons.py -u  ${combined_name}.mapping_to_cons -q ${combined_name}.fastq -c ${combined_name}_clean.ids -f ${combined_name}.final.contigs.masked -l ${minClustSize} -m ${maxClustSize}  single -d DB_NAME""")
      genericBlock(o, """vcf_generator.py -t ${THREADS} -d  ${DB_NAME}""")
      genericBlock(o, """vcfmod.py -t ${THREADS} -d  ${DB_NAME}""")
      print >> o, "cd ${OLDDIR}"



if __name__ == "__main__":
	parser = argparse.ArgumentParser()

	parser.add_argument('-c', '--config', required = True, help = "configuration file in yaml format" )
	parser.add_argument('-d', '--database', required = True, help = "Name of the seanome sqlite database" )
	parser.add_argument('-t', '--threads',  required = False, default = 1, type = int, help = "total number of threads to utilize (default: 1)")
	parser.add_argument('-j', '--jobs',  required = False, default = 1, type = int, help = "Number of jobs to run in parallel.  This will divide the number of threads, and undersubscribe in case of uneven division (default: 1)")
	parser.add_argument('-w', '--workdir', required = False, default =".", help = "Output directory to run job in (deafult: current directory)")

	subparsers = parser.add_subparsers(dest='action', help='Available commands')

	parser_sub = subparsers.add_parser('multiple')
	parser_sub.set_defaults(func = generateMulti)

	parser_sub = subparsers.add_parser('single')
	parser_sub.set_defaults(func = generateSingle)
	args = parser.parse_args()
	args.database = os.path.basename(args.database)
	args.func(args)
