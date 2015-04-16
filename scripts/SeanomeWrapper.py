#!/usr/bin/python
import sys
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

findcsr:
  - minlen: 150
  - minsim: 0.94
   
"""

MAIN_SCRIPT="""primary_script.sh"""
SAMPLE_SCRIPT="""%s_sample_%s_step_%s.sh"""
SAMPLE_RUNNER="""sample_runner_script_%s.sh"""

def printBashHeader(o):
   print >> o, "#!/bin/bash"


def bashargsparse(o, mlen = 150, msim = 0.94, onlyparams = False):  
   print >> o, "minlen=%s"%(mlen)
   print >> o, "minsim=%s"%(msim)
   print >> o, "minClustSize=3"
   print >> o, "maxClustSize=200"
   if onlyparams:
      return
   print >> o, """while getopts ":c:z:l:s:" opt; do"""
   print >> o, """case $opt in"""
   print >> o, """c) minClustSize="${OPTARG}" """
   print >> o, """;;"""
   print >> o, """z) maxClustSize="${OPTARG}" """
   print >> o, """;;"""
   print >> o, """l) minlen="${OPTARG}" """
   print >> o, """;;"""
   print >> o, """s) minsim="${OPTARG}" """
   print >> o, """;;"""
   print >> o, """\?)"""
   print >> o, """echo "Invalid option: -${OPTARG}" >&2"""
   print >> o, """exit 1"""
   print >> o, """;;"""
   print >> o, """:)"""
   print >> o, """echo "Option -${OPTARG} requires an argument." >&2"""
   print >> o, """exit 1"""
   print >> o, """;;"""
   print >> o, """esac"""
   print >> o, """done"""


def genericBlock(oscript, line):
   print >> oscript, """echo -e "\\n"; date; echo -e "START %s" """%(line.replace('"','\\"'))
   print >> oscript, line
   print >> oscript, """echo -e "\\n"; date; echo -e "END %s\\n" """%(line.replace('"','\\"'))


def printPearLine(oscript):
   genericBlock(oscript, """pear -f ${input_forward_path} -r ${input_reverse_path} -o ${NAME} -j ${THREADS} """)


def filterCMDWithPID(oscript, command):
      print >> oscript, """echo -e "\\n"; date; echo -e "%s" """%(command)
      print >> oscript, """%s &"""%(command)
      print >> oscript, """pidArr+=($!)"""


def printQualFilterAndMerge(oscript, threads):
   command = """fastq_quality_filter -q 20 -p 75 -i ${NAME}.%(type)s.fastq -o ${NAME}.%(type)s.cleaned.fastq -Q 33 """
   if threads >= 3:
      print >> oscript, """pidArr=()"""
      print >> oscript, """echo -e "\\n"; date; echo -e "START: Quality filter" """
      filterCMDWithPID(oscript, command%dict(type="assembled"))
      filterCMDWithPID(oscript, command%dict(type="unassembled.forward"))
      filterCMDWithPID(oscript, command%dict(type="unassembled.reverse"))
      print >> oscript, """wait ${pidArr[@]}"""
      print >> oscript, """echo -e "\\n"; date; echo -e "END: Quality filter"\\n"""
   elif threads == 2:
      print >> oscript, """pidArr=()"""
      print >> oscript, """echo -e "\\n"; date; echo -e "START: Quality filter" """
      filterCMDWithPID(oscript, command%dict(type="assembled"))
      filterCMDWithPID(oscript, command%dict(type="unassembled.forward"))
      print >> oscript, """wait ${pidArr[@]}"""
      filterCMDWithPID(oscript, command%dict(type="unassembled.reverse"))
      print >> oscript, """wait ${pidArr[@]}"""
      print >> oscript, """echo -e "\\n"; date; echo -e "END: Quality filter"\\n"""
   else:
      genericBlock(oscript, """%s"""%(command%dict(type="assembled")))
      genericBlock(oscript, """%s"""%(command%dict(type="unassembled.forward")))
      genericBlock(oscript, """%s"""%(command%dict(type="unassembled.reverse")))

   genericBlock(oscript, """cat ${NAME}.assembled.cleaned.fastq ${NAME}.unassembled.forward.cleaned.fastq  ${NAME}.unassembled.reverse.cleaned.fastq > ${NAME}.fastq""")
   print >> oscript, """merged_fastq="${NAME}.fastq" """
   print >> oscript, """merged_fasta="${NAME}.fasta" """


def renameFastq(oscript, addSuffix):
   if not addSuffix:
      print >> oscript, """echo -e "\\n"; date; echo -e "START fastq_rename ${merged_fastq} > ${merged_fastq}_ && mv ${merged_fastq}_ ${merged_fastq}" """
      print >> oscript, """fastq_rename ${merged_fastq} > ${merged_fastq}_ && mv ${merged_fastq}_ ${merged_fastq}"""
      print >> oscript, """echo -e "\\n"; date; echo -e "END fastq_rename ${merged_fastq} > ${merged_fastq}_ && mv ${merged_fastq}_ ${merged_fastq}" """
   else:
      print >> oscript, """echo -e "\\n"; date; echo -e "START fastq_rename ${merged_fastq} ${NAME} > ${merged_fastq}_ && mv ${merged_fastq}_ ${merged_fastq}" """
      print >> oscript, """fastq_rename ${merged_fastq} ${NAME} > ${merged_fastq}_ && mv ${merged_fastq}_ ${merged_fastq}"""
      print >> oscript, """echo -e "\\n"; date; echo -e "END fastq_rename ${merged_fastq} ${NAME} > ${merged_fastq}_ && mv ${merged_fastq}_ ${merged_fastq}" """


def buildCommonStepOne(args, oscript, pairs, ident, threads, prefix = False):
   printBashHeader(oscript)
   print >> oscript, "input_forward_path='%s'"%(pairs[0]) 
   print >> oscript, "input_reverse_path='%s'"%(pairs[1]) 
   print >> oscript, "NAME='%s'"%(ident) 
   print >> oscript, "THREADS=%s\n"%(threads) 
   print >> oscript, "mkdir ${NAME}"
   print >> oscript, "cd ${NAME}"
   printPearLine(oscript)
   printQualFilterAndMerge(oscript, threads)
   renameFastq(oscript, prefix)


def commonMainScript(o, args, threads, clustered, basic = False):
   printBashHeader(o)
   print >> o , "OLDDIR=`pwd`"
   print >> o, "THREADS=%s"%(threads) 
   print >> o, "DB_NAME='%s'\n"%(args.database) 
   if not basic:
      print >> o, "mkdir -p %s"%(args.workdir)	
      print >> o, "cd %s\n"%(args.workdir)	      
      print >> o, "mkdir csr"
      if clustered:
         print >> o, "mkdir csr/cluster_bams"


def advanceNotice(o, args, name):      
   minclust = "${minClustSize}"
   maxclust = "${maxClustSize}"
   if args.advance:
      print >> o, """echo -e "Enter min cluster coverage for library %s, followed by [ENTER]:" """%(name)
      print >> o, """read minsize"""
      print >> o, """echo -e "Enter max cluster coverage for library %s, followed by [ENTER]:" """%(name)
      print >> o, """read maxsize"""
      minclust = "${minsize}"
      maxclust = "${maxsize}"
   return minclust, maxclust


def buildMakeSamScript(childname, ident, threads, runJelly = True):
   with open(childname, "w") as oscript:
      printBashHeader(oscript)
      bashargsparse(oscript)
      print >> oscript , """cd "%s" """%(ident)
      genericBlock(oscript, """make_sam_with_cons.py -u  %(short)s.mapping_to_cons -q  %(short)s.fastq -c  %(short)s_clean.ids -f  %(short)s.final.contigs.masked -l ${minClustSize} -m ${maxClustSize} multiple -o "%(short)s" """%dict(short = ident) )
      if runJelly:
         genericBlock(oscript, """jellyfish count -m 21 -s 512M -t %(threads)s -C %(ref)s -o %(short)s.jelly"""%dict(threads=threads, ref="%s_pseudo_ref_parts.fasta"%(ident), short=ident ) )
         genericBlock(oscript, """jellyfish dump  %(short)s.jelly > %(short)s.kmer.counts; rm %(short)s.jelly"""%dict(short=ident) )
         genericBlock(oscript, """filter_kmer_counts.py %(short)s.kmer.counts %(short)s 2; rm %(short)s.kmer.counts %(short)s_kmer.filter.tmp"""%dict(short=ident))
      print >> oscript, """cd .."""
      if runJelly:
         return ["%s.bam"%(ident), "%s.bam.bai"%(ident), "%s_pseudo_ref_parts.fasta"%(ident), "%s_pseudo_ref.fasta"%(ident), "%s_kmer.filter"%(ident) ]
      else:
         return ["%s.bam"%(ident), "%s.bam.bai"%(ident), "%s_pseudo_ref_parts.fasta"%(ident), "%s_pseudo_ref.fasta"%(ident), None ]


def buildChildRunner(args, children_scripts, passalongs, multi, threads = 1, largesingle = False):
   o = open(SAMPLE_RUNNER%(1), "w")
   printBashHeader(o)
   if multi:
      oo = open(SAMPLE_RUNNER%(2), "w") 
      printBashHeader(oo)
      bashargsparse(oo, onlyparams = args.advance)
   if args.jobs != 1:
      proc = []
      for c, cname in enumerate(children_scripts):
         if c % args.jobs == 0:
            if c != 0:
               print >> o, "wait ${pidArr[@]}"
               proc = []
            print >> o, """pidArr=()"""	
         genericBlock(o, """bash %s &\npidArr+=($!)"""%(cname))
         proc.append(1)
      if proc:
         print >> o, "wait ${pidArr[@]}"
      if multi:
         proc = []
         for c, d in enumerate(passalongs):
            childname = SAMPLE_SCRIPT%(d[0], c, 2)
            d.extend(buildMakeSamScript(childname, d[0], threads, not largesingle))
            if c % args.jobs == 0:
               if c != 0:
                  print >> o, "wait ${pidArr[@]}"
                  proc = []
               print >> o, """pidArr=()"""	
            minclust, maxclust = advanceNotice(oo, args, d[0])
            genericBlock(oo, """bash %s -c %s -z %s &\npidArr+=($!)"""%(childname, minclust, maxclust))
            proc.append(1)
         if proc:
            print >> o, "wait ${pidArr[@]}"
   else:	
      for cname in children_scripts:
         genericBlock(o, """bash %s"""%(cname))
      if multi:
         for c, d in enumerate(passalongs):
            childname = SAMPLE_SCRIPT%(d[0], c, 2)
            d.extend(buildMakeSamScript(childname, d[0], threads, not largesingle))
            minclust, maxclust = advanceNotice(oo, args, d[0])
            genericBlock(oo, """bash %s -c %s -z %s"""%(childname, minclust, maxclust))
   o.close()
   if multi:
      oo.close()
   return passalongs


def clusterSection(oscript):
   genericBlock(oscript, """vsearch -sortbylength ${NAME}.masked.fasta --output ${NAME}.masked.sorted.fasta --threads ${THREADS}""")
   genericBlock(oscript, """vsearch --cluster_smallmem ${NAME}.masked.sorted.fasta --strand plus --id 0.95 --msaout ${NAME}_1.msa --userout ${NAME}_1.out --userfields query+target+caln+qstrand+tstrand --mincols 80 --maxdiffs 10 --threads ${THREADS}""")
   genericBlock(oscript, """muso.py  -t ${THREADS} -i1 ${NAME}_1.msa -o1 ${NAME}_mod_1.cons -i2 ${NAME}_1.out -o2 ${NAME}_mod_1.out -g 1 -m 3 -n 2000""")
   genericBlock(oscript, """vsearch -sortbylength ${NAME}_mod_1.cons --output ${NAME}_mod_1.sorted.fasta -threads ${THREADS}""")
   genericBlock(oscript, """vsearch --cluster_smallmem ${NAME}_mod_1.sorted.fasta --strand both --id 0.95 --msaout ${NAME}_2.msa --userout ${NAME}_2.out --userfields query+target+caln+qstrand+tstrand --mincols 80 --maxdiffs 10 --threads ${THREADS}""")
   genericBlock(oscript, """muso.py  -t ${THREADS} -i1 ${NAME}_2.msa -o1 ${NAME}.final.contigs -i2 ${NAME}_2.out -o2 ${NAME}.final.out -g 2 -m 3 -n 2000""")
   genericBlock(oscript, """trackOverlaps.py -i1 ${NAME}_mod_1.out  -i2 ${NAME}.final.out  -o ${NAME}.mapping_to_cons""")


def maskingSection(oscript, skipmasking = False):
   if not skipmasking: 
      genericBlock(oscript, """vsearch -maskfasta ${NAME}.final.contigs -hardmask  -output ${NAME}.final.contigs.masked -qmask dust -threads ${THREADS}""")
      genericBlock(oscript, """sed -i -E '/>/!s/[acgt]/N/g' ${NAME}.final.contigs.masked""")
      genericBlock(oscript, """grep -v '>'  ${NAME}.final.contigs.masked |  sed -e "1i>${NAME}" | seqret stdin ${NAME}.temp.pseudo.fasta""")
      genericBlock(oscript, """build_lmer_table -l 12 -sequence ${NAME}.temp.pseudo.fasta -freq freqs_${NAME}""")
      genericBlock(oscript, """RepeatScout -sequence ${NAME}.temp.pseudo.fasta -output repeats_${NAME}.fa -freq freqs_${NAME} -l 12""")
      genericBlock(oscript, """usearch -usearch_local repeats_${NAME}.fa  -db ${NAME}.final.contigs.masked -id 0.80 -strand both -userfields query+target+id+qcov+tcov+qlo+qhi+tlo+thi -userout ${NAME}_bad.ids -query_cov 0.75  -threads ${THREADS} --maxaccepts 0 --maxrejects 0""")
      genericBlock(oscript, """maskSeqs.py -i ${NAME}_bad.ids  -f ${NAME}.final.contigs.masked -o ${NAME}.final.contigs.masked_ && mv ${NAME}.final.contigs.masked_ ${NAME}.final.contigs.masked""")
   else:
      genericBlock(oscript, """cp ${NAME}.final.contigs ${NAME}.final.contigs.masked""")
   genericBlock(oscript, """grep ">" ${NAME}.final.contigs.masked | sed 's/>//' |  sed 's/<unknown description>//' > ${NAME}_clean.ids""")         


def singleMergedInput(args, threads, passalongs):
   with open(SAMPLE_RUNNER%(2), "w") as o:
      commonMainScript(o, args, threads, False)
      print >> o, "cd csr"
      print >> o, """NAME="combined" """
      genericBlock(o, """cat %s > ${NAME}.fastq """%( " " .join( [ """ "../%(ident)s/%(ident)s.fastq" """%dict(ident=d[0]) for d in passalongs] ) ) )
      genericBlock(o, """seqtk seq -A ${NAME}.fastq > ${NAME}.fasta""")
      genericBlock(o, """vsearch --derep_fulllength  ${NAME}.fasta -minseqlength 1 -output ${NAME}.dedup.fasta -uc ${NAME}.uc --threads ${THREADS}""")
      genericBlock(o, """vsearch -maskfasta ${NAME}.dedup.fasta --hardmask --output ${NAME}.masked.fasta --threads ${THREADS}""")
      clusterSection(o)
      genericBlock(o, """update_mapping.py -i ${NAME}.uc -m ${NAME}.mapping_to_cons""")
      maskingSection(o, args.skipmasking)
      genericBlock(o, """coverageInformation.py -s "${NAME}"  -c ${NAME}_clean.ids -m ${NAME}.mapping_to_cons""")
   return "combined"


def orderingFunction(oo, passalongs, dummyorder = False):
   """
   Initially generate sometype of similarity measure between the different samples.  With this information,
   we need to get it so that we can work with it dynamically in the bash script.  Initally we populate
   two assosiative arrays with information for each sample.  We then iterate line by line, over the output
   from the python script.  this provides us with a processing order.  We use the associative arrays to get us the correct
   input files to use as we move through the samples, in the defined ordering.
   """

   if not dummyorder:
      genericBlock(oo, """findBestOrder_quick.py kmer_counts _kmer.filter search.order""")
   else:
      genericBlock(oo, """echo -e "%s" > search.order """%("\\n".join([s[0] for s in passalongs])) )

   print >> oo, "declare -A mapping_full;"
   print >> oo, "declare -A mapping_part;"
   for d in passalongs:
      print >> oo, """mapping_part["%s"]="%s" """%(d[0], "../%s/%s"%(d[0], d[3]))
      print >> oo, """mapping_full["%s"]="%s" """%(d[0], "../%s/%s"%(d[0], d[4]))              
   print >> oo, """ordering=()"""
   print >> oo, """i=0"""
   print >> oo, """while read line"""
   print >> oo, """do"""
   print >> oo, """ ordering[$i]=${line}"""
   print >> oo, """ i=$(($i+1))"""
   print >> oo, """done < search.order"""
   genericBlock(oo, """Seanome.py -t ${THREADS} -d ${DB_NAME} seed_csr -i1 ${mapping_part[${ordering[0]}]} -n1 ${ordering[0]} -i2 ${mapping_part[${ordering[1]}][0]} -n2 ${ordering[1]}  -l ${minlen} -s ${minsim} """)
   if len(passalongs) > 2:
      print >> oo, """for i in {2..%s}"""%((len(passalongs)-1))
      print >> oo, """do"""
      genericBlock(oo, """ Seanome.py -t ${THREADS} -d ${DB_NAME} find_csr -g ${mapping_full[${ordering[$i]}][1]} -l ${minlen} -s ${minsim} """%dict(parent = d[0], ref = d[4]))
      print >> oo, """done"""


def sampleRunner3_Multi(args, parameters, threads, passalongs, largesingle = False):
   with open(SAMPLE_RUNNER%(3), "w") as oo:
      commonMainScript(oo, args, threads, True)
      bashargsparse(oo,mlen = parameters.get('findcsr', {}).get('minlen', "150") , msim = parameters.get('findcsr', {}).get('minsim', "0.94"))
      print >> oo, "cd csr"
      for d in passalongs:
         genericBlock(oo, """mv ../%(parent)s/%(bam)s ../%(parent)s/%(bamidx)s cluster_bams/"""%dict(parent = d[0], bam = d[1], bamidx = d[2]) )
            
      print >> oo, """mkdir kmer_counts"""
      for d in passalongs:
         genericBlock(oo, """mv ../%(parent)s/%(cnt)s kmer_counts/"""%dict(parent = d[0], cnt = d[5]) )         
      if len(passalongs) >= 2:
         orderingFunction(oo, passalongs, largesingle)
      else:
         print >> sys.stderr, "Only 1 library is present.. Require at least 2 libraries!"
         sys.exit(0)
      genericBlock(oo, """Seanome.py -t ${THREADS} -d ${DB_NAME} consensus""")
      genericBlock(oo, """Seanome.py -t ${THREADS} -d ${DB_NAME} inferSAM  -s cluster_bams""")
      genericBlock(oo, """Seanome.py -t ${THREADS} -d ${DB_NAME} trimAL""")
      genericBlock(oo, """Seanome.py -t ${THREADS} -d ${DB_NAME} cleanSAM""")
      genericBlock(oo, """Seanome.py -t ${THREADS} -d ${DB_NAME} generateVCF""")
      genericBlock(oo, """Seanome.py -t ${THREADS} -d ${DB_NAME} updateVCF""")
      print >> oo, """cd .."""


def buildSampleScripts(args, parameters, threads, single):
   children_scripts = []
   cid = 0
   passalongs = []		
   for ident, pairs in parameters['libraries'].iteritems():
      ident = ident.replace("_","-").replace(" ","-")
      scriptname = os.path.join(args.workdir, SAMPLE_SCRIPT%(ident, cid, 1))
      cid += 1
      children_scripts.append(scriptname)
      with open(scriptname, "w") as oscript:
         passalongs.append([ident])       
         buildCommonStepOne(args, oscript, pairs, ident, threads, single)
         if not single:
            genericBlock(oscript, """seqtk seq -A ${merged_fastq} > ${merged_fasta}""")
            genericBlock(oscript, """vsearch -maskfasta ${merged_fasta} --hardmask --output ${NAME}.masked.fasta --threads ${THREADS}""")
            clusterSection(oscript)
            maskingSection(oscript, args.skipmasking)
            genericBlock(oscript, """coverageInformation.py -s "${NAME}"  -c ${NAME}_clean.ids -m ${NAME}.mapping_to_cons""")
         print >> oscript, """cd .."""
   return children_scripts, passalongs


def generateMulti(args, largesingle = False):
   threads = max(1, int(math.floor( float(args.threads) / float(args.jobs))))
   parameters = yaml.load(open(args.config))
   children_scripts, passalongs = buildSampleScripts(args, parameters, threads, False)
   buildChildRunner(args, children_scripts, passalongs, True, threads, largesingle)
   threads = args.threads   
   with open(MAIN_SCRIPT, "w") as o:
      commonMainScript(o, args, threads, True)  
      print >> o, "bash %s"%(SAMPLE_RUNNER%(1))
      print >> o, "bash %s"%(SAMPLE_RUNNER%(2))
      print >> o, "bash %s"%(SAMPLE_RUNNER%(3))
      print >> o, "cd ${OLDDIR}"
   sampleRunner3_Multi(args, parameters, threads, passalongs, largesingle)
		

def sampleRunner3_Single(args, parameters, threads, comboname):
   with open(SAMPLE_RUNNER%(3), "w") as oo:
      commonMainScript(oo, args, threads, False)
      print >> oo, "cd csr"
      print >> oo, """NAME="%s" """%(comboname)
      bashargsparse(oo,mlen = parameters.get('findcsr', {}).get('minlen', "150") , msim = parameters.get('findcsr', {}).get('minsim', "0.94"))
      minclust, maxclust = advanceNotice(oo, args, "${NAME}")   
      genericBlock(oo, """make_sam_with_cons.py -u  ${NAME}.mapping_to_cons -q ${NAME}.fastq -c ${NAME}_clean.ids -f ${NAME}.final.contigs.masked -l %s -m %s single -d ${DB_NAME} -t ${THREADS} -s ${minlen} """%(minclust, maxclust))
      genericBlock(oo, """Seanome.py -t ${THREADS} -d ${DB_NAME} generateVCF""")
      genericBlock(oo, """Seanome.py -t ${THREADS} -d ${DB_NAME} updateVCF""")
      print >> oo, """cd .."""


def generateSingle(args):
   threads = max(1, int(math.floor( float(args.threads) / float(args.jobs))))
   parameters = yaml.load(open(args.config))
   if len(parameters['libraries']) > 2:
      generateMulti(args, largesingle = True)
      return
   children_scripts, passalongs = buildSampleScripts(args, parameters, threads, True)
   buildChildRunner(args, children_scripts, passalongs, False, threads)
   threads = args.threads
   with open(MAIN_SCRIPT, "w") as o:
      commonMainScript(o, args, threads, False)
      print >> o, "bash %s"%(SAMPLE_RUNNER%(1))
      print >> o, "bash %s"%(SAMPLE_RUNNER%(2))
      print >> o, "bash %s"%(SAMPLE_RUNNER%(3))
      print >> o, "cd ${OLDDIR}"
   comboname = singleMergedInput(args, threads, passalongs)
   sampleRunner3_Single(args, parameters, threads, comboname)


if __name__ == "__main__":
	parser = argparse.ArgumentParser()

	parser.add_argument('-c', '--config', required = True, help = "configuration file in yaml format" )
	parser.add_argument('-d', '--database', required = True, help = "Name of the seanome sqlite database" )
	parser.add_argument('-t', '--threads',  required = False, default = 1, type = int, help = "total number of threads to utilize (default: 1)")
	parser.add_argument('-j', '--jobs',  required = False, default = 1, type = int, help = "Number of jobs to run in parallel.  This will divide the number of threads, and undersubscribe in case of uneven division (default: 1)")
	#parser.add_argument('-w', '--workdir', required = False, default =".", help = "Output directory to run job in (deafult: current directory)")
        #parser.add_argument("-a", "--advance", action = "store_true", required = False, help = "generates coverage information and requires the use to provide input")
        parser.add_argument("-s", "--skipmasking", action = "store_true", required = False, help = "Disable repeat masking")
        parser.set_defaults(workdir=".")
        parser.set_defaults(advance=False)

	subparsers = parser.add_subparsers(dest='action', help='Available commands')
	parser_sub = subparsers.add_parser('multiple')
	parser_sub.set_defaults(func = generateMulti)
	parser_sub = subparsers.add_parser('single')
	parser_sub.set_defaults(func = generateSingle)

	args = parser.parse_args()
	args.database = os.path.basename(args.database)
	args.func(args)
