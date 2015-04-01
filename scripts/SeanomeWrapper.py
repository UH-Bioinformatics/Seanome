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
    
"""

MAIN_SCRIPT="""primary_script.sh"""
CHILD_SCRIPT="""%s_child_%s_step_%s.sh"""
CHILD_RUNNER="""child_runner_script_%s.sh"""


def clusterSizeHdr(o):
   print >> o, "minClustSize=3"
   print >> o, "maxClustSize=200"
   print >> o, """if [ -z "$1" ]; then """
   print >> o, "minClustSize=3"
   print >> o, "else"
   print >> o, "minClustSize=${1}"
   print >> o, """fi"""
   print >> o, """if [ -z "$2" ]; then """
   print >> o, "maxClustSize=200"
   print >> o, "else"
   print >> o, "maxClustSize=${2}"
   print >> o, """fi"""


def minlenAndminSim(o, parameters, shift = False):
   if shift:
      larg = "$3"
      sarg = "$4"
   else:
      larg = "$1"
      sarg = "$2"

   mlen = parameters.get('findcsr', {}).get('minlen', "150")
   msim = parameters.get('findcsr', {}).get('minsim', "0.94")
   print >> o, "minlen=%s"%(mlen)
   print >> o, "minsim=%s"%(msim)
   print >> o, """if [ -z "%s" ]; then """%(larg)
   print >> o, "minlen=%s"%(mlen)
   print >> o, "else"
   print >> o, "minlen=%s"%(larg)
   print >> o, """fi"""
   print >> o, """if [ -z "%s" ]; then """%(sarg)
   print >> o, "minsim=%s"%(msim)
   print >> o, "else"
   print >> o, "minsim=%s"%(sarg)
   print >> o, """fi"""


def printPearLine(oscript):
   genericBlock(oscript, """pear -f ${input_forward_path} -r ${input_reverse_path} -o ${NAME} -j ${THREADS} """)


def printQualFilterAndMerge(oscript, threads):
   command = """fastq_quality_filter -q 20 -p 75 -i ${NAME}.%(type)s.fastq -o ${NAME}.%(type)s.cleaned.fastq -Q 33 """
   if threads >= 3:
      print >> oscript, """pidArr=()"""
      print >> oscript, """echo -e "\\n"; date; echo -e "START: Quality filter" """

      print >> oscript, """echo -e "\\n"; date; echo -e "%s" """%(command%dict(type="assembled"))
      print >> oscript, """%s &"""%(command%dict(type="assembled"))
      print >> oscript, """pidArr+=($!)"""

      print >> oscript, """echo -e "\\n"; date; echo -e "%s" """%(command%dict(type="unassembled.forward"))
      print >> oscript, """%s &"""%(command%dict(type="unassembled.forward"))
      print >> oscript, """pidArr+=($!)"""

      print >> oscript, """echo -e "\\n"; date; echo -e "%s" """%(command%dict(type="unassembled.reverse"))
      print >> oscript, """%s &"""%(command%dict(type="unassembled.reverse"))
      print >> oscript, """pidArr+=($!)"""
      print >> oscript, """wait ${pidArr[@]}"""
      print >> oscript, """echo -e "\\n"; date; echo -e "END: Quality filter"\\n"""
   elif threads == 2:
      print >> oscript, """pidArr=()"""
      print >> oscript, """echo -e "\\n"; date; echo -e "START: Quality filter" """
      print >> oscript, """echo -e "\\n"; date; echo -e "%s" """%(command%dict(type="assembled"))
      print >> oscript, """%s &"""%(command%dict(type="assembled"))
      print >> oscript, """pidArr+=($!)"""

      print >> oscript, """echo -e "\\n"; date; echo -e "%s" """%(command%dict(type="unassembled.forward"))
      print >> oscript, """%s &"""%(command%dict(type="unassembled.forward"))
      print >> oscript, """pidArr+=($!)"""
      print >> oscript, """wait ${pidArr[@]}"""

      print >> oscript, """echo -e "\\n"; date; echo -e "%s" """%(command%dict(type="unassembled.reverse"))
      print >> oscript, """%s"""%(command%dict(type="unassembled.reverse"))
      print >> oscript, """echo -e "\\n"; date; echo -e "END: Quality filter"\\n"""
   else:
      genericBlock(oscript, """%s"""%(command%dict(type="assembled")))
      genericBlock(oscript, """%s"""%(command%dict(type="unassembled.forward")))
      genericBlock(oscript, """%s"""%(command%dict(type="unassembled.reverse")))

   genericBlock(oscript, """cat ${NAME}.assembled.cleaned.fastq ${NAME}.unassembled.forward.cleaned.fastq  ${NAME}.unassembled.reverse.cleaned.fastq > ${NAME}.fastq""")
   print >> oscript, """merged_fastq="${NAME}.fastq" """
   print >> oscript, """merged_fasta="${NAME}.fasta" """


def renameFastq(oscript, addPrefix):
   if not addPrefix:
      print >> oscript, """echo -e "\\n"; date; echo -e "START fastq_rename ${merged_fastq} > ${merged_fastq}_ && mv ${merged_fastq}_ ${merged_fastq}" """
      print >> oscript, """fastq_rename ${merged_fastq} > ${merged_fastq}_ && mv ${merged_fastq}_ ${merged_fastq}"""
      print >> oscript, """echo -e "\\n"; date; echo -e "END fastq_rename ${merged_fastq} > ${merged_fastq}_ && mv ${merged_fastq}_ ${merged_fastq}" """
   else:
      print >> oscript, """echo -e "\\n"; date; echo -e "START fastq_rename ${merged_fastq} ${NAME} > ${merged_fastq}_ && mv ${merged_fastq}_ ${merged_fastq}" """
      print >> oscript, """fastq_rename ${merged_fastq} ${NAME} > ${merged_fastq}_ && mv ${merged_fastq}_ ${merged_fastq}"""
      print >> oscript, """echo -e "\\n"; date; echo -e "END fastq_rename ${merged_fastq} ${NAME} > ${merged_fastq}_ && mv ${merged_fastq}_ ${merged_fastq}" """


def genericBlock(oscript, line):
   print >> oscript, """echo -e "\\n"; date; echo -e "START %s" """%(line.replace('"','\\"'))
   print >> oscript, line
   print >> oscript, """echo -e "\\n"; date; echo -e "END %s\\n" """%(line.replace('"','\\"'))


def buildMakeSamScript(childname, ident, threads):
   with open(childname, "w") as oscript:
      print >> oscript, "#!/bin/bash"
      print >> oscript, """: ${1?"Usage: $0 min_coverage max_coverage"}"""
      print >> oscript, """: ${2?"Usage: $0 min_coverage max_coverage"}"""
      print >> oscript, """minClustSize=${1}"""
      print >> oscript, """maxClustSize=${2}"""
      print >> oscript , """cd "%s" """%(ident)
      genericBlock(oscript, """make_sam_with_cons.py -u  %(short)s.mapping_to_cons -q  %(short)s.fastq -c  %(short)s_clean.ids -f  %(short)s.final.contigs.masked -l ${minClustSize} -m ${maxClustSize} multiple -o "%(short)s" """%dict(short = ident) )
      genericBlock(oscript, """jellyfish count -m 21 -s 512M -t %(threads)s -C %(ref)s -o %(short)s.jelly"""%dict(threads=threads, ref="%s_pseudo_ref_parts.fasta"%(ident), short=ident ) )
      genericBlock(oscript, """jellyfish dump  %(short)s.jelly > %(short)s.kmer.counts; rm %(short)s.jelly"""%dict(short=ident) )
      genericBlock(oscript, """filter_kmer_counts.py %(short)s.kmer.counts %(short)s 2; rm %(short)s.kmer.counts %(short)s_kmer.filter.tmp"""%dict(short=ident))
      print >> oscript, """cd .."""
      return ["%s.bam"%(ident), "%s.bam.bai"%(ident), "%s_pseudo_ref_parts.fasta"%(ident), "%s_pseudo_ref.fasta"%(ident), "%s_kmer.filter"%(ident) ]


def buildCommonStepOne(args, oscript, pairs, ident, threads, prefix = False):
   print >> oscript, "#!/bin/bash"
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
   print >> o, "#!/bin/bash"
   print >> o , "OLDDIR=`pwd`"
   print >> o, "THREADS=%s"%(threads) 
   print >> o, "DB_NAME='%s'\n"%(args.database) 
   if not basic:
      print >> o, "mkdir -p %s"%(args.workdir)	
      print >> o, "cd %s\n"%(args.workdir)	      
      print >> o, "mkdir csr"
      if clustered:
         print >> o, "mkdir csr/cluster_bams"


def buildChildRunner(args, children_scripts, passalongs, multi, threads = 1):
   o = open(CHILD_RUNNER%(1), "w")
   print >> o, "#!/bin/bash"
   if multi:
      oo = open(CHILD_RUNNER%(2), "w") 
      print >> oo, "#!/bin/bash"
      clusterSizeHdr(oo)
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
            childname = CHILD_SCRIPT%(d[0], c, 2)
            d.extend(buildMakeSamScript(childname, d[0], threads))
            if c % args.jobs == 0:
               if c != 0:
                  print >> o, "wait ${pidArr[@]}"
                  proc = []
               print >> o, """pidArr=()"""	
            minclust, maxclust = advanceNotice(oo, args, d[0])
            genericBlock(oo, """bash %s %s %s &\npidArr+=($!)"""%(childname, minclust, maxclust))
            proc.append(1)
         if proc:
            print >> o, "wait ${pidArr[@]}"
   else:	
      for cname in children_scripts:
         genericBlock(o, """bash %s"""%(cname))
      if multi:
         for c, d in enumerate(passalongs):
            childname = CHILD_SCRIPT%(d[0], c, 2)
            d.extend(buildMakeSamScript(childname, d[0], threads))
            minclust, maxclust = advanceNotice(oo, args, d[0])
            genericBlock(oo, """bash %s %s %s"""%(childname, minclust, maxclust))
   o.close()
   if multi:
      oo.close()
   return passalongs


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


def clusterSection(oscript):
   genericBlock(oscript, """vsearch -sortbylength ${NAME}.masked.fasta --output ${NAME}.masked.sorted.fasta --threads ${THREADS}""")
   genericBlock(oscript, """vsearch --cluster_smallmem ${NAME}.masked.sorted.fasta --strand plus --id 0.95 --msaout ${NAME}_1.msa --userout ${NAME}_1.out --userfields query+target+caln+qstrand+tstrand --mincols 80 --maxdiffs 10 --threads ${THREADS}""")
   genericBlock(oscript, """muso.py  -t ${THREADS} -i1 ${NAME}_1.msa -o1 ${NAME}_mod_1.cons -i2 ${NAME}_1.out -o2 ${NAME}_mod_1.out -g 1 -m 3 -n 2000""")
   genericBlock(oscript, """vsearch -sortbylength ${NAME}_mod_1.cons --output ${NAME}_mod_1.sorted.fasta -threads ${THREADS}""")
   genericBlock(oscript, """vsearch --cluster_smallmem ${NAME}_mod_1.sorted.fasta --strand both --id 0.95 --msaout ${NAME}_2.msa --userout ${NAME}_2.out --userfields query+target+caln+qstrand+tstrand --mincols 80 --maxdiffs 10 --threads ${THREADS}""")
   genericBlock(oscript, """muso.py  -t ${THREADS} -i1 ${NAME}_2.msa -o1 ${NAME}.final.contigs -i2 ${NAME}_2.out -o2 ${NAME}.final.out -g 2 -m 3 -n 2000""")
   genericBlock(oscript, """trackOverlaps.py -i1 ${NAME}_mod_1.out  -i2 ${NAME}.final.out  -o ${NAME}.mapping_to_cons""")


def maskingSection(oscript):
   genericBlock(oscript, """vsearch -maskfasta ${NAME}.final.contigs -hardmask  -output ${NAME}.final.contigs.masked -qmask dust -threads ${THREADS}""")
   genericBlock(oscript, """sed -i -E '/>/!s/[acgt]/N/g' ${NAME}.final.contigs.masked""")
   genericBlock(oscript, """grep -v '>'  ${NAME}.final.contigs.masked |  sed -e "1i>${NAME}" | seqret stdin ${NAME}.temp.pseudo.fasta""")
   genericBlock(oscript, """build_lmer_table -l 12 -sequence ${NAME}.temp.pseudo.fasta -freq freqs_${NAME}""")
   genericBlock(oscript, """RepeatScout -sequence ${NAME}.temp.pseudo.fasta -output repeats_${NAME}.fa -freq freqs_${NAME} -l 12""")
   genericBlock(oscript, """usearch -usearch_local repeats_${NAME}.fa  -db ${NAME}.final.contigs.masked -id 0.80 -strand both -userfields query+target+id+qcov+tcov+qlo+qhi+tlo+thi -userout ${NAME}_bad.ids -query_cov 0.75  -threads ${THREADS} --maxaccepts 0 --maxrejects 0""")
   genericBlock(oscript, """maskSeqs.py -i ${NAME}_bad.ids  -f ${NAME}.final.contigs.masked -o ${NAME}.final.contigs.masked_ && mv ${NAME}.final.contigs.masked_ ${NAME}.final.contigs.masked""")
   genericBlock(oscript, """grep ">" ${NAME}.final.contigs.masked | sed 's/>//' |  sed 's/<unknown description>//' > ${NAME}_clean.ids""")         


def singleMergedInput(args, threads, passalongs):
   comboname = "_".join(passalongs)
   with open(CHILD_RUNNER%(2), "w") as o:
      commonMainScript(o, args, threads, False)
      print >> o, "cd csr"
      print >> o, """NAME="%s" """%(comboname)
      genericBlock(o, """cat %s > ${NAME}.fastq """%( " " .join( [ """ "../%(ident)s/%(ident)s.fastq" """%dict(ident=d) for d in passalongs] ) ) )
      genericBlock(o, """seqtk seq -A ${NAME}.fastq > ${NAME}.fasta""")
      genericBlock(o, """vsearch --derep_fulllength  ${NAME}.fasta -minseqlength 1 -output ${NAME}.dedup.fasta -uc ${NAME}.uc --threads ${THREADS}""")
      genericBlock(o, """vsearch -maskfasta ${NAME}.dedup.fasta --hardmask --output ${NAME}.masked.fasta --threads ${THREADS}""")
      clusterSection(o)
      genericBlock(o, """update_mapping.py -i ${NAME}.uc -m ${NAME}.mapping_to_cons""")
      maskingSection(o)
      genericBlock(o, """coverageInformation.py -s "${NAME}"  -c ${NAME}_clean.ids -m ${NAME}.mapping_to_cons""")
   return comboname


def generateMulti(args):
   threads = max(1, int(math.floor( float(args.threads) / float(args.jobs))))
   parameters = yaml.load(open(args.config))
   children_scripts = []
   cid = 0
   passalongs = []		
   for ident, pairs in parameters['libraries'].iteritems():
      ident = ident.replace("_","-").replace(" ","-")
      scriptname = os.path.join(args.workdir, CHILD_SCRIPT%(ident, cid, 1))
      cid += 1
      children_scripts.append(scriptname)
      with open(scriptname, "w") as oscript:
         buildCommonStepOne(args, oscript, pairs, ident, threads, False)
         passalongs.append([ident])
         genericBlock(oscript, """seqtk seq -A ${merged_fastq} > ${merged_fasta}""")
         genericBlock(oscript, """vsearch -maskfasta ${merged_fasta} --hardmask --output ${NAME}.masked.fasta --threads ${THREADS}""")
         clusterSection(oscript)
         maskingSection(oscript)
         genericBlock(oscript, """coverageInformation.py -s "${NAME}"  -c ${NAME}_clean.ids -m ${NAME}.mapping_to_cons""")
         print >> oscript, """cd .."""

   buildChildRunner(args, children_scripts, passalongs, True, threads)
   threads = args.threads   
   with open(MAIN_SCRIPT, "w") as o:
      commonMainScript(o, args, threads, True)  

      print >> o, "bash %s"%(CHILD_RUNNER%(1))
      print >> o, "bash %s"%(CHILD_RUNNER%(2))
      print >> o, "bash %s"%(CHILD_RUNNER%(3))
      with open(CHILD_RUNNER%(3), "w") as oo:
         commonMainScript(oo, args, threads, True)
         minlenAndminSim(oo, parameters)
         print >> oo, "cd csr"
         for d in passalongs:
            genericBlock(oo, """mv ../%(parent)s/%(bam)s ../%(parent)s/%(bamidx)s cluster_bams/"""%dict(parent = d[0], bam = d[1], bamidx = d[2]) )
            
         print >> oo, """mkdir kmer_counts"""
         for d in passalongs:
            genericBlock(oo, """mv ../%(parent)s/%(cnt)s kmer_counts/"""%dict(parent = d[0], cnt = d[5]) )         
         # TODO: Need to insert an ordering function gives us the order in which seeds and and other things are processed...
         if len(passalongs) >= 2:
            genericBlock(oo, """findBestOrder_quick.py kmer_counts _kmer.filter search.order""")
            print >> oo, "declare -A mapping;"
            for d in passalongs:
               print >> oo, """mapping["%s"]="%s" """%(d[0], "../%s/%s"%(d[0], d[3]))
               
            print >> oo, """ordering=()"""
            print >> oo, """i=0"""
            print >> oo, """while read line"""
            print >> oo, """do"""
            print >> oo, """ ordering[$i]=${line}"""
            print >> oo, """ i=$(($i+1))"""
            print >> oo, """done < search.order"""
            genericBlock(oo, """Seanome.py -t ${THREADS} -d ${DB_NAME} seed_csr -i1 ${mapping[${ordering[0]}]} -n1 ${ordering[0]} -i2 ${mapping[${ordering[1]}]} -n2 ${ordering[1]}  -l ${minlen} -s ${minsim} """)
            if len(passalongs) > 2:
               print >> oo, """for i in {2..%s}"""%((len(passalongs)-1))
               print >> oo, """do"""
               genericBlock(oo, """ Seanome.py -t ${THREADS} -d ${DB_NAME} find_csr -g ${mapping[${ordering[$i]}]} -l ${minlen} -s ${minsim} """%dict(parent = d[0], ref = d[4]))
               print >>oo, """done"""
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
      print >> o, "cd ${OLDDIR}"
		

def generateSingle(args):
   threads = max(1, int(math.floor( float(args.threads) / float(args.jobs))))
   parameters = yaml.load(open(args.config))
   children_scripts = []
   cid = 0
   passalongs = []		
   for ident, pairs in parameters['libraries'].iteritems():
      ident = ident.replace("_","-").replace(" ","-")	
      scriptname = os.path.join(args.workdir, CHILD_SCRIPT%(ident, cid, 1))
      cid += 1
      children_scripts.append(scriptname)
      with open(scriptname, "w") as oscript:
         buildCommonStepOne(args, oscript, pairs, ident, threads, True)
	 print >> oscript, """cd .."""
      passalongs.append( ident  )

   buildChildRunner(args, children_scripts, passalongs, False, threads)
   threads = args.threads
   with open(MAIN_SCRIPT, "w") as o:
      commonMainScript(o, args, threads, False)
      print >> o, "bash %s"%(CHILD_RUNNER%(1))
      print >> o, "bash %s"%(CHILD_RUNNER%(2))
      print >> o, "bash %s"%(CHILD_RUNNER%(3))
      comboname = singleMergedInput(args, threads, passalongs)
      with open(CHILD_RUNNER%(3), "w") as oo:
         commonMainScript(oo, args, threads, False)
         print >> oo, "cd csr"
         print >> oo, """NAME="%s" """%(comboname)
         clusterSizeHdr(oo)
         minclust, maxclust = advanceNotice(oo, args, "${NAME}")   
         minlenAndminSim(oo, parameters, True)    
         genericBlock(oo, """make_sam_with_cons.py -u  ${NAME}.mapping_to_cons -q ${NAME}.fastq -c ${NAME}_clean.ids -f ${NAME}.final.contigs.masked -l %s -m %s single -d ${DB_NAME} -t ${THREADS} -s ${minlen} """%(minclust, maxclust))
         genericBlock(oo, """vcf_generator.py -t ${THREADS} -d  ${DB_NAME}""")
         genericBlock(oo, """vcfmod.py -t ${THREADS} -d  ${DB_NAME}""")
      print >> o, "cd ${OLDDIR}"


if __name__ == "__main__":
	parser = argparse.ArgumentParser()

	parser.add_argument('-c', '--config', required = True, help = "configuration file in yaml format" )
	parser.add_argument('-d', '--database', required = True, help = "Name of the seanome sqlite database" )
	parser.add_argument('-t', '--threads',  required = False, default = 1, type = int, help = "total number of threads to utilize (default: 1)")
	parser.add_argument('-j', '--jobs',  required = False, default = 1, type = int, help = "Number of jobs to run in parallel.  This will divide the number of threads, and undersubscribe in case of uneven division (default: 1)")
	#parser.add_argument('-w', '--workdir', required = False, default =".", help = "Output directory to run job in (deafult: current directory)")
        parser.add_argument("-a", "--advance", action = "store_true", required = False, help = "generates coverage information and requires the use to provide input")
        parser.set_defaults(workdir=".")

	subparsers = parser.add_subparsers(dest='action', help='Available commands')
	parser_sub = subparsers.add_parser('multiple')
	parser_sub.set_defaults(func = generateMulti)
	parser_sub = subparsers.add_parser('single')
	parser_sub.set_defaults(func = generateSingle)

	args = parser.parse_args()
	args.database = os.path.basename(args.database)
	args.func(args)
