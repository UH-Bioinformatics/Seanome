#/bin/bash

: ${1?"Usage: $0 ids_lst bamsuffix reference"}
: ${2?"Usage: $0 ids_lst bamsuffix reference"}
: ${3?"Usage: $0 ids_lst bamsuffix reference"}

# we are assuming this jar file is on the user path
vscan=`which VarScan.jar`
samtools mpileup -f ${3}  -b ${1} | java -jar ${vscan} mpileup2snp --min-var-freq 0.01 --output-vcf 1   2> /dev/null

rm *_${2} 2> /dev/null
rm ${1} 2> /dev/null

