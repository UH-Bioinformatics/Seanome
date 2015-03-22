#!/usr/bin/python

import os
import sys
import argparse
import logging
import stat

from utils.inferSAM import inferSAM
from utils.sqlitedb import buildsqlitedb
from utils.threadpool import ProducerConsumer
from utils.findCSR import find_csr
from utils.findSeedCSR import find_seed_csr
from utils.makeCons import generateConsensusSequences
from utils.trimall import runTrimAL
from utils.cleanSAM import cleanSAM
from utils.vcf_generator import generateVCFfiles
from utils.modvcf import updateVCFfiles

#https://stackoverflow.com/questions/1889597/deleting-directory-in-python/1889686#1889686
# to handle the case of read only files in the remove tree for shutil.rmtree
def remove_readonly(func, path, excinfo):
    os.chmod(path, stat.S_IWRITE)
    func(path)

def runInstance(myInstance):
   dryRunInstance(myInstance)
   myInstance.run()


def dryRunInstance(myInstance):
   logging.warning(myInstance.dryRun())


def makeDirOrdie(dirPath):
   if not os.path.isdir(dirPath):
      os.makedirs(dirPath)
   else:
       # TODO; Uncomment after testing done
      #logging.error("Split fasta directory %s already exists " % hmmerOutputDir)
      #sys.exit()
      pass
   return dirPath

#    logging.debug("Finding seed common shared regions in %s: \n Starting... "%(args.input))



def main(argv):
    parser = argparse.ArgumentParser(description="Seanome description", epilog="Seanome long text description")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s '+version)
    parser.add_argument('-t', '--threads', type=int, default = 1)
    parser.add_argument("-q", "--quiet", action = "store_true", required = False, help = "Silence subprocess output")
    parser.add_argument('-d', '--database', required=True, help="Seanome sqlite database")
    subparsers = parser.add_subparsers(dest='action', help='Available commands')

    # seed_csr
    parser_mask = subparsers.add_parser('seed_csr')
    parser_mask.add_argument('-i1', '--input1',  required=True, help="pseudo references parts 1")
    parser_mask.add_argument('-n1', '--name1',  required=True, help="short name for references parts 1")   
    parser_mask.add_argument('-i2', '--input2',  required=True, help="pseudo references parts 2")
    parser_mask.add_argument('-n2', '--name2',  required=True, help="short name for references parts 2")
    parser_mask.add_argument('-l', '--min_csr_len', type=int, default=150, help=" Minimum common shared region length")
    parser_mask.add_argument('-s', '--min_csr_sim', type=float, default=0.88, help=" Minimum common shared region similarity")
    parser_mask.set_defaults(func = find_seed_csr)

    # csr
    parser_mask = subparsers.add_parser('find_csr')
    parser_mask.add_argument('-g', '--genome', required=True, help="Reference genome") # pseudo genome 
    parser_mask.add_argument('-l', '--min_csr_len', type=int, default=150, help=" Minimum common shared region length")
    parser_mask.add_argument('-s', '--min_csr_sim', type=float, default=0.88, help=" Minimum common shared region similarity")
    parser_mask.set_defaults(func = find_csr)

    # Build transitive SAM from MSA
    parser_mask = subparsers.add_parser('inferSAM')
    parser_mask.add_argument('-s', '--split_sam_dir', required = True, help = 'Split SAM Files directory')
    parser_mask.set_defaults(func=inferSAM)
    
    # consensus
    parser_mask = subparsers.add_parser('consensus')
    parser_mask.set_defaults(func=generateConsensusSequences)

    # trimAL
    parser_mask = subparsers.add_parser('trimAL')
    parser_mask.set_defaults(func=runTrimAL)

    # cleanSAM
    parser_mask = subparsers.add_parser('cleanSAM')
    parser_mask.set_defaults(func=cleanSAM)

    # cleanSAM
    parser_mask = subparsers.add_parser('generateVCF')
    parser_mask.set_defaults(func=generateVCFfiles)

    # cleanSAM
    parser_mask = subparsers.add_parser('updateVCF')
    parser_mask.set_defaults(func=updateVCFfiles)

    # Parse arguments
    args = parser.parse_args()
    #logging.debug("Initial ARGS are:")
    #logging.debug(args)

    con = buildsqlitedb(args.database)
    args.func(args, con)
    con.commit()
    con.close()
    


if __name__ == "__main__":
   # test
   version = "alpha 0.01"
   FORMAT = "%(asctime)-15s  %(message)s"
   logging.basicConfig(format=FORMAT, level=logging.DEBUG)
   main(sys.argv)
