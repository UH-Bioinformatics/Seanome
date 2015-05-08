#!/usr/bin/env python2.7

import argparse
import sqlite3

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extracts the vcff and dat files after vcfmod.py')
    parser.add_argument('-d', '--database', help='Specified Seanome Database', required = True)

    args = parser.parse_args()
    con = sqlite3.connect(args.database, check_same_thread=False)
    cur = con.cursor()
    cur.execute(""" SELECT B.name, A.vcf, A.json FROM trimmed_modvcf AS A JOIN files AS B ON (A.fileID = B.id); """)
    for dat in cur:
        if not dat:
            continue
        name = dat[0]
        vcf = dat[1]
        json = dat[2]
        with open("%s.mod.vcf"%(name), "w") as o:
            o.write(vcf)
        with open("%s.mod.dat"%(name), "w") as o:
            o.write(json)            
    con.commit()
    con.close()
