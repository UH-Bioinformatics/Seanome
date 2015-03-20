#!/usr/bin/python

import argparse
import sqlite3

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extracts the vcf files after vcf_generator.py')
    parser.add_argument('-d', '--database', help='Specified Seanome Database', required = True)

    args = parser.parse_args()
    con = sqlite3.connect(args.database, check_same_thread=False)
    cur = con.cursor()
    cur.execute(""" SELECT B.name, A.vcf FROM trimmed_vcf AS A JOIN files AS B ON (A.fileID = B.id); """)
    for dat in cur:
        if not dat:
            continue
        name = dat[0]
        vcf = dat[1]
        with open("%s.vcf"%(name), "w") as o:
            o.write(vcf)
    con.commit()
    con.close()
