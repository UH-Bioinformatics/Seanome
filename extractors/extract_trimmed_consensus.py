#!/usr/bin/env python2.7

import argparse
import sqlite3

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extracts the trimmed sam files')
    parser.add_argument('-d', '--database', help='Specified Seanome Database', required = True)

    args = parser.parse_args()
    con = sqlite3.connect(args.database, check_same_thread=False)
    cur = con.cursor()
    cur.execute(""" SELECT B.name, A.sequence FROM trimmed_consensus AS A JOIN files AS B ON (A.fileID = B.id); """)
    for dat in cur:
        if not dat:
            continue
        name = dat[0]
        seq = dat[1]
        with open("%s.trimmed.con.fasta"%(name), "w") as o:
            print >> o, ">Consensus\n%s"%(seq)
    con.commit()
    con.close()
