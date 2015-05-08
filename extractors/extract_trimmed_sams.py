#!/usr/bin/env python2.7

import argparse
import sqlite3

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extracts the trimmed sam files')
    parser.add_argument('-d', '--database', help='Specified Seanome Database', required = True)

    args = parser.parse_args()
    con = sqlite3.connect(args.database, check_same_thread=False)
    cur = con.cursor()
    cur.execute(""" SELECT B.name, A.sam FROM trimmed_inferSAM AS A JOIN files AS B ON (A.fileID = B.id); """)
    for dat in cur:
        if not dat:
            continue
        name = dat[0]
        sam = dat[1]
        with open("%s.trimmed.sam"%(name), "w") as o:
            o.write(sam)
    con.commit()
    con.close()
