#!/usr/bin/python

import argparse
import sqlite3

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extracts the trimmed sam files')
    parser.add_argument('-d', '--database', help='Specified Seanome Database', required = True)

    args = parser.parse_args()
    con = sqlite3.connect(args.database, check_same_thread=False)
    cur = con.cursor()
    cur.execute("""SELECT 
                         B.name, A.IDs, A.SEQS 
                   FROM 
                        (SELECT fileID, group_concat(seqID, '\t') AS IDs, group_concat(sequence, '\t') AS SEQS FROM trimmed_csr group by fileID) AS A 
                   JOIN 
                        files AS B 
                   ON 
                        (A.fileID = B.id); """)
    for dat in cur:
        if not dat:
            continue
        name = dat[0]
        seqs = "\n".join([">%s\n%s"%(str(i), str(s)) for i, s in zip(dat[1].split("\t"), dat[2].split("\t"))])
        with open("%s.trimmed.csr.fasta"%(name), "w") as o:
            o.write(seqs)
    con.commit()
    con.close()
