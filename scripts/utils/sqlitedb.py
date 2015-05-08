import sqlite3


QUERY_CSR_AS_SEQS = """SELECT fileID, group_concat(seqID, '\t') AS IDs, group_concat(sequence, '\t') AS SEQS FROM csr GROUP BY fileID"""
QUERY_TRIMMED_CSR_AS_SEQS = """SELECT fileID, group_concat(seqID, '\t') AS IDs, group_concat(sequence, '\t') AS SEQS FROM trimmed_csr GROUP BY fileID"""


def buildsqlitedb(dbname, ismulti = True):
    if ismulti:
        return buildMultieDb(dbname)
    else:
        return buildSingleDB(dbname)


def commontables(dbname):
    con = sqlite3.connect(dbname, check_same_thread=False)
    con.execute("""PRAGMA foreign_keys = ON;""")

    con.execute("""CREATE TABLE IF NOT EXISTS files(id INTEGER PRIMARY KEY, name TEXT);""")
    con.execute("""CREATE INDEX IF NOT EXISTS file_name_idx ON files(name ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS groups(id INTEGER PRIMARY KEY, fileID INTEGER, groupid TEXT,  sampleid Text,
                   species TEXT, coverage INTEGER, trimmed_coverage INTEGER DEFAULT 0, FOREIGN KEY(fileID) REFERENCES files(id) );""")
    con.execute("""CREATE INDEX IF NOT EXISTS groups_fileid_idx ON groups(fileID ASC);""")
    con.execute("""CREATE INDEX IF NOT EXISTS groups_species_idx ON groups(species ASC);""")
    con.execute("""CREATE INDEX IF NOT EXISTS groups_group_idx ON groups(groupid ASC);""")
    con.execute("""CREATE INDEX IF NOT EXISTS groups_sampleid_idx ON groups(sampleid ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_consensus(id INTEGER PRIMARY KEY, fileID INTEGER, sequence TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_con_fileid_idx ON trimmed_consensus(fileID ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_inferSAM(id INTEGER PRIMARY KEY, fileID INTEGER, sam TEXT, bam BLOB, bamidx BLOB, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_infersam_fileid_idx ON trimmed_inferSAM(fileID ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_vcf(id INTEGER PRIMARY KEY, fileID INTEGER, vcf TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_vcf_fileid_idx ON trimmed_vcf(fileID ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_modvcf(id INTEGER PRIMARY KEY, fileID INTEGER, vcf TEXT, json TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_modvcf_fileid_idx ON trimmed_modvcf(fileID ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS fst(id INTEGER PRIMARY KEY, fileID INTEGER, groupA INTEGER, groupB INTEGER, pos INTEGER, value REAL, FOREIGN KEY(fileID) REFERENCES files(id), FOREIGN KEY(groupA) REFERENCES groups(id), FOREIGN KEY(groupB) REFERENCES groups(id) );""")
    con.execute("""CREATE INDEX IF NOT EXISTS fst_name_idx ON fst(fileID ASC);""")
    con.execute("""CREATE INDEX IF NOT EXISTS fst_grpA_idx ON fst(groupA ASC);""")
    con.execute("""CREATE INDEX IF NOT EXISTS fst_grpB_idx ON fst(groupB ASC);""")

    return con


def buildSingleDB(dbname):
    con = commontables(dbname)
    con.commit()
    return con


def buildMultieDb(dbname):
    con = commontables(dbname)

    con.execute("""CREATE TABLE IF NOT EXISTS csr(id INTEGER PRIMARY KEY, fileID INTEGER, seqID TEXT, 
                   sequence TEXT, seed BOOLEAN DEFAULT 0, FOREIGN KEY(fileID) REFERENCES files(id));""")

    con.execute("""CREATE TABLE IF NOT EXISTS consensus(id INTEGER PRIMARY KEY, fileID INTEGER, sequence TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS con_fileid_idx ON consensus(fileID ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_csr(id INTEGER PRIMARY KEY, fileID INTEGER, seqID TEXT, sequence TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_csr_fileid_idx ON trimmed_csr(fileID ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_logs(id INTEGER PRIMARY KEY, fileID INTEGER, positions TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_log_fileid_idx ON trimmed_logs(fileID ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS inferSAM(id INTEGER PRIMARY KEY, fileID INTEGER, sam TEXT, bam BLOB, bamidx BLOB, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS infersam_fileid_idx ON inferSAM(fileID ASC);""")

    con.commit()
    return con
