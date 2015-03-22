import os

CONSENSUS_NAME = "Consensus"


def removeFiles(flist):
    for f in flist:
        try:
            os.remove(f)
        except:
            pass
