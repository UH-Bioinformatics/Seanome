import os

def removeFiles(flist):
    for f in flist:
        try:
            os.remove(f)
        except:
            pass
