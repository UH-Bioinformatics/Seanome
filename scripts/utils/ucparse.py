import sys
import os
from collections import defaultdict

class parser(object):
    def __init__(self, infile):
        self.infile = infile
        self.clusters = defaultdict(list)
        self.reps ={}


    def getClusters(self, val=None):
        if val:
            return self.clusters[val]
        else:
            return self.clusters


    def getReps(self):
        return self.reps


    def parse(self):
        pass


class uclustUserParser(parser):
    def __init__(self, infile):
        super(self.__class__, self).__init__(infile)


    def parse(self, index, use, cutoff = 0, maxcutoff = 0):
        if(not os.path.exists(self.infile)):
            print >> sys.stderr, "usearch did not generate an output.  Make sure usearch executed properly."
            return False
        maxv = 0
        with open(self.infile) as cfile:
            for l in cfile:
                l = l.strip().split()
                if l[1] not in use:
                    continue
                if l[5] == '*':
                    self.clusters[l[1]].append( (l[0], l[2], l[3], index[l[0]], False,) )
                else:
                    tmp = list(index[l[5]])
                    tmp[0] = l[0]
                    tmp = tuple(tmp)
                    self.clusters[l[1]].append( (l[0], l[2], l[3], tmp, False,) )
                maxv = max(maxv, len( self.clusters[l[1]]) )

        if maxv < cutoff:
            print >> sys.stderr, "Reducing the minimum coverage value to %s from %s, due to no cluster containing the original minimum coverage"%(maxv, cutoff)
            #print >> sys.stdout, "Reducing the minimum coverage value to %s from %s, due to no cluster containing the original minimum coverage"%(maxv, cutoff)
            cutoff = maxv
        for k in self.clusters.keys():
            if cutoff >= 0 and len(self.clusters[k]) < cutoff:
                self.clusters.pop(k, None)
            elif maxcutoff > 0 and len(self.clusters[k]) > maxcutoff:
                self.clusters.pop(k, None)
        return True


