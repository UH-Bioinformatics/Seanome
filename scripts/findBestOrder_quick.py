#!/usr/bin/python

# TODO: Clean the directories in alignment that are not optimal
# capture control+c and stop the execution of the threads?

# when running via django, we run into an error since MPLCONFIGDIR points into the root home..
# give it some temporary location to spawn user generated files..
import os

#os.environ['MPLCONFIGDIR'] = "/tmp/"

#from scipy import cluster
import numpy as np
import os
import sys
import math



def computeBrayCurtis(parent_dict, probe):
    distance = [0.0, 0.0]
    visited = set()

    for l in probe:
        kmer, cnt = l.strip().split("\t")
        p = float(parent_dict.get(kmer, 0.0))
        q = float(cnt)
        visited.add(kmer)
        distance[0] +=  abs(p - q) 
        distance[1] +=  abs(p + q) 
    for k in parent_dict.iterkeys():
        if k in visited:
            continue
        p = float(parent_dict[k])
        distance[0] +=  abs(p)
        distance[1] +=  abs(p)
        
    distance = distance[0] / distance[1]
    return distance
    


def computeBhattacharyya(parent_dict, probe):
    distance = 0.0
    for l in probe:
        kmer, cnt = l.strip().split("\t")
        if not kmer in parent_dict:
            continue
        
        #p = float(parent_dict[kmer]) / 100.0
        #q = float(cnt) / 100.0

        p = float(parent_dict[kmer])
        q = float(cnt)
        distance += math.sqrt( p * q)
    #distance = math.sqrt(1 - distance)
    #print distance
    try:
        distance = -1 * math.log(distance)
    except:
        distance = float("inf")
    return distance
    

def computeMinkowski(parent_dict, probe, power):
    visited = set()
    distance = 0.0
    for l in probe:
        kmer, cnt = l.strip().split("\t")
        visited.add(kmer)
        distance += abs(float(parent_dict.get(kmer, 0.0)) - float(cnt)) ** power
    for k in parent_dict.iterkeys():
        if k in visited:
            continue
        distance += (abs(float(parent_dict[k])) ** power)

    distance = distance ** (1.0/power)
    #math.sqrt(distance)
    return distance    


def computeNormalizedEuclidean(parent_dict, probe):
    elements =0
    visited = set()
    distance = 0.0
    for l in probe:
        kmer, cnt = l.strip().split("\t")
        visited.add(kmer)
        distance += (float(parent_dict.get(kmer, 0.0)) - float(cnt)) ** 2
        elements += 1
    for k in parent_dict.iterkeys():
        if k in visited:
            continue
        distance += (float(parent_dict[k]) ** 2)
        elements += 1
    distance = math.sqrt(distance) / math.sqrt(float(elements))
    return distance    


def computeEuclidean(parent_dict, probe):
    return computeMinkowski(parent_dict, probe, 2)


def computeManhattan(parent_dict, probe):
    return computeMinkowski(parent_dict, probe, 1)


def computeCanberra(parent_dict, probe):
    distance = 0.0
    visited = set()
    for l in probe:
        kmer, cnt = l.strip().split("\t")
        p = float(parent_dict.get(kmer, 0.0))
        q = float(cnt)
        visited.add(kmer)
        distance +=  (abs(p - q) / (abs(p) + abs(q)))
    for k in parent_dict.iterkeys():
        if k in visited:
            continue
        p = float(parent_dict[k])
        distance +=  (abs(p) / (abs(p))) 
    return distance
    

def returnBestPutativePerm(path, suffix, measure = "score"):
    ''' 
    Returns the most likely permutation 
    '''

    files = [ (f.replace(sys.argv[2], ""), os.path.join(sys.argv[1], f), ) for f in os.listdir(sys.argv[1]) if f.endswith(sys.argv[2])]
    files.sort(key = lambda x: x[0])
    namesList = [ f[0] for f in files]
    tempNamesHash=dict([ (f[0], i) for i, f in enumerate(files)])

    distMatrix =  np.matrix([[float("nan") for i in xrange(len(files))] for j in xrange(len(files))])    
    for i in xrange(len(files)):
        parent = files[i]
        parent_dict = dict((l.strip().split("\t") for l in open(parent[1])) )
        for j in xrange(i+1, len(files)):
            with open(files[j][1]) as probe:
                #distance = computeCanberra(parent_dict, probe)
                #distance = computeBrayCurtis(parent_dict, probe)

                #distance = computeEuclidean(parent_dict, probe)
                #distance = computeManhattan(parent_dict, probe)
                #distance = computeNormalizedEuclidean(parent_dict, probe)
                #distance = computeMinkowski(parent_dict, probe, 4)
                #print files[j][1]
                #print files[i][1]
                distance = computeBhattacharyya(parent_dict, probe)

            distMatrix[i,j]=distance 
            distMatrix[j,i]=distance

    order=[] 
    #print tempNamesHash
    #print distMatrix
    #print namesList
    minValue = distMatrix[np.isnan(distMatrix) == False].min()
    #print minValue
    order.append(namesList[np.where(distMatrix == minValue)[0][0,0]])
    order.append(namesList[np.where(distMatrix == minValue)[1][0,0]])
    tempNamesHash.pop(order[0])
    tempNamesHash.pop(order[1])

    while len(tempNamesHash) > 0 : 
        #if measure  == "score":
        #    measureSoFar = 0
        #else:
        measureSoFar =  float("inf")
        index = -1
        # order from smallest to largest distance measure
        for i in tempNamesHash.keys():
            avgDist = float(sum([distMatrix[namesList.index(x),tempNamesHash[i]] for x in order]))/float(len(order) )

            if avgDist < measureSoFar or np.isinf(avgDist):
                index = tempNamesHash[i]
                measureSoFar = avgDist
        order.append(namesList[index]) 
        tempNamesHash.pop(namesList[index])
    
    return order


if len(sys.argv) != 4:
    print "USAGE: %s <directory> <kmer counter suffix> <outfile>"%(sys.argv[0])
    sys.exit(1)
    

order = returnBestPutativePerm(sys.argv[1], sys.argv[2], "score")
with open(sys.argv[3], "w") as o:
    o.write("\n".join(order))
#print "best initial tree is: %s" % order



