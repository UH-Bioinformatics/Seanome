#!/usr/bin/env python2.7

import argparse
import sys
import yaml
import json
from collections import defaultdict


def generateCoverage(cleanids, mapping):
    valid = set([i.strip() for i in open(cleanids)])
    coverage = defaultdict(int)
    for l in open(mapping):        
        l = l.strip().split()
        if l[1] in valid:
            coverage[l[1]] += 1
    return coverage


def generateCounterPerCoverage(coverage):
    buckets = defaultdict(int)
    for c in coverage.itervalues():
        buckets[c] += 1
    return buckets


if '__main__' == __name__:
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cleanids', required = True, help = "clean id listing file" )
    parser.add_argument('-s', '--shortname', required = True, help = "short name representing this library" )
    parser.add_argument('-m', '--mapping', required = True, help = "mapping file from trackoverlap" )

    args = parser.parse_args()

    coverage = generateCoverage(args.cleanids, args.mapping)
    with open("%s.coverage.yaml"%(args.shortname), "w") as o:
        yaml.dump(dict(coverage), o)
    
    buckets = generateCounterPerCoverage(coverage)
    with open("%s.coverage.json"%(args.shortname), "w") as o:
        o.write(json.dumps([args.shortname, sorted(buckets.items(), key = lambda x: x[0]) ], separators=(',',':')))
