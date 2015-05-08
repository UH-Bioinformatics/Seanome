#!/usr/bin/env python2.7

import argparse
import json
import copy
from collections import defaultdict

def parse(infile):
    mapping = defaultdict(list)
    with open(infile) as cfile:
        for line in cfile:
            if line.startswith("C"):
                continue
            if line.startswith("S"):
                continue
            line = line.strip().split()
            mapping[line[-1]].append(line[-2])
    return mapping


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Add the deduplicated sequences to the mapping file, along with who they were duplicates of")
    parser.add_argument('-i', '--input',  required = True, help = "vsearch/usearch uc file")
    parser.add_argument('-m', '--map',  required = True, help = "mapping_to_cons file")

    args = parser.parse_args()
 
    existing = {}
    with open(args.map) as infile:
        for line in infile:
            parts =line.strip().split()
            existing[parts[0]] = parts
    collection = parse(args.input)
    with open(args.map, "a") as o:
        for k, v in collection.iteritems():
            row = copy.copy(existing[k])
            for dat in v:
                row[0] = dat
                row[-1] = k
                print >> o, "\t".join(row)
