#!/usr/bin/env python

import sys
import os
import re
from Bio.PDB import PDBParser


def calculate_BDT(pdb_file, bindres_file, pred_file, d_zero):
    pdbid_chain = os.path.basename(pred_file)[:5]
    if re.match(r'[a-z]', pdbid_chain[-1]):
        pdbid_chain = pdbid_chain[:4] + pdbid_chain[4].upper()
    
    p = PDBParser()
    s = p.get_structure("test", pdb_file)
    
    cas = []
    chain = s[0][pdbid_chain[-1]]
    for res in chain.get_residues():
        if res.id[0] == ' ':
            cas.append(res['CA'])
    
    bind_record = []
    with open(bindres_file) as fp:
        for line in fp:
            if re.match(r'[01]', line):
                bind_record = map(int, list(line.rstrip()))
    binding_indexes = set([i for i, record in enumerate(bind_record) if record == 1])
    
    preds = []
    with open(pred_file) as fp:
        for line in fp:
            pred = line.rstrip().split('\t')[2]
            if pred == "1":
                preds.append(1)
            else:
                preds.append(0)
    pred_indexes = [i for i, pred in enumerate(preds) if pred == 1]
    
    
    N = max(len(binding_indexes), len(pred_indexes))
    Ss = []
    for i in pred_indexes:
        if i in binding_indexes:
            Ss.append(1.0)
        else:
            min_distance_index = (float("inf"), -1)
            for j in binding_indexes:
                distance = cas[i] - cas[j]
                if distance < min_distance_index[0]:
                    min_distance_index = (distance, j)
            Ss.append(1 / (1 + (min_distance_index[0] / d_zero)**2))
    BDT = sum(Ss) / N
    return BDT


if __name__ == "__main__":
    if len(sys.argv) <= 1 or sys.argv[1] in {"-h", "--help"}:
        print "Usage: python calculate_bdt.py <pdb_file> <bindres_file> <pred_file> <d_zero>"
        sys.exit()

    pdb_file = sys.argv[1]
    bindres_file = sys.argv[2]
    pred_file = sys.argv[3]
    d_zero = float(sys.argv[4])
    print calculate_BDT(pdb_file, bindres_file, pred_file, d_zero)
