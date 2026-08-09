#!/usr/bin/env python

import argparse
import pickle
import networkx as nx
import numpy as np
import json

from rna_gpf.parsing import gpmap_to_dict

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Input genotype-phenotype map "
                        "file", required=True)
    parser.add_argument("-g", "--genotypes", help="genotype file", required=True)
    parser.add_argument("-r", "--ref", help="Reference gp map", type=str, required=False)
    parser.add_argument("-u", "--unfolded", help="Phenotype to ignore, e.g "
                        "unfolded", type=str, required=False)
    parser.add_argument("-o", "--output", help="Output file", required=True)

    args = parser.parse_args()
    
    gpm = pickle.load(open(args.file, "rb"))
    with open(args.genotypes, "r") as gt_file:
        genotypes = [line.strip() for line in gt_file]

    ref_gpm = pickle.load(open(args.ref, "rb"))

    outdir = {}

    alphabet_trans = {"J": "G", "K": "A", "L": "C", "M": "U"}

    map_gts = lambda gt: "".join([alphabet_trans[l] for l in gt])

    for gt in genotypes:
        try:
            ph = gpm.map(gt)
        except KeyError:
            ph = args.unfolded

        gt_ref = map_gts(gt)

        try:
            ref_ph = ref_gpm.map(gt_ref)
        except KeyError:
            ref_ph = args.unfolded

            if ph == args.unfolded:
                outdir["both_unf"] = outdir.get("both_unf", 0) + 1
            else:
                outdir["ref_unf"] = outdir.get("ref_unf", 0) + 1
            continue

        if ref_ph == ph:
            outdir["folded_match"] = outdir.get("folded_match", 0) + 1
        elif ph == args.unfolded:
            outdir["que_unf"] = outdir.get("que_unf", 0) + 1
        else:
            outdir["folded_mismatch"] = outdir.get("folded_mismatch", 0) + 1

    with open(args.output, "w") as out:
        json.dump(outdir, out, indent=4)
