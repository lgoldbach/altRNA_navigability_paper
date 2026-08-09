#!/usr/bin/env python

import argparse
import pickle
import networkx as nx
import numpy as np

from rna_gpf.parsing import gpmap_to_dict

def parse_pairs(dotbracket):
    """Return the set of base pairs (i, j), i < j, in a pseudoknot-free
    dot-bracket string, via a simple stack."""
    stack = []
    pairs = set()
    for i, c in enumerate(dotbracket):
        if c == "(":
            stack.append(i)
        elif c == ")":
            j = stack.pop()
            pairs.add((j, i))
    return pairs


def f1_score(reference, query):
    """
    reference: phenotype A, dot-bracket string
    query: phenotype B, dot-bracket string
    Returns the base-pair F1 score (float in [0, 1]).

    Edge cases: if both structures have zero base pairs (both fully
    unpaired), they're identical and F1 = 1.0. If exactly one has zero
    base pairs (so there's no possible overlap), F1 = 0.0.
    """
    pairs_a = parse_pairs(reference)
    pairs_b = parse_pairs(query)

    tp = len(pairs_a & pairs_b)
    fn = len(pairs_a - pairs_b)
    fp = len(pairs_b - pairs_a)

    if tp == 0:
        return 1.0 if (not pairs_a and not pairs_b) else 0.0

    sensitivity = tp / (tp + fn)
    ppv = tp / (tp + fp)
    return 2 * sensitivity * ppv / (sensitivity + ppv)


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Input genotype-phenotype map "
                        "file", required=True)
    parser.add_argument("-r", "--ref", help="Reference gp map", type=str, required=False)
    parser.add_argument("-u", "--unfolded", help="Phenotype to ignore, e.g "
                        "unfolded", type=str, required=False)
    parser.add_argument("-o", "--output", help="Output file", required=True)

    args = parser.parse_args()
    
    gpm = pickle.load(open(args.file, "rb"))
    ref_gpm = pickle.load(open(args.ref, "rb"))

    ph_match_counts = {ph: 0 for ph in gpm.phenotype_set}
    ph_unfolded_counts = {ph: 0 for ph in gpm.phenotype_set}
    ph_diff_counts = {ph: 0 for ph in gpm.phenotype_set}

    alphabet_trans = {"J": "G", "K": "A", "L": "C", "M": "U"}

    map_gts = lambda gt: "".join([alphabet_trans[l] for l in gt])

    f1_scores = []
    ref_phs = []
    query_phs = []

    folded_only_f1s = []
    for gt in gpm.genotypes:
        ph = gpm.map(gt)

        gt_ref = map_gts(gt)

        try:
            ref_ph = ref_gpm.map(gt_ref)
        except KeyError:
            ref_ph = args.unfolded
        
        f1 = f1_score(ref_ph, ph)
        if ref_ph != args.unfolded and ph != args.unfolded:
            folded_only_f1s.append(f1)

        f1_scores.append(f1)
        ref_phs.append(ref_ph)
        query_phs.append(ph)

    with open(args.output, "w") as out:
        out.write(f"{np.mean(f1_scores)}\n\n")
        out.write(f"{np.mean(folded_only_f1s)}\n\n")
        for f1, ref_ph, qu_ph in zip(f1_scores, ref_phs, query_phs):
            out.write(f"{ref_ph} {qu_ph} {f1}\n")
    