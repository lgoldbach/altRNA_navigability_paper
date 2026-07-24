#!/usr/bin/env python

import argparse
from rna_folding.parsing import dict_to_gpmap
from rna_folding.gp_map import GenotypePhenotypeGraph
import pickle

import numpy as np


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", help="Output file", type=str)
    parser.add_argument("-i", "--input", help="input gp_map", type=str)
    parser.add_argument("-u", "--unfolded", help="unfolded phenotype", type=str)
    parser.add_argument("-a", "--alphabet", help="The sequence alphabet, e.g. 'AUGC'", type=str, required=True)
    parser.add_argument("-t", "--genotypes", help="Genotype list", type=str, required=True)
    parser.add_argument("-r", "--ranking", help="A file containing phenotypes that defines the ranking. IMPORTANT: The last phenotype has to be the unfolded", type=str)

    
    args = parser.parse_args()
    with open(args.ranking, "r") as f:
        phenotypes = [line.strip().split()[0] for line in f if line.strip()]

    phenotypes.append(args.unfolded)  # add unfolded to the end

    ph_to_rank = dict([(ph, i) for i, ph in enumerate(phenotypes)])

    r_unf = len(phenotypes)-1  # last rank to unfolded phenotype

    flat_gp_map = {}
    with open(args.input, "r") as gp_map:
        for line in gp_map:
            line_ = line.strip().split(" ")
            ph = line_[0]

            try:
                rank = ph_to_rank[ph]
            except KeyError:  # the phenotype is not part of the ranking
                rank = r_unf
                ph_to_rank[ph] = rank  # assign it the lowest rank (unfolded ph)
                
            for gt in line_[1:]:
                if gt not in flat_gp_map:  # O(1) operation because dict
                    flat_gp_map[gt] = rank  # assign phenotype
                elif rank < flat_gp_map[gt]:
                    flat_gp_map[gt] = rank  # update to better ranked ph
        

    with open(args.genotypes, "r") as gt_file:
        genotypes = [line.strip() for line in gt_file]

    # # turn rank into phenotype and genotype number into sequence
    flat_gp_map_str = {}
    for gt in flat_gp_map:
        ph = phenotypes[flat_gp_map[gt]]  # translate rank to ph
        seq = genotypes[int(gt)]  # genotype number to str
        flat_gp_map_str[seq] = ph  # assign phenotype instead of rank now
    
    # Turn into gp graph object (required to get neutral components)
    print(args.alphabet, flush=True)
    gp_map = GenotypePhenotypeGraph.read_from_dict(flat_gp_map_str, alphabet=args.alphabet)
    
    # turn into ph_to_gt map
    ph_to_gt = dict([(ph, []) for ph in phenotypes])
    for i, gt in enumerate(genotypes):
        ph_to_gt[flat_gp_map_str[gt]].append(i)  # append number genotype to the ph it maps to
    
    # save to file
    dict_to_gpmap(ph_to_gt, args.output)

