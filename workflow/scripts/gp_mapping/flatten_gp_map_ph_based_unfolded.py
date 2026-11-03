#!/usr/bin/env python

import argparse
from rna_gpf.parsing import dict_to_gpmap
from rna_gpf.utils import load_phenotype_and_metric_from_file
from rna_gpf.gp_map import GenotypePhenotypeGraph
import pickle

import numpy as np


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", help="Output file", type=str)
    parser.add_argument("-i", "--input", help="input gp_map", type=str)
    parser.add_argument("-u", "--unfolded", help="unfolded phenotype", type=str)
    parser.add_argument("-p", "--unf_probs", help="per phenotype unfolded prob", type=str)
    parser.add_argument("-c", "--nc_cutoff", help="Neutral components smaller than the size cutoff will be mapped to unfolded phenotype", type=int, required=True)
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

    per_ph_unf_prob = dict(zip(*load_phenotype_and_metric_from_file(args.unf_probs)))

    # # turn rank into phenotype and genotype number into sequence
    flat_gp_map_str = {}
    for gt in flat_gp_map:
        ph = phenotypes[flat_gp_map[gt]]  # translate rank to ph
        seq = genotypes[int(gt)]  # genotype number to str

        if np.random.sample() < per_ph_unf_prob[ph]:
            flat_gp_map_str[seq] = args.unfolded  # assign unfolded
        else:
            flat_gp_map_str[seq] = ph

    # Turn into gp graph object (required to get neutral components)
    gp_map = GenotypePhenotypeGraph.read_from_dict(flat_gp_map_str, alphabet=args.alphabet)

    folded_phenotypes = [ph for ph in phenotypes if ph != args.unfolded]
    # get all the neutral components and their genotypes keys by phenotype and
    # neutral component id.
    nc_sizes = gp_map.get_neutral_components(return_boundaries=False, phenotypes=folded_phenotypes)
    count_nc_gt_below_cut = 0
    count_nc_below_cut = 0

    # set all neutral components below the cutoff to unfolded
    for ph in nc_sizes:  # first layer of dictionary contains phenotype keys
        for nc in nc_sizes[ph]:  # loop through all neut. comp. of this phen.
            if len(nc_sizes[ph][nc]) < args.nc_cutoff:  # if below cutoff
                count_nc_below_cut += 1
                for gt in nc_sizes[ph][nc]:  # set all the genotypes as unfold.
                    flat_gp_map_str[gt] = args.unfolded
                    count_nc_gt_below_cut += 1

    
    # turn into ph_to_gt map
    ph_to_gt = dict([(ph, []) for ph in phenotypes])
    for i, gt in enumerate(genotypes):
        ph_to_gt[flat_gp_map_str[gt]].append(i)  # append number genotype to the ph it maps to
    
    # save to file
    dict_to_gpmap(ph_to_gt, args.output)

