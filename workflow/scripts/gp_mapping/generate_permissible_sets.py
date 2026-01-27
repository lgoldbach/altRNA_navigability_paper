#!/usr/bin/env python

import argparse

from rna_gpf.base_pairing import BasePairing
from rna_gpf.mapping_functions import gp_mapper, nussinov


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input file with genotypes")
    parser.add_argument("-o", "--output", help="File output for phenotypes")
    parser.add_argument("-m", "--min_loop_size", required=True, type=int, default=1, help="Minimum size for loop")
    parser.add_argument("-s", "--suboptimal", type=int, required=True,
                        help="Create all suboptimal structures with number of base-pairs in the range"
                             "of max - s, where s is an integer. Without the flag, only one structure "
                             "is computed.")
    parser.add_argument("-z", "--structures_max", required=False, type=int, help="Limit on how many suboptimal structures to generate")
    parser.add_argument("-p", "--base_pairing", required=False, type=int, default=-1, help="Which base-pairing to choose. I.e. from the base-pairing simple graphs, which one to pick "
                        "e.g. for 4 bases there are 11 possible base-pairings, so possible input is any number between 1 and 11, If given -1 then it uses canonical base-pairing and AUGC bases")
    parser.add_argument("-a", "--alphabet", required=False, type=str, default="AUGC", help="Which bases do the genotypes contain, e.g. 'AUGC' for canonical RNA")
    parser.add_argument("-g", "--graph_path", required=True, type=str,
                        help="Path to folder containing the base-pairing "
                        "graphs files, e.g. graph4.adj. Check base_pairing.py "
                        "for info on where these graphs come from.")

    args = parser.parse_args()
    
    pairing = BasePairing(bases=args.alphabet,
                          graph_path=args.graph_path, 
                          id=args.base_pairing)
    
    mapping = lambda seq: nussinov(seq, 
                                   base_pairing=pairing, 
                                   min_loop_size=args.min_loop_size, 
                                   suboptimal=args.suboptimal,
                                   structures_max=args.structures_max)

    # generate g-p map and save to output file
    gp_mapper(input=args.input, output=args.output, 
              mapping_function=mapping)
