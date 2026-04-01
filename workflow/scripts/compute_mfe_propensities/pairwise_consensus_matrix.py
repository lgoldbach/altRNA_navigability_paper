#!/usr/bin/env python

import argparse
import pickle
from rna_gpf.analysis import pairwise_consensus_matrix



if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input phenotype-genotype map",
                        required=True)
    parser.add_argument("-r", "--reference",
                        help="One-to-one reference gp-map", required=True),
    parser.add_argument("-o", "--output", help="Consensus matrix",
                        required=True)
    parser.add_argument("-p", "--phenotypes", help="List of phenotypes in "
                        "order corresponding to consesnus matrix indices",
                        required=True)

    args = parser.parse_args()    

    pg_map = pickle.load(open(args.input, "rb"))
    gp_map_ref = pickle.load(open(args.reference, "rb"))

    # get unique list of phenotype from the reference gp map
    # the values are lists that contain only one phenotype, hence the 
    # "complicated" approach
    phenotypes = list(set([i[0] for i in gp_map_ref.values()]))
    
    phenotypes.remove("............")

    A = pairwise_consensus_matrix(phenotypes=phenotypes,
                                  pg_map=pg_map, 
                                  ref_gp_map=gp_map_ref)
    
    pickle.dump(A, open(args.output, "wb"))

    with open(args.phenotypes, "w") as file:
        for ph in phenotypes:
            file.write(ph + "\n")
