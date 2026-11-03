#!/usr/bin/env python

import argparse

from rna_gpf.gp_map import GenotypePhenotypeGraph
import pickle

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Input genotype-phenotype map "
                        "file", required=True)
    parser.add_argument("-g", "--genotypes", help="Reference file that "
                        "contains the ith genotype on the ith line. Used to "
                        "map the numeric genotypes back to the sequences ", 
                        required=True)
    parser.add_argument("-a", "--alphabet", help="Define the alphabet used in "
                        "the genotype map, e.g. -a AUGC", required=True)
    parser.add_argument("-i", "--ignore_phenotype", help="Phenotype to ignore", required=False, default=None)
    parser.add_argument("-o", "--output", help="File output for robustness",
                        required=True)
    

    args = parser.parse_args()
    
    gpm = GenotypePhenotypeGraph.read_from_ph_to_gt_file(path=args.file, 
                                              genotype_ref_path=args.genotypes,
                                              alphabet=args.alphabet,
                                              ignore_phenotype=args.ignore_phenotype)

    pickle.dump(gpm, open(args.output, "wb"))
