#!/usr/bin/env python

import argparse

from rna_gpf.mapping_functions import gp_mapper, viennaRNA_mfe


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input file with genotypes")
    parser.add_argument("-o", "--output", help="File output for phenotypes")

    args = parser.parse_args()

    # generate g-p map and save to output file
    gp_mapper(input=args.input, output=args.output, 
              mapping_function=viennaRNA_mfe)
