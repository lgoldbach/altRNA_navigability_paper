#!/usr/bin/env python

import argparse
import numpy as np

from rna_folding.mapping_functions import viennaRNA_mfe
from rna_folding.parsing import genotype_file_to_numpy


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input file with genotypes")
    parser.add_argument("-o", "--output", help="File output for phenotypes")

    args = parser.parse_args()

    genotypes = genotype_file_to_numpy(args.input)
    mfes = {}

    # get mfe values for all mfe phenotypes and store in dict keyed by phenotypes
    for gt in genotypes:
        gt = np.random.choice(genotypes)
        ph, mfe = viennaRNA_mfe(gt, return_mfe=True)
        if "(" in ph:  # only store mfe of phenotypes with at least one bp
            mfe = np.round(mfe, 4)
            if ph in mfes:
                mfes[ph].append(mfe)
            else:
                mfes[ph] = [mfe]

    with open(args.output, "w") as out:
        for ph in mfes:
            out.write(ph + " " + str(np.mean(mfes[ph])) + " ".join([str(num) for num in mfes[ph]]) + "\n")

