#!/usr/bin/env python

import argparse
import numpy as np
from rna_gpf.parsing import read_nc_to_gt_file


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", help="Output file", type=str)
    parser.add_argument("-f", "--nc_to_gt", help="Neutral network to genotype file", type=int)
    parser.add_argument("-n", "--num", help="Number of genotypes to set unfolded", type=int)

    args = parser.parse_args()

    nc_to_gt = read_nc_to_gt_file(args.nc_to_gt)

    with open(args.output, "w") as f:
        for g in genotypes:
            f.write("".join(g + ("\n",)))
    f.close()
