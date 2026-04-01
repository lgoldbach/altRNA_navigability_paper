#!/usr/bin/env python

import argparse
import pickle
from rna_gpf.parsing import gpmap_to_dict, gpmap_pgdict



if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input phenotype-genotype map",
                        required=True)
    parser.add_argument("-g1", "--genotypes", help="Genotype file",
                        required=True)
    parser.add_argument("-r", "--reference",
                        help="one-to-one reference gp-map", required=True),
    parser.add_argument("-g2", "--genotypes_ref", help="Genotype file",
                        required=True)
    parser.add_argument("-d", "--dead", help="String for the dead phenotype",
                        required=True)
    parser.add_argument("-o1", "--output", help="Dictionary phenotype-genotype",
                        required=True)
    parser.add_argument("-o2", "--output_ref", help="Dictionary ref phenotype-genotype",
                        required=True)

    args = parser.parse_args()    


translation = {"J": "G", "K": "A", "L": "C", "M": "U"}


# map nussinov to AUGC
# translation = {"L": "A", "J": "U", "M": "G", "K": "C"}
translation = {"J": "G", "K": "A", "L": "C", "M": "U"}

D = gpmap_pgdict(gpmap_file=args.input, genotype_file=args.genotypes)
ref_gpmap = gpmap_to_dict(args.reference, args.genotypes_ref)

# Map to AUGC alphabet and remove gt that map to unfolded and ph that don't
# appear in ref. This speeds up the process later on.
ref_phenos = set([i[0] for i in ref_gpmap.values()])
D_augc_folded = {}
for ph in ref_phenos:  # only consider phenotypes from reference GP map
    D_augc_folded[ph] = []
    for g in D[ph]:
        new_g = "".join([translation[s] for s in g])  # map to AUGC alphabet.
        if new_g in ref_gpmap:  # ignore where reference maps to dead phenotype
            D_augc_folded[ph].append(new_g)

# remove entries for dead phenotypes.
for g in ref_gpmap:
    if ref_gpmap[g] == args.dead:
        del ref_gpmap[g]

pickle.dump(D_augc_folded, open(args.output, "wb"))
pickle.dump(ref_gpmap, open(args.output_ref, "wb"))
