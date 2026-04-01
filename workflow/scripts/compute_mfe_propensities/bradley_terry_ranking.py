#!/usr/bin/env python

import argparse
import pickle
import numpy as np


from rna_gpf.analysis import infer_bradley_terry_scores

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input consensus matrix, 2D "
                        "numpy array in pickle format", required=True)
    parser.add_argument("-o", "--output", help="Bradley terry scores in csv "
                        "format ", type=str, required=True)
    parser.add_argument("-p", "--phenotypes", help="phenotypes ", required=True)

    args = parser.parse_args()
    
    with open(args.phenotypes, "r") as file:
        phenotypes = np.array([line.strip() for line in file])    

    A = pickle.load(open(args.input, "rb"))
    print(A, np.sum(A), A[np.where(A>0)], flush=True)
    # A = np.array([[0, 2, 0, 1], [3, 0, 5, 0], [0, 3, 0, 1], [4, 0, 3, 0]])

    scores = infer_bradley_terry_scores(A)
    np.set_printoptions(suppress=True)
    
    scores_sorted, phenotypes_ranked = zip(*[(r, p) for r, p in sorted(zip(scores, phenotypes), key=lambda pair: pair[0], reverse=True)])

    with open(args.output, "w") as file:
        for p, r in zip(phenotypes_ranked, scores_sorted):
            file.write(f"{p} {r}\n")
