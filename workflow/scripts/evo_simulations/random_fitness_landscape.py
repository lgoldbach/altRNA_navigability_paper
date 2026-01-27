#!/usr/bin/env python

import argparse
import pickle
import networkx as nx

import matplotlib.pyplot as plt
import datetime
import time
import numpy as np


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--phenotypes", help="Input phenotypes ", 
                        type=str, required=True)
    parser.add_argument("-l", "--low_f", help="Lower fitness limit of an open"
                        "fitness interval, for all f: min_f < f < max_f", 
                        type=float, required=True)
    parser.add_argument("-u", "--upp_f", help="Upper fitness limit of an open"
                        "fitness interval, for all f: min_f < f < max_f", 
                        type=float, required=True)
    parser.add_argument("-d", "--lethal_ph", help="Lethal phenotype to which "
                        "low_f fitness limit will be assigned, i.e. it will "
                        "have the lowest fitness of all phenotypes exlusively", 
                        required=False)      
    parser.add_argument("-g", "--global_peak", help="Phenotype which is going"
                        "to have the highest exclusive fitness (upp_f). ",
                        required=False)                 
    parser.add_argument("-o", "--output", help="Output file name", 
                        type=str, required=True)

    args = parser.parse_args()
    
    phenotypes = []
    with open(args.phenotypes, "r") as f:
        for line in f:
            phenotypes.append(line.strip())
    
    rng = np.random.default_rng()
    ph_to_f = {}
    # randomly assign fitness from the open interval (low_f, upp_f), i.e. for
    # all f: low_f < f < upp_f
    for ph in phenotypes:
        # make sure f is never =low_f because numpy uses half open intervals [,)
        # f=low_f is reserved for lethal phenotypes which is applied below
        f = args.low_f
        while f == args.low_f:
            f = rng.uniform(args.low_f, args.upp_f)
        ph_to_f[ph] = f

    # assign lowest fitness to respective phenotype if applicable
    if args.lethal_ph:
        ph_to_f[args.lethal_ph] = args.low_f

    # check if global peak phenotype was given and assign highest fitness
    if args.global_peak:
        ph_to_f[args.global_peak] = args.upp_f
    # else:  # if not global peak given, pick one at random.
    #     global_peak = args.lethal_ph
    #     while global_peak == args.lethal_ph:  # don't allow lethal ph to be global peak
    #         global_peak = np.random.choice(phenotypes)
    #     ph_to_f[global_peak] = args.upp_f

    with open(args.output, "w") as f:
        for ph in ph_to_f:
            f.write(ph + " " + str(ph_to_f[ph]) + "\n")
