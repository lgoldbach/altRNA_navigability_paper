#!/usr/bin/env python

import argparse
import pickle
import networkx as nx
import numpy as np
import datetime
import time

from rna_gpf.adaptive_walks import kimura_fixation_from_fitness, pairwise_transition_prob_dict, kimura_fixation, productive_adaptive_walk_w_T


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--gp_map", help="Input gp map.", required=True)
    parser.add_argument("-f", "--fl", help="Fitness landscape "
                        "file", required=True)
    parser.add_argument("-s", "--sample_size_walks", help="How many random walks to sample ",
                        type=int, required=False)
    parser.add_argument("-n", "--population_size", help="Populaton size", 
                        type=int, required=False)
    parser.add_argument("-m", "--max_steps", help="Maximum number of steps "
                        "per walk", type=int, required=False)
    parser.add_argument("-r", "--seed", help="random seed", type=int,
                        required=False)
    parser.add_argument("-a", "--avoid", help="Phenotypes to avoid as starting "
                        "nodes", required=False)
    parser.add_argument("-l", "--walk_lengths", help="file for output data",
                        required=True)
    parser.add_argument("-p", "--paths", help="File where all paths are saved", required=True)
    
    args = parser.parse_args()

    if args.seed:
        rng = np.random.default_rng(seed=args.seed)
    else:
        rng = np.random.default_rng(seed=1996)

    G = pickle.load(open(args.gp_map, "rb"))

    phenotypes = G.phenotype_set
    
    # load fl. There is no target phenotype yet. All phenotype have fitness
    # in the interval [0, 1)
    ph_to_fitness_original = {}
    with open(args.fl, "r") as f:
            for line in f:
                data = line.strip().split(" ")
                phenotype = data[0]
                fitness = float(data[1])
                ph_to_fitness_original[phenotype] = fitness
    
    # every phenotype is the target at least once
    for target_ph in phenotypes:
        print(target_ph, flush=True)
        # read in fitness landscape
        ph_to_fitness = ph_to_fitness_original.copy()
        ph_to_fitness[target_ph] = 1  # set fitness of target to 1

        print("A", flush=True)
        # fixartion probability
        fix_prob = lambda x, y: kimura_fixation_from_fitness(x, y, N=args.population_size)
        print("B", flush=True)

        a = datetime.datetime.now()
        # precompute transition probabilities between any pairs of phenotypes
        T = pairwise_transition_prob_dict(f_map=ph_to_fitness, func=fix_prob)
        b = datetime.datetime.now()
        c = b-a
        print(c.seconds/60, c.seconds, flush=True)

        all_nodes = set(G.genotypes)
        # get list of target nodes. Do not start walks from there (would be redundant)
        non_starting_nodes = [g for i, g in enumerate(G.genotypes) if G.phenotypes[i]==target_ph]

        # phenotypes to avoid as starting nodes, e.g. lethal ones.
        if args.avoid:
            # get all lethal nodes
            lethal_nodes = [g for i, g in enumerate(G.genotypes) if G.phenotypes[i]==args.avoid]
            non_starting_nodes += lethal_nodes  # also disallow lethal nodes as starting nodes
        # extract nodes that are not target (or lethal if applicable)
        potential_starting_nodes = list(all_nodes.difference(set(non_starting_nodes)))
        
        start_gt = rng.choice(potential_starting_nodes, size=min(args.sample_size_walks, len(potential_starting_nodes)), replace=False)

        print(start_gt, flush=True)
        adaptive_walk_lengths = []  # store adaptive walk lenghts for each phenotype
        paths = []  # store whole paths of genotypes 

        for g in start_gt:
            # store adaptive walks by target phenotype
            path = productive_adaptive_walk_w_T(G, g,
                                    fitness_function=ph_to_fitness, 
                                    T=T,
                                    max_steps=args.max_steps,
                                    rng=rng)
            
            if ph_to_fitness[G.map(path[-1])] == 1:  # walk reached target
                adaptive_walk_lengths.append(len(path))
            else:
                adaptive_walk_lengths.append(-1)  # walk didn't reach target
    
            paths.append(path)  # save path
        

        # Write adaptive walk path each paths of genotypes into a single line each 
        with open(args.paths, "a") as file:
            file.write(target_ph + "\n")
            for path in paths:
                file.write(" ".join(path) + "\n")
        
        # record the adaptive walk lenghts
        with open(args.walk_lengths, "a") as file:
            file.write(target_ph + "\n")
            for lengths in adaptive_walk_lengths:
                file.write(str(lengths) + "\n")

        