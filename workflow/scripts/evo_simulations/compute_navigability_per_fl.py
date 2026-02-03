#!/usr/bin/env python

import argparse
import numpy as np
import pickle
from rna_gpf.parsing import load_phenotype_and_metric_from_file
from rna_gpf.analysis import get_peaks



if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-w", "--walk_lengths", help="Files that contain " \
    "adaptive walk lengths for each phenotype. One file per fitness landscape", 
    nargs="+", required=True)       
    parser.add_argument("-f", "--fitness_landscapes", help="Fitness landscapes", nargs="+", required=True)
    parser.add_argument("-o", "--output", help="Output file name",
                        type=str, required=True)

    args = parser.parse_args()


    fl_success_count = {}
    peak_ph = {}
    for i, (path, fl) in enumerate(zip(args.walk_lengths, args.fitness_landscapes)):  # i enumerates fitn. landsc.
        path_count = 0
        fl_success_count[i] = 0
    
        with open(path, "r") as file:
            for j, line_ in enumerate(file):
                path_count += 1
                line = line_.strip().split()
                walk_length = int(line[0])
                if walk_length > 0:  # -1 would indicate unsuccessful walk
                    fl_success_count[i] += 1  # count as success
        fl_success_count[i] /= path_count  # get fraction of successful paths
        
        ph, f = load_phenotype_and_metric_from_file(fl)
        ph_to_f = dict(zip(ph, f))
        peak_ph[i] = ph[np.argmax(f)]

    with open(args.output, "w") as f:
        for i in fl_success_count:
            f.write(f"{peak_ph[i]} {str(fl_success_count[i])} \n")  # new line after every fl block
    
        

                        
                        



    
                    

            
