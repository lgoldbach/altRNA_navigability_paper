#!/usr/bin/env python

import argparse
import numpy as np
import pickle
from rna_gpf.parsing import load_phenotype_and_metric_from_file
from rna_gpf.analysis import get_peaks



if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--nc_graph", help="Neutral component graph (pickle)", required=True)       
    parser.add_argument("-f", "--fitness_landscapes", help="Fitness landscapes", nargs="+", required=True)
    parser.add_argument("-o", "--output", help="Output file name",
                        type=str, required=True)

    args = parser.parse_args()

    nc_graph = pickle.load(open(args.nc_graph, "rb"))

    peak_ph = {}
    ruggedness = {}

    for fl in args.fitness_landscapes:
        ph, f = load_phenotype_and_metric_from_file(fl)
        ph_to_f = dict(zip(ph, f))
        peak_ph[fl] = ph[np.argmax(f)]

        ruggedness[fl] = 0
        peaks_nc, peaks_f = get_peaks(nc_graph, ph_to_f)
        for nc in peaks_nc:
            ruggedness[fl] += nc_graph.nodes[nc]["size"]

    
    with open(args.output, "w") as f:
        for i in ruggedness:
            f.write(f"{peak_ph[i]} {str(ruggedness[i])} \n")  # new line after every fl block
    
        

                        
                        



    
                    

            
