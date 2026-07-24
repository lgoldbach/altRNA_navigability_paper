#!/usr/bin/env python

import argparse
import pickle
import networkx as nx
import matplotlib.pyplot as plt
import datetime
import time
import numpy as np

from rna_gpf.parsing import gpmap_to_dict



if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Input genotype-phenotype map "
                        "file", required=True)
    parser.add_argument("-i", "--ignore", help="Phenotype to ignore, e.g "
                        "unfolded", type=str, required=False)
    parser.add_argument("-g", "--nc_graph", help="File output for neutral "
                        "component graph. Should end in .pickle",
                        required=True)
    parser.add_argument("-m", "--nc_to_g", help="Output file for neutral component"
                        "to genotype map", required=True)
    

    args = parser.parse_args()
    
    a = datetime.datetime.now()
    print("Start loading gpm", a, flush=True)
    gpm = pickle.load(open(args.file, "rb"))
    b = datetime.datetime.now()
    c = b-a
    print(f"Done in {np.round(c.seconds/3600, 2), np.round(c.seconds/60, 2), c.seconds}, {b}", flush=True)

    if args.ignore:
        phenotypes = [ph for ph in gpm.phenotype_set if ph != args.ignore]
        print(len(gpm.phenotype_set), args.ignore, type(args.ignore), type(gpm.phenotype_set[0]))
    else:
        phenotypes = None
    print(len(phenotypes), flush=True)
    a = datetime.datetime.now()
    print("Starting neutral_components", a, flush=True)
    ncs, boundaries = gpm.get_neutral_components(phenotypes=phenotypes,
                                                return_boundaries=True)

    b = datetime.datetime.now()
    c = b-a
    print(f"Done in {np.round(c.seconds/3600, 2), np.round(c.seconds/60, 2), c.seconds}, {b}", flush=True)

    print(f"Number of nc: {len([nc for ph in ncs for nc in ncs[ph]])}", flush=True)

    a = datetime.datetime.now()
    print("Starting building nc graph", a, flush=True)
    nc_graph = nx.Graph()
    # add a node for each neutral component and its phenotype value
    for ph in ncs:
        for nc in ncs[ph]:
            nc_graph.add_node(nc, phenotype=ph, size=len(ncs[ph][nc]))

    for (i, j) in boundaries:
        try:
            nc_i = gpm.neutral_components[i] # get nc id
            nc_j = gpm.neutral_components[j]  # get nc id

            # edge weights are added in .5 increments because they will be
            # double counted, once from each direction
            if not nc_graph.has_edge(nc_i, nc_j):
                nc_graph.add_edge(nc_i, nc_j, weight=.5)  # add edge to nc graph
            else:
                # increase weight if the edge already exists. In undirected graph 
                # the order of nodes does not matter when referencing the edge,
                # (ni, nj) or (nj, ni) are the same edge
                nc_graph.edges[(nc_i, nc_j)]["weight"] += .5
        except KeyError:
            # this is a boundary to a genotype that was marked as "ignored",
            # e.g the unfolded phenotype
            pass

    # b = datetime.datetime.now()
    # c = b-a
    # print(f"Done in {np.round(c.seconds/3600, 2), np.round(c.seconds/60, 2), c.seconds}, {b}", flush=True)
        
    # a = datetime.datetime.now()
    print("Dump nc graph", a, flush=True)
    pickle.dump(nc_graph, open(args.nc_graph, "wb"))

    with open(args.nc_to_g, "w") as f:
        for ph in ncs:
            for nc in ncs[ph]:
                f.write(str(nc) + " " + str(ph) + " " + " ".join(ncs[ph][nc]) + "\n")

    # b = datetime.datetime.now()
    # c = b-a
    # print(f"Done in {np.round(c.seconds/3600, 2), np.round(c.seconds/60, 2), c.seconds}, {b}", flush=True)
    # # arr = nx.to_numpy_array(ph_graph)
    # # attr = nx.get_node_attributes(ph_graph, "phenotype")
    # # pickle.dump((arr, attr), open(args.output, "wb"))

    # labels = {node: str(node) + " " + nc_graph.nodes[node]["phenotype"] for node in nc_graph}
    # pos = nx.spring_layout(nc_graph)
    
    # weights = [nc_graph[u][v]['weight']*3 for u,v in nc_graph.edges()]

    # color_by_ph = {ph: np.random.choice(range(256), size=3)/256 for ph in ncs}
    # colors = [color_by_ph[nc_graph.nodes[node]["phenotype"]] for node in nc_graph.nodes]

    # nx.draw(nc_graph, pos=pos, labels=labels, node_color=colors, node_size=[s*1000 for s in nx.get_node_attributes(nc_graph, "size").values()], width=weights)

    # edge_labels = {(u, v): int(nc_graph[u][v]['weight']) for u,v in nc_graph.edges()}
    # nx.draw_networkx_edge_labels(nc_graph, pos=pos, edge_labels=edge_labels)
    # plt.savefig("neutral_component_graph.pdf", dpi=30)
    
    # plt.clf()

    # gpm.add_hamming_edges()
    # labels = nx.get_node_attributes(gpm, 'neutral_component')
    # nx.draw(gpm, labels=labels)
    # plt.savefig("neutral_components.pdf", dpi=30)
    
    # plt.clf()

    # labels = {node: node + " " + str(gpm.nodes[node]["neutral_component"]) + " " + gpm.nodes[node]["phenotype"] for node in gpm if "neutral_component" in gpm.nodes[node]}
    # nx.draw(gpm, labels=labels)
    # plt.savefig("gpm_ph.pdf", dpi=30)
   