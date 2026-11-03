# Code for "The promiscuity and degeneracy of ten unique RNA base-pairing rules can determine the success of adaptive RNA evolution" by Goldbach, Martinm, Wagner 2026 (doi.org/10.1371/journal.pcbi.1014914).

This repository contains all the code and data needed for the modeling, simulations and data analysis performed in the paper.

### Requirements

All requirements can be found in the **requirements.txt**, are compatible with **Python version 3.11** and should be installable with standard package managers. The only package you should have to **install manually** is [rna_gpf](https://github.com/lgoldbach/rna_gpf), which is contains the basic functionality for dealing with RNA secondary structure genotype-phenotype-fitness maps.

### Project structure and workflow

All operations are carried out as Snakemake workflows. Snakemake is a popular declarative workflow manager. The Snakefile (*workflow/Snakefile*) declares which files should be produced and the the rules (*workflow/rules/*) contain the input-output logic for producing these files step-by-step using the associated Python scripts (*workflow/scripts/*) and resources (*resources/*). The only resources ("data") needed for this project are the adjacency graphs of the 10 RNA base-pairing rules. The configurations file (*config/config.yaml*) contains all the parameters for this project, like the RNA sequence length, RNA folding constraints, number of fitness landscapes to be sampled, etc..

### Building the GP maps

For simplicity, the current main Snakefile (workflow/Snakefile) only instructs the generation of the GP maps of the ten unique base-pairing rules and the reference ViennaRNA GP map. Code for further analyses is contained in the additional rules and scriptes inside "workflow/". To create the 11 GP maps, run a suitable snakemake command from the root, e.g.: 

`snakemake --executor slurm --jobs 11 --latency-wait 60` 

for execution on a compute cluster, which will run the workflow executing 11 jobs in parallel at a time (latency-wait is set to 60 seconds to account for the file-system latency that some clusters have). Note that we have declared cluster resource requirements for all Snakemake rules. These resources may or may not be compatible with your cluster. Note that the permissible sets that will be produced each require 2-20 GB of memory.

If you are unfamiliar with Snakemake, we recommend reading up on its documentation first and starting with running smaller subsets of the workflow first.

Feel free to reach out to Leander Goldbach if you have any questions.
