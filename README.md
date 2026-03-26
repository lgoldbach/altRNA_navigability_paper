# Code for paper: Goldbach, Martin, Wagner 2026

This repository contains all the code and data needed for the modeling, simulations and data analysis performed in the paper.


### Requirements

All requirements can be found in the **requirements.txt**, are compatible with **Python version 3.11** and should be installable with standard package managers. The only package you should have to **install manually** is [rna_gpf](https://github.com/lgoldbach/rna_gpf), which is contains the basic functionality for dealing with RNA secondary structure genotype-phenotype-fitness maps.


### Project structure and workflow

All operations are carried out as part of one big Snakemake workflow. Snakemake is a popular declarative workflow manager. The Snakefile (*workflow/Snakefile*) declares which files should be produced and the the rules (*workflow/rules/*) contain the input-output logic for producing these files step-by-step using the associated Python scripts (*workflow/scripts/*) and resources (*resources/*). The only resources ("data") needed for this project are the adjacency graphs of the 10 RNA base-pairing rules. The configurations file (*config/config.yaml*) contains all the parameters for this project, like the RNA sequence length, RNA folding constraints, number of fitness landscapes to be sampled, etc..


### Executing the workflow and dealing with large files

The whole workflow takes a long time and includes thousands of separate operations, some of which require a lot of memory (>10GB). We have executed the workflow on a slurm cluster with many nodes (many of which have sufficiently high memory). If you want to rerun the whole workflow, we recommend doing the same, using for example:

`snakemake --executor slurm  --jobs 50 --latency-wait 60`

which will run the workflow executing 50 jobs in parallel at a time (latency-wait is set to 60 seconds to account for the file-system latency that some clusters have). Note that we have declared cluster resource requirements for all Snakemake rules. These resources may or may not be compatible with your cluster.


If you are unfamiliar with Snakemake, we recommend reading up on its documentation first and starting with running smaller subsets of the workflow first (A small test-workflow will be added to the repo soon).

Feel free to reach out to Leander Goldbach if you have any questions.
