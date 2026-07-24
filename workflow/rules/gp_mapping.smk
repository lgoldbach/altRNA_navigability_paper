import os

rule generate_alt_genotype_space:
    output:
        "<output>/alt_genotypes.txt"
    params:
        alphabet=config["alt_alphabet"],
        seq_len=config["sequence_length"]
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/gp_mapping/build_genotype_space.py " 
        "-a {params.alphabet} "
        "-l {params.seq_len} "
        "-o {output}"

rule generate_viennaRNA_genotype_space:
    output:
        "<output>/gp_map_viennaRNA/genotypes.txt"
    params:
        alphabet=config["viennaRNA"]["alphabet"],
        seq_len=config["sequence_length"]
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/gp_mapping/build_genotype_space.py " 
        "-a {params.alphabet} "
        "-l {params.seq_len} "
        "-o {output}"

rule make_permissible_sets:
    input:
        "<output>/alt_genotypes.txt"
    output:
        "<output>/gp_map_{bp}/permissible_sets.txt"  # regex to constrain to numbers and prevent matching subdirectory gp_map.txt
    params:
        graph_path=config["bp_rule_graphs"],
        min_loop_size=config["folding_params"]["min_loop"],
        suboptimal=config["folding_params"]["hbonds_below_max"],
        structures_max=config["folding_params"]["permis_set_max"],
        alphabet=config["alt_alphabet"]
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"]
    shell:
        "workflow/scripts/gp_mapping/generate_permissible_sets.py "
        "-i {input} "
        "-o {output} "
        "-m {params.min_loop_size} "
        "-s {params.suboptimal} "
        "-z {params.structures_max} "
        "-g {params.graph_path} "
        "-a {params.alphabet} "
        "-p {wildcards.bp} "

rule extract_phenotype_list:
    input:
        "<output>/gp_map_{bp}/gp_map.pickle"
    output:
        "<output>/gp_map_{bp}/phenotypes.txt"
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    wildcard_constraints:
        bp='[0-9]*'
    run:
        import pickle

        gpm = pickle.load(open(str(input), "rb"))
        with open(str(output), "w") as f:
            for ph in gpm.phenotype_set:
                f.write(ph + "\n")

rule extract_phenotype_list_ref_unfolded:
    input:
        "<output>/gp_map_{bp}/gp_map_ref_unfolded.pickle"
    output:
        "<output>/gp_map_{bp}/phenotypes_ref_unfolded.txt"
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    wildcard_constraints:
        bp='[0-9]*'
    run:
        import pickle

        gpm = pickle.load(open(str(input), "rb"))
        with open(str(output), "w") as f:
            for ph in gpm.phenotype_set:
                f.write(ph + "\n")

rule extract_phenotype_list_random_unfolded:
    input:
        "<output>/gp_map_{bp}/gp_map_unfolded{n_unfolded}.pickle"
    output:
        "<output>/gp_map_{bp}/phenotypes_unfolded{n_unfolded}.txt"
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    wildcard_constraints:
        bp='[0-9]*'
    run:
        import pickle

        gpm = pickle.load(open(str(input), "rb"))
        with open(str(output), "w") as f:
            for ph in gpm.phenotype_set:
                f.write(ph + "\n")        

rule extract_phenotype_list_nc_unfolded:
    input:
        "<output>/gp_map_{bp}/gp_map_nc_unfolded_n_unf{n_unfolded}.pickle",
    output:
        "<output>/gp_map_{bp}/phenotypes_nc_unfolded_n_unf{n_unfolded}.txt",
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    run:
        import pickle

        gpm = pickle.load(open(str(input), "rb"))
        with open(str(output), "w") as f:
            for ph in gpm.phenotype_set:
                f.write(ph + "\n")   


rule pick_random_unfolded_genotypes:
    output:
        "<output>/unfolded_genotypes_{n_unfolded}.txt"
    input:
        "<output>/alt_genotypes.txt"
    params:
        total_seqs = len(config["alt_alphabet"]) ** config["sequence_length"]
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    run:
        import numpy as np
        gt_list = np.arange(0, params["total_seqs"])
        rand_g = np.random.choice(gt_list, size=int(wildcards["n_unfolded"]), replace=False)
        with open(output[0], "w") as f:
            for g in rand_g:
                f.write(str(g) + "\n")


rule flatten_gp_map_nc_cutoff_random_unfolded:
    input:
        gp_map="<output>/gp_map_{bp}/permissible_sets.txt",
        genotypes="<output>/alt_genotypes.txt",
        ranking="<output>/gp_map_" + config["mfe_prop_ref_gp_map"] + "/mfe_propensities.csv",
        unfolded_genotypes="<output>/unfolded_genotypes_{n_unfolded}.txt",
    output:
        "<output>/gp_map_{bp}/gp_map_unfolded{n_unfolded}.txt"
    params:
        nc_size_cutoff=config["folding_params"]["nc_size_cutoff"],
        unfolded=config["folding_params"]["unfolded"],
        alphabet=config["alt_alphabet"]
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/gp_mapping/flatten_gp_map_nc_cutoff.py "
        "-i {input.gp_map} "
        "-r {input.ranking} "
        "-g {input.unfolded_genotypes} "
        "-t {input.genotypes} "
        "-c {params.nc_size_cutoff} "
        "-u {params.unfolded} "
        "-a {params.alphabet} "
        "-o {output} "

rule flatten_gp_map_nc_cutoff_ref_unfolded:
    input:
        gp_map="<output>/gp_map_{bp}/permissible_sets.txt",
        genotypes="<output>/alt_genotypes.txt",
        ranking="<output>/gp_map_" + config["mfe_prop_ref_gp_map"] + "/mfe_propensities.csv",
        unfolded_genotypes="<output>/ref_unfolded_genotypes.txt",
    output:
        "<output>/gp_map_{bp}/gp_map_ref_unfolded.txt"
    params:
        nc_size_cutoff=config["folding_params"]["nc_size_cutoff"],
        unfolded=config["folding_params"]["unfolded"],
        alphabet=config["alt_alphabet"]
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/gp_mapping/flatten_gp_map_nc_cutoff.py "
        "-i {input.gp_map} "
        "-r {input.ranking} "
        "-g {input.unfolded_genotypes} "
        "-t {input.genotypes} "
        "-c {params.nc_size_cutoff} "
        "-u {params.unfolded} "
        "-a {params.alphabet} "
        "-o {output} "

rule flatten_gp_map_nc_unfolded:
    input:
        gp_map="<output>/gp_map_{bp}/permissible_sets.txt",
        genotypes="<output>/alt_genotypes.txt",
        ranking="<output>/gp_map_" + config["mfe_prop_ref_gp_map"] + "/mfe_propensities.csv",
        unfolded_genotypes="<output>/gp_map_{bp}/unfolded_genotypes_nc_sample_n_unf{n_unfolded}.txt",
    output:
        "<output>/gp_map_{bp}/gp_map_nc_unfolded_n_unf{n_unfolded}.txt"
    params:
        nc_size_cutoff=config["folding_params"]["nc_size_cutoff"],
        unfolded=config["folding_params"]["unfolded"],
        alphabet=config["alt_alphabet"]
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/gp_mapping/flatten_gp_map_nc_cutoff.py "
        "-i {input.gp_map} "
        "-r {input.ranking} "
        "-g {input.unfolded_genotypes} "
        "-t {input.genotypes} "
        "-c {params.nc_size_cutoff} "
        "-u {params.unfolded} "
        "-a {params.alphabet} "
        "-o {output} "


rule build_gpmap_python_object_nc_unfolded:
    input:
        gp_map="<output>/gp_map_{bp}/gp_map_nc_unfolded_n_unf{n_unfolded}.txt",
        genotypes="<output>/alt_genotypes.txt"
    output:
        "<output>/gp_map_{bp}/gp_map_nc_unfolded_n_unf{n_unfolded}.pickle",
    params:
        alphabet=config["alt_alphabet"],
        ignore_phenotype=config["unfolded"]  # ignore unviable phenotype
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/gp_mapping/build_gpmap_pickle.py "
        "-f {input.gp_map} "
        "-g {input.genotypes} "
        "-a {params.alphabet} "
        "-i {params.ignore_phenotype} "
        "-o {output}"


rule build_gpmap_python_object_ref_unfolded:
    input:
        gp_map="<output>/gp_map_{bp}/gp_map_ref_unfolded.txt",
        genotypes="<output>/alt_genotypes.txt"
    output:
        "<output>/gp_map_{bp}/gp_map_ref_unfolded.pickle",
    params:
        alphabet=config["alt_alphabet"],
        ignore_phenotype=config["unfolded"]  # ignore unviable phenotype
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/gp_mapping/build_gpmap_pickle.py "
        "-f {input.gp_map} "
        "-g {input.genotypes} "
        "-a {params.alphabet} "
        "-i {params.ignore_phenotype} "
        "-o {output}"

rule flatten_gp_map_no_unfolded:
    input:
        gp_map="<output>/gp_map_{bp}/permissible_sets.txt",
        genotypes="<output>/alt_genotypes.txt",
        ranking="<output>/gp_map_" + config["mfe_prop_ref_gp_map"] + "/mfe_propensities.csv",
    output:
        "<output>/gp_map_{bp}/gp_map_no_added_unfolded.txt"
    params:
        alphabet=config["alt_alphabet"],
        unfolded=config["folding_params"]["unfolded"],
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/gp_mapping/flatten_gp_map_no_added_unfolded.py "
        "-i {input.gp_map} "
        "-r {input.ranking} "
        "-t {input.genotypes} "
        "-u {params.unfolded} "
        "-a {params.alphabet} "
        "-o {output} "


rule build_gpmap_python_object_random_unfolded:
    input:
        gp_map="<output>/gp_map_{bp}/gp_map_unfolded{n_unfolded}.txt",
        genotypes="<output>/alt_genotypes.txt"
    output:
        "<output>/gp_map_{bp}/gp_map_unfolded{n_unfolded}.pickle"
    params:
        alphabet=config["alt_alphabet"],
        ignore_phenotype=config["unfolded"]  # ignore unviable phenotype
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/gp_mapping/build_gpmap_pickle.py "
        "-f {input.gp_map} "
        "-g {input.genotypes} "
        "-a {params.alphabet} "
        "-i {params.ignore_phenotype} "
        "-o {output}"

rule build_gpmap_python_object_no_added_unfolded:
    input:
        gp_map="<output>/gp_map_{bp}/gp_map_no_added_unfolded.txt",
        genotypes="<output>/alt_genotypes.txt"
    output:
        "<output>/gp_map_{bp}/gp_map_no_added_unfolded.pickle"
    params:
        alphabet=config["alt_alphabet"],
        ignore_phenotype=None  # ignore unviable phenotype
    resources:
        mem_mb_per_cpu=60000,
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/gp_mapping/build_gpmap_pickle.py "
        "-f {input.gp_map} "
        "-g {input.genotypes} "
        "-a {params.alphabet} "
        "-i {params.ignore_phenotype} "
        "-o {output}"

rule viennaRNA_mfe_gp_map:
    input:
        "<output>/gp_map_viennaRNA/genotypes.txt"
    output:
        "<output>/gp_map_viennaRNA/gp_map.txt"
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/gp_mapping/viennaRNAmfe_gp_mapping.py "
        "-i {input} "
        "-o {output} "

rule compute_phenotype_distribution:
    input:
        gp_map="{path}/gp_map.txt"
    output:
        "{path}/phenotype_distribution.txt"
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "compute_phenotype_distribution.py "
        "-f {input} "
        "-o {output} "

rule build_nc_graph:
    input:
        "{path}/gp_map.pickle"
    output:
        nc_graph="{path}/nc_graph.pickle",
        nc_to_genotype="{path}/nc_to_gt.txt"
    params:
        ignore=config["unfolded"]
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/gp_mapping/nc_graph.py "
        "-f {input} "
        "-i {params.ignore} "
        "-g {output.nc_graph} "
        "-m {output.nc_to_genotype} "


rule build_nc_graph_no_added_unfolded:
    input:
        "<output>/gp_map_{bp}/gp_map_no_added_unfolded.pickle"
    output:
        nc_graph="<output>/gp_map_{bp}/nc_graph_no_added_unfolded.pickle",
        nc_to_genotype="<output>/gp_map_{bp}/nc_to_gt_no_added_unfolded.txt",
    params:
        ignore=None
    resources:
        mem_mb_per_cpu=60000,
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/gp_mapping/nc_graph.py "
        "-f {input} "
        "-i {params.ignore} "
        "-g {output.nc_graph} "
        "-m {output.nc_to_genotype} "

rule nc_graph_to_nc_sizes_txt:
    input:
        "<output>/gp_map_{bp}/nc_graph.pickle"
    output:
        "<output>/gp_map_{bp}/neutral_component_sizes.txt"
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    run:
        import pickle
        # load nc graph and write nc sizes into a txt file
        nc_graph = pickle.load(open(input[0], "rb"))
        with open(output[0], "w") as f:
            for node in nc_graph:
                f.write(f"{str(node)} {str(nc_graph.nodes[node]['size'])}\n")   

rule sample_random_ncs_for_unfolded:
    input: 
        nc_to_gt="<output>/gp_map_{bp}/nc_to_gt_no_added_unfolded.txt",
    output:
        "<output>/gp_map_{bp}/unfolded_genotypes_nc_sample_n_unf{n_unfolded}.txt"
    params:
        unfolded_num=config["folding_params"]["n_unfolded"]
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/gp_mapping/sample_unfolded_nc.py "
        "-i {input} "
        "-n {params.unfolded_num} "
        "-o {output} "
 