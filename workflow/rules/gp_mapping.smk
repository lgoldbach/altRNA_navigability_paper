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
        "workflow/scripts/gp_map_{bp}/permissible_sets.txt"  # regex to constrain to numbers and prevent matching subdirectory gp_map.txt
    params:
        graph_path=config["bp_rule_graphs"],
        min_loop_size=config["folding_params"]["min_loop"],
        suboptimal=config["folding_params"]["hbonds_below_max"],
        structures_max=config["folding_params"]["permis_set_max"],
        alphabet=config["alt_alphabet"]
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
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
    run:
        import pickle

        gpm = pickle.load(open(str(input), "rb"))
        with open(str(output), "w") as f:
            for ph in gpm.phenotype_set:
                f.write(ph + "\n")

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
        nc_to_genotype="{path}/nc_to_genotype.txt"
    params:
        ignore=config["unfolded"]
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "nc_graph.py "
        "-f {input} "
        "-i {params.ignore} "
        "-g {output.nc_graph} "
        "-m {output.nc_to_genotype} "

rule nc_graph_to_nc_sizes_txt:
    input:
        "{path}/nc_graph.pickle"
    output:
        "{path}/neutral_component_sizes.txt"
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    run:
        import pickle
        # load nc graph and write nc sizes into a txt file
        nc_graph = pickle.load(open(input[0], "rb"))
        with open(output[0], "w") as f:
            for node in nc_graph:
                f.write(str(nc_graph.nodes[node]["size"]) + "\n")   
