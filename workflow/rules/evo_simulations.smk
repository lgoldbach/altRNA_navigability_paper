

rule build_random_fitness_landscape:
    input:
        "<output>/gp_map_viennaRNA/mfe_prop/mfe_prop_ranking.txt"
    output:
        "{path}/fl_{fl_id}/fl.txt"
    params:
        min_f=config["fl_params"]["min_f"],
        max_f=config["fl_params"]["max_f"],
        lethal_ph=config["unfolded"]
    shell:
        "workflow/scripts/evo_simulations/random_fitness_landscape.py "
        "--phenotypes {input} "
        "--low_f {params.min_f} "
        "--upp_f {params.max_f} "
        "--lethal_ph {params.lethal_ph} "
        "-o {output} "

rule run_adaptive_walks_productive_trace_paths:
    input:
        fl="{bp_graph{bp}/ranking{rank}/fitness_landscapes}/fl_{fl_id}/fl.txt",
        gp_map="bp_graph{bp}/ranking{rank}/gp_map.pickle"
    output:
        walk_lengths="bp_graph{bp}/ranking{rank}/fitness_landscapes/fl_{fl_id}/adaptive_walks/pop_size{pop_size}/walk_lengths_seed{seed}.csv",
        paths="bp_graph{bp}/ranking{rank}/fitness_landscapes/fl_{fl_id}/adaptive_walks/pop_size{pop_size}/paths_seed{seed}.csv"
    params:        
        max_steps=None,  # int
        sample_size_walks=None,  # int
        avoid_phenotype=None  # immutable
    shell:
        "workflow/scripts/evo_simulations/adaptive_walks_productive_trace_path.py "
        "-i {input.gp_map} "
        "-f {input.fl} "
        "-l {output.walk_lengths} "
        "-p {output.paths} "
        "-n {wildcards.pop_size} "
        "-m {params.max_steps} "
        "-a {params.avoid_phenotype} "
        "-s {params.sample_size_walks} "
        "-r {wildcards.seed} "