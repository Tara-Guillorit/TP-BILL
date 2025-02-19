import argparse
import os
from pathlib import Path
import json
import matplotlib.pyplot as plt
import numpy as np

import utils.plots as plots
from utils.similarity import merge_samples
from utils.stats import pairwise_similarity, count_types, count_types_by_group, len_by_type, unique_variants
from utils.misc import parse_range, build_vcf_path, build_vcf_path_test
from utils.read_vcf import parse_vcf_noerror
from utils.read_ORF import list_interval_with_dico


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Frequency evolution analyse options")

    # Commande :  python3 iteration_stats.py -f 0.15 -d 500 -ec 20000 270000 -lc 3000 -s 1 -sa 1-10 -it 5 -o results/iteration_stats/P90

    parser.add_argument('-f', '--freq', type=float, help="Minimum allele frequency to filter", required=False, default=0.15)
    parser.add_argument('-d', '--depth', type=int, help="Mininum reads depth to filter", required=False, default=400)

    parser.add_argument('-ec', '--edgecut', type=int, nargs='+', help="Bounds to cut huge repetitive insertions", required=False, default=[20000, 270000])
    parser.add_argument('-lc', '--lencut', type=int, help="Maximum variant lenght outside of the edge cut bounds", required=False, default=3000)

    parser.add_argument('-s', '--similarity', type=float, help="Minimum similarity required to group variants", required=False, default=1.)

    parser.add_argument('-sa', '--samples', type=str, help="Samples of interest, range (1 to 10 => 1-10)", required=False, default="1-10")
    parser.add_argument('-it', '--iteration', type=int, help="Iteration to analyse, only one (1=P15 to 5=P90)", required=True)

    parser.add_argument('-o', '--output', type=str, help="Directory to write the output", required=True)

    args = parser.parse_args()
    iterations = [15, 30, 50, 65, 90]
    args.iteration = iterations[args.iteration-1]
    args.samples = parse_range(args.samples)
    args_dict = vars(args).copy()

    args.output = Path(args.output)

    if args.output.exists():
        print("ERROR : le dossier output " + str(args.output) + " existe déjà")
        exit(1)
    else:
        args.output.mkdir(parents=True)

    with open(args.output / "input.json", "w") as fp:
        json.dump(args_dict, fp, indent=4) 

    figs_dir = args.output / "figs"
    figs_dir.mkdir()

    filter_dir = figs_dir / "filter"
    filter_dir.mkdir()

    descript_dir = figs_dir / "descript"
    descript_dir.mkdir()


    ############ Getting data ##############

    sample_labels = [f"P{args.iteration}-{s}" for s in args.samples]
    cold_labels = [f"P{args.iteration}-{s}" for s in args.samples if s <= 5]
    heat_labels = [f"P{args.iteration}-{s}" for s in args.samples if s > 5]

    data = [parse_vcf_noerror(build_vcf_path(s, args.iteration)) for s in args.samples]
    filtered_data = [[v for v in d if (args.edgecut[0] < v['pos'] and v['pos'] < args.edgecut[1]) or abs(v['svlen']) < args.lencut] for d in data]
    filtered_data = [[v for v in d if v['af'] >= args.freq and v['depth'] >= args.depth] for d in data]

    all_variant = sum(data, [])
    all_filtered = sum(filtered_data, [])

    
    ############ Distribution of af, dr, dv ##############

    af_distrib = [v["af"] for v in all_variant]
    af_filtered = [v["af"] for v in all_filtered]

    depth_distrib = [v['depth'] for v in all_variant]
    depth_filtered = [v['depth'] for v in all_filtered]

    plots.distrib_af_depth(af_distrib, af_filtered, depth_distrib, depth_filtered, args.iteration, None, filter_dir)


    ############ Distribution of lengths at genome bounds ##############

    len_by_pos = [(v['pos'], abs(v['svlen'])) for v in all_variant]
    len_by_pos_filtered = [(v['pos'], abs(v['svlen'])) for v in all_filtered]

    len_by_step, len_by_step_filter = plots.len_by_pos(len_by_pos, len_by_pos_filtered, 300000, 1000, 5000, args.iteration, None, filter_dir)


    ############ Filtering and grouping data ##############

    data = filtered_data
    grouped_by_sample = merge_samples(data, sample_labels, sim_thresold=args.similarity)


    ############ Filtering and grouping data ##############

    types_counts = count_types([grp[0][1]["svtype"] for grp in grouped_by_sample])
    types_by_sample = count_types_by_group(grouped_by_sample, sample_labels)
    types_by_choc = count_types_by_group(grouped_by_sample, [cold_labels, heat_labels])

    fig, ax1 = plt.subplots()
    plots.barplot(types_counts, ax1)
    ax1.set_ylabel('Nombre de variants')
    ax1.set_xlabel('Types de variants')
    ax1.set_title(f'Nombre de variant par type dans le passage {args.iteration}')
    plt.savefig(descript_dir / "types.pdf")

    fig, ax2 = plt.subplots()
    plots.grouped_barplot(types_by_sample, sample_labels, ax2)
    ax2.set_ylabel('Nombre de variants')
    ax2.set_xlabel('Échantillons du passage')
    ax2.set_title(f'Nombre de variants par échantillon dans le passage {args.iteration}')
    plt.savefig(descript_dir / "types_by_sample.pdf")

    fig, ax3 = plt.subplots()
    plots.grouped_barplot(types_by_choc, ["Choc Froid", "Choc Chaud"], ax3)
    ax3.set_ylabel('Nombre de variants')
    ax3.set_title(f'Nombre de variants par choc dans le passage {args.iteration}')
    plt.savefig(descript_dir / "type_by_choc.pdf")

    types_dict = {
        "types": types_counts,
        "types_by_sample": types_by_sample,
        "types_by_choc": types_by_choc
    }


    ############ Filtering and grouping data ##############

    len_distrib = len_by_type([(grp[0][1]["svtype"], np.mean([x[1]["svlen"] for x in grp])) for grp in grouped_by_sample])

    fig, ax4 = plt.subplots()
    plots.boxplot(len_distrib, False, ax4)
    ax4.set_ylabel("Longueurs observées")
    ax4.set_title(f"Longueurs des variants dans le passage {args.iteration} (sans outliers)")
    plt.savefig(descript_dir / "lengths.pdf")
    
    fig, ax5 = plt.subplots()
    plots.boxplot(len_distrib, True, ax5)
    ax5.set_ylabel("Longueurs observées")
    ax5.set_title(f"Longueurs des variants dans le passage {args.iteration}")
    plt.savefig(descript_dir / "lengths_outliers.pdf")


    ############ Similarity matrix ##############

    sims = pairwise_similarity(grouped_by_sample, sample_labels)

    fig, ax6 = plt.subplots()
    plots.heatmap(sims, sample_labels, ax6)
    ax6.set_title(f"Proportion de variants partagés entre les échantillons du passage {args.iteration}")
    plt.savefig(descript_dir / "sim_heatmap.pdf")


    ############ Unique variants and the ORF they potentially modify ##############

    uniques = unique_variants(grouped_by_sample, cold_labels, heat_labels, "Cold", "Heat")
    orfs = list_interval_with_dico("ORF.fasta")

    for v in uniques:
        matched = []
        for orf in orfs:
            found = False
            for inter in orf["location"]:
                found = found or (min(v["end"], int(inter[1])) - max(v["pos"], int(inter[0])) > 0)
            if found:
                matched.append(orf["locus_tag"])
        v["orfs"] = matched
    
    
    ############ Write results to json ##############

    tl = [[x if not np.isnan(x) else None for x in row] for row in sims]

    results = {
        "types": types_dict,
        "lengths": len_distrib,
        "pairwise_sim": tl,
        "af": {
            "nofilt": af_distrib,
            "filt": af_filtered
        },
        "depth": {
            "nofilt": depth_distrib,
            "filt": depth_filtered
        },
        "len_by_pos": {
            "nofilt": len_by_step,
            "filt": len_by_step_filter
        }
    }

    results["iteration"] = args.iteration
    results["samples"] = args.samples
    results["similarity"] = args.similarity
    results["filters"] = {
        "freq": args.freq,
        "depth": args.depth,
        "edgecut": args.edgecut,
        "lencut": args.lencut,
    }

    with open(args.output / "results.json", 'w') as fp:
        json.dump(results, fp, indent=4)

    with open(args.output / "uniques.json", 'w') as fp:
        json.dump(uniques, fp, indent=4)
