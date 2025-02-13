import argparse
import os
from pathlib import Path
import json
import matplotlib.pyplot as plt
import numpy as np

import utils.plots as plots
from utils.similarity import merge_samples
from utils.stats import pairwise_similarity, count_types, count_types_by_group, len_by_type
from utils.misc import parse_range, build_vcf_path, build_vcf_path_test
from utils.read_vcf import parse_vcf_noerror
from utils.read_ORF import list_interval_with_dico


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Frequency evolution analyse options")

    # Commande :  python3 sample_stats.py -f 0.1 -d 250 -ec 10000 270000 -lc 2000 -s 1 -sa 10 -it 1-5 -o results/sample_stats/10

    parser.add_argument('-f', '--freq', type=float, help="Minimum allele frequency to filter", required=False, default=0.05)
    parser.add_argument('-d', '--depth', type=int, help="Mininum depth to filter", required=False, default=150)
    parser.add_argument('-ec', '--edgecut', type=int, nargs='+', help="Bounds to cut huge repetitive insertions", required=False, default=[10000, 270000])
    parser.add_argument('-lc', '--lencut', type=int, help="Maximum variant lenght outside of the edge cut bounds", required=False, default=2000)
    parser.add_argument('-s', '--similarity', type=float, help="Minimum similarity required to group variants", required=False, default=1.)
    parser.add_argument('-sa', '--sample', type=int, help="Sample of interest, only one (1 to 10)", required=True)
    parser.add_argument('-it', '--iterations', type=str, help="Iterations to analyse, range (1=P15 to 5=P90 => 1-5)", required=False, default="1-5")
    parser.add_argument('-o', '--output', type=str, help="Directory to write the output", required=True)
    args = parser.parse_args()

    iterations = [15, 30, 50, 65, 90]
    args.iterations = [iterations[i-1] for i in parse_range(args.iterations)]
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


    # Getting data
    it_labels = [f"P{i}-{args.sample}" for i in args.iterations]
    before_labels = [f"P{i}-{args.sample}" for i in args.iterations if i < 30]
    after_labels = [f"P{i}-{args.sample}" for i in args.iterations if i >= 30]

    data = [parse_vcf_noerror(build_vcf_path(args.sample, i)) for i in args.iterations]
    data = [[v for v in d if v['af'] > args.freq and max(v['depth'])> args.depth] for d in data]
    data = [[v for v in d if (args.edgecut[0] < v['pos'] < args.edgecut[1]) or v['svlen'] < args.lencut] for d in data]
 
    grouped_by_it = merge_samples(data, it_labels, sim_thresold=args.similarity)


    # Distribution of variations types
    types_counts = count_types([grp[0][1]["svtype"] for grp in grouped_by_it])
    types_by_it = count_types_by_group(grouped_by_it, it_labels)
    types_by_time = count_types_by_group(grouped_by_it, [before_labels, after_labels])

    fig, ax1 = plt.subplots()
    plots.barplot(types_counts, ax1)
    ax1.set_ylabel('Nombre de variants')
    ax1.set_xlabel('Types de variants')
    ax1.set_title(f'Nombre de variant par type dans l\'échantillon {args.sample}')
    plt.savefig(figs_dir / "types.pdf")

    fig, ax2 = plt.subplots()
    plots.grouped_barplot(types_by_it, it_labels, ax2)
    ax2.set_ylabel('Nombre de variants')
    ax2.set_xlabel('Échantillons du passage')
    ax2.set_title(f'Nombre de variants par passage dans l\'échantillon {args.sample}')
    plt.savefig(figs_dir / "types_by_it.pdf")

    fig, ax3 = plt.subplots()
    plots.grouped_barplot(types_by_time, ["Avant", "Après"], ax3)
    ax3.set_ylabel('Nombre de variants')
    ax3.set_title(f'Nombre de variants avant et après choc thermique dans l\'échantillon {args.sample}')
    plt.savefig(figs_dir / "type_by_time.pdf")

    types_dict = {
        "types": types_counts,
        "types_by_it": types_by_it,
        "types_by_time": types_by_time
    }


    # Distribution of lengths
    len_distrib = len_by_type([(grp[0][1]["svtype"], np.mean([x[1]["svlen"] for x in grp])) for grp in grouped_by_it])

    fig, ax4 = plt.subplots()
    plots.boxplot(len_distrib, False, ax4)
    ax4.set_ylabel("Longueurs observées")
    ax4.set_title(f"Longueurs des variants dans l'échantillon {args.sample} (sans outliers)")
    plt.savefig(figs_dir / "lengths.pdf")
    
    fig, ax5 = plt.subplots()
    plots.boxplot(len_distrib, True, ax5)
    ax5.set_ylabel("Longueurs observées")
    ax5.set_title(f"Longueurs des variants dans le passage {args.sample}")
    plt.savefig(figs_dir / "lengths_outliers.pdf")


    # Similarity matrix
    sims = pairwise_similarity(grouped_by_it, it_labels)

    fig, ax6 = plt.subplots()
    plots.heatmap(sims, it_labels, ax6)
    ax6.set_title(f"Proportion de variants partagés entre les passage de l'échantillon {args.sample}")
    plt.savefig(figs_dir / "sim_heatmap.pdf")


    # Variants uniques to each group
    uniques = []

    for g in grouped_by_it:
        present_in = 0
        counts = len(g)
        for v in g:
            if v[0] in before_labels:
                present_in = present_in | 1
            elif v[0] in after_labels:
                present_in = present_in | 2
            else:
                raise ValueError("Unexpected label")
        
        if present_in == 1 or present_in == 2:
            to_add = g[0][1].copy()
            to_add["n_found"] = counts
            to_add["samples_found"] = [v[0] for v in g]
            to_add["depth"] = [v[1]["depth"] for v in g]
            to_add["id"] = [v[1]["id"] for v in g]
            to_add["choc"] = "before" if present_in == 1 else "after"
            if "alt" in to_add:
                to_add["alt"] = str(to_add["alt"][0])
            uniques.append(to_add)

    
    # Uniques in ORFS
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
    
    # Write final results
    tl = [[x if not np.isnan(x) else None for x in row] for row in sims]

    results = {
        "types": types_dict,
        "lengths": len_distrib,
        "pairwise_sim": tl,
    }

    results["iterations"] = args.iterations
    results["sample"] = args.sample
    results["similarity"] = args.similarity
    results["filters"] = {
        "freq": args.freq,
        "depth": args.depth,
        "edgecut": args.edgecut,
        "lencut": args.lencut
    }

    with open(args.output / "results.json", 'w') as fp:
        json.dump(results, fp, indent=4)

    with open(args.output / "uniques.json", 'w') as fp:
        json.dump(uniques, fp, indent=4)
