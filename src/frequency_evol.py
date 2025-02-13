import argparse
import os
from pathlib import Path
import json
import matplotlib.pyplot as plt

from utils.plots import plot_lines
from utils.similarity import find_similar_variant
from utils.misc import parse_range, build_vcf_path, build_vcf_path_test
from utils.read_vcf import parse_vcf_noerror


def read_by_sample(sample, its):
    result = []
    for i in its:
        p = build_vcf_path_test(sample, i)
        result.append(parse_vcf_noerror(p))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Frequency evolution analyse options")

    # Mutation :  {'pos': 143509, 'id': 'Sniffles2.DEL.318S0', 'svtype': 'DEL', 'svlen': -24, 'end': 143533, 'af': 0.905, 'depth': [1772.0, 1807.0, 1805.0, 1819.0, 1875.0], 'n_found': 3}
    # Commande :  python3 frequency_evol.py -p 143509 -l -24 -t 'DEL' -s 1 -sa 1-10 -it 1-5 -o results/frequency_evol/ORF78

    parser.add_argument('-p', '--pos', type=int, help="Position de départ", required=True)
    parser.add_argument('-l', '--len', type=int, help="Taille de la mutation", required=True)
    parser.add_argument('-t', '--type', type=str, help="Type de la mutation", required=True)
    parser.add_argument('-s', '--similarity', type=float, help="Proportion de similarité requise", required=False, default=1.)
    parser.add_argument('-a', '--alt', type=str, help="Séquence alternative en cas d'insertion", required=False, default="")
    parser.add_argument('-sa', '--samples', type=str, help="Samples of interest, range (1 to 10 => 1-10)", required=False, default="1-10")
    parser.add_argument('-it', '--iterations', type=str, help="Iterations to analyse, range (1=P15 to 5=P90 => 1-5)", required=False, default="1-5")
    parser.add_argument('-o', '--output', type=str, help="Directory to write the output", required=True)
    args = parser.parse_args()

    args.samples = parse_range(args.samples)
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


    # list of all variants for each iteration (value), for each sample (key)
    data = {s: read_by_sample(s, args.iterations) for s in args.samples}
    searched = {'pos': args.pos, 'svlen': args.len, 'svtype': args.type, 'alt': args.alt}

    # list of occurence of the variant of interest for each iteration (value), for each sample (key)
    variant_occs = {}
    for s in args.samples:
        sample_its = data[s]
        variant_occs[s] = [find_similar_variant(it, searched, args.similarity) for it in sample_its]

    #print("Occurences")
    #print_dict(variant_occs)

    # list of frequency of the variant of interest for each iteration (value), for each sample (key)
    variant_freqs = {}
    for s in args.samples:
        variant_freqs[s] = [occ['af'] if occ is not None else 0 for occ in variant_occs[s]]

    #print("frequences")
    #print_dict(variant_freqs)

    results = {}
    results["variants"] = {s: {
        "occurences": [occ['id'] if occ is not None else None for occ in variant_occs[s]],
        "frequencies": variant_freqs[s],
        } for s in args.samples}
    results["searched"] = searched
    results["iterations"] = args.iterations
    results["samples"] = args.samples
    results["similarity"] = args.similarity

    with open(args.output / "results.json", 'w') as fp:
        json.dump(results, fp, indent=4)


    colors = {
        'autumn': [s for s in args.samples if s > 5],
        'winter': [s for s in args.samples if s <= 5]
    }

    fig, ax = plt.subplots()
    plot_lines(variant_freqs, args.iterations, colors, ax, False)
    ax.set_xlabel("Passages de l'expérience")
    ax.set_ylabel("Fréquence allélique de la mutation")
    #plt.show()

    plt.savefig(figs_dir / "frequency_evol_by_sample.pdf")
