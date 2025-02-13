import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np
from pathlib import Path


def barplot(counts, ax):
    """ Affiche un barplot avec les valeur exactes en haut des barres
    """
    ax.bar(list(counts.keys()), list(counts.values()))

    # display exact value on top
    ax.set_ylim(0, ax.get_ylim()[1] * 1.1)
    xlocs = ax.get_xticks()
    for i, v in enumerate(list(counts.values())):
        ax.text(xlocs[i] - 0.05 * (len(str(v)) / 2), v + 0.02 * ax.get_ylim()[1], str(v))


def grouped_barplot(grouped_counts, labels, ax):
    """ Affiche un barplot avec chaque bar divisé en groupes
    """
    bottom = np.zeros(len(labels))

    for group, group_count in grouped_counts.items():
        p = ax.bar(labels, group_count, label=group, bottom=bottom)
        bottom += group_count
        ax.bar_label(p, label_type='center')
    
    ax.legend()


def boxplot(distribs, showfliers, ax):
    """ Affiche un boxplot avec une box pour chaque entré de distribs (dictionnaire {label : distrib})
    """
    ax.boxplot(list(distribs.values()), showfliers=showfliers)
    ax.set_xticks([y + 1 for y in range(len(distribs))], labels=list(distribs.keys()))
    ax.yaxis.grid(True)


def heatmap(sim, sample_labels, ax):
    """ Affiche une heatmap représentant la matrice sim, et affiche les valeurs arondies pour chaque case
    """
    im = ax.imshow(sim)

    ax.set_xticks(range(len(sample_labels)), labels=sample_labels, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(len(sample_labels)), labels=sample_labels)

    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("", rotation=-90, va="bottom")

    # Loop over data dimensions and create text annotations.
    for i in range(len(sample_labels)):
        for j in range(len(sample_labels)):
            text = ax.text(j, i, round(sim[i, j], 2), ha="center", va="center", color="w")
    


def plot_lines(Y, x_labels, colors, ax, plot_nulls=False):
    """ Plot each list in y as a line plot
    - y is a dict like {sample label (str) : sample frequencies (list), ...}
    - x_labels is a list of label for each index in the sample frequencies (usually the iterations : P15, ..., P90)
    - colors is a dict with color scheme as key (str) and 
    """
    x = [i + 1 for i in range(len(x_labels))]
    ax.set_xticks(x, labels=x_labels)

    colors_index = defaultdict(int)

    for sample, y in Y.items():
        if plot_nulls or (sum(y) > 0):
            color = None
            for c, v in colors.items():
                if sample in v:
                    color = plt.get_cmap(c)((colors_index[c]/ (len(v)-1)) if len(v) > 1 else 0)
                    colors_index[c] += 1
            ax.plot(x, y, color=color, label=sample)
    ax.legend()


def distrib_af_depth(af, af_filt, depth, depth_filt, it, sample, figdir):
    fig, ax1 = plt.subplots()
    ax1.hist(af, bins=50, label="Non filtré")
    ax1.hist(af_filt, bins=50, label="Filtré")
    ax1.set_xlabel('Fréquence allélique (AF)')
    ax1.set_ylabel('Nombre de variants')
    if it is not None:
        ax1.set_title(f"Fréquence alléliques du passage {it}")
    else:
        ax1.set_title(f"Fréquence alléliques de l'échantillon {sample}")
    ax1.legend()
    plt.savefig(figdir / "af.pdf")

    fig, ax2 = plt.subplots()
    ax2.hist(depth, bins=50, label="Non filtré")
    ax2.hist(depth_filt, bins=50, label="Filtré")
    ax2.set_xlabel('Profondeur de read (DR + DV)')
    ax2.set_ylabel('Nombre de variants')
    if it is not None:
        ax2.set_title(f"Profondeur de read du passage {it}")
    else:
        ax2.set_title(f"Profondeur de read de l'échantillon {sample}")
    ax2.legend()
    plt.savefig(figdir / "depth.pdf")

    fig, ax3 = plt.subplots()
    ax3.plot(depth, af, 'o', label="Non filtré")
    ax3.plot(depth_filt, af_filt, 'o', label="Filtré")
    ax3.set_yscale('log')
    ax3.set_xscale('log')
    ax3.set_xlabel('Profondeur de read (DR + DV)')
    ax3.set_ylabel('Fréquence allélique (AF)')
    if it is not None:
        ax3.set_title(f'Fréquence allélique sur Profondeur read dans le passage {it}')
    else:
        ax3.set_title(f'Fréquence allélique sur Profondeur read dans l\'échantillon {sample}')
    ax3.legend()
    plt.savefig(figdir / "depth_af.pdf")


def len_by_pos(lens, lens_filt, genome_size, steps, max_y, it, sample, figdir):
    step_size = genome_size / steps
    x = [(step_size) * i for i in range(0, steps)]
    y = [[] for i in range(0, steps)]
    y_fil = [[] for i in range(0, steps)]

    for l in lens:
        step = int(l[0] / (step_size))
        y[step].append(l[1])
    y = [np.mean(y[i]) if len(y[i]) > 0 else 0 for i in range(len(y))]

    for l in lens_filt:
        step = int(l[0] / (step_size))
        y_fil[step].append(l[1])
    y_fil = [np.mean(y_fil[i]) if len(y_fil[i]) > 0 else 0 for i in range(len(y_fil))]

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))

    ax1.bar(x[0:int(0.2*steps)], y[0:int(0.2*steps)], width=3000, label="Non filtré")
    ax1.bar(x[0:int(0.2*steps)], y_fil[0:int(0.2*steps)], width=3000, label="Filtré")
    ax1.set_ylim(ax1.get_ylim()[0], max_y)
    ax1.set_xlabel(f'Position de début de la mutation')
    ax1.set_title(f"Début du génome")
    ax1.legend()

    ax2.bar(x[int(0.8*steps):steps], y[int(0.8*steps):steps], width=3000, label="Non filtré")
    ax2.bar(x[int(0.8*steps):steps], y_fil[int(0.8*steps):steps], width=3000, label="Filtré")
    ax2.set_ylim(ax2.get_ylim()[0], max_y)
    ax2.set_xlabel(f'Position de début de la mutation')
    ax2.set_title(f"Fin du génome")
    ax2.legend()

    if it is not None:
        plt.suptitle(f"Taille moyenne des variants aux bornes du génome dans le passage {it}")
    else:
        plt.suptitle(f"Taille moyenne des variants aux bornes du génome dans l'échantillon {sample}")

    plt.tight_layout()
    plt.savefig(figdir / "len_by_pos.pdf")

    len_by_step = [(x[i], y[i]) for i in range(len(x))]
    len_by_step_filter = [(x[i], y_fil[i]) for i in range(len(x))]
    return len_by_step, len_by_step_filter
