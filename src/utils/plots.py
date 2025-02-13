import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np


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
    ax.boxplot(list(distribs.values()), tick_labels=list(distribs.keys()), showfliers=showfliers)
    ax.set_xticks([y + 1 for y in range(len(distribs))])
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
