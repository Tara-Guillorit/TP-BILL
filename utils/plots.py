import matplotlib.pyplot as plt
from collections import defaultdict


def plot_count_by_type(counts, file=None):
    """ Construit un barplot du nombre de variants par type

    Args: 
        counts (dict): dictionaire avec les type pour clés et leur nombre d'occurences pour valeur
        file (str): path du fichier ou sauvegarder le plot, l'affiche si None

    Return:
        None : sauvegarde le plot dans <file> ou l'affiche
    """
    fig, ax = plt.subplots()
    ax.bar(list(counts.keys()), list(counts.values()))

    ax.set_ylabel('Nombre de variants')
    ax.set_title('Nombre de variants par type')

    # display exact value on top
    ax.set_ylim(0, ax.get_ylim()[1] * 1.1)
    xlocs, xlabs = plt.xticks()
    for i, v in enumerate(list(counts.values())):
        plt.text(xlocs[i] - 0.05 * (len(str(v)) / 2), v + 0.02 * ax.get_ylim()[1], str(v))

    if file:
        plt.savefig(file)
    else:
        plt.show()


def plot_grouped_count_by_type(grouped_counts, labels, file=None):
    """ Construit un barplot du nombre de variants par type, chaque bar est segmenté (chaud / froid par exemple)

    Args: 
        grouped_counts (dict): dictionaire avec les groupes (pour chaque couleur) pour clés et une liste de nombre (pour chaque barre) d'occurence pour valeur
        labels (list): liste de labels correspondants à chaque indice des listes d'occurences
        file (str): path du fichier ou sauvegarder le plot, l'affiche si None

    Return:
        None : sauvegarde le plot dans <file> ou l'affiche
    """
    fig, ax = plt.subplots()
    bottom = np.zeros(len(labels))

    for group, group_count in grouped_counts.items():
        p = ax.bar(labels, group_count, label=group, bottom=bottom)
        bottom += group_count
        ax.bar_label(p, label_type='center')

    ax.set_ylabel('Nombre de variants')
    ax.set_title('Nombre de variants par type')
    ax.legend()

    if file:
        plt.savefig(file)
    else:
        plt.show()


def plot_len_by_type(lengths, showfliers, file=None):
    """ Construit un des boxplots avec les distributions de tailles par type

    Args: 
        counts (dict): dictionaire avec les type pour clés et leur distributions de taille pour valeur 
        showfliers (boolean): si vrai, affiche les outliers (par défaut ceux qui dépasse 1.5 * Q1 ou Q3)
        file (str): path du fichier ou sauvegarder le plot, l'affiche si None

    Return:
        None : sauvegarde le plot dans <file> ou l'affiche
    """
    fig, ax = plt.subplots()

    ax.boxplot(list(lengths.values()), tick_labels=list(lengths.keys()), showfliers=showfliers)
    ax.set_xticks([y + 1 for y in range(len(lengths))])
    ax.set_ylabel('Tailles observées')
    ax.yaxis.grid(True)
    ax.set_title("Distribution de taille des variants")

    if file:
        plt.savefig(file)
    else:
        plt.show()


def variant_heatmap(pairwise_sim, sample_labels, file=None):
    """ Construit une heatmap à partir d'une matrice de similarité entre échantillons

    Args:
        pairwise_sim (np.array): matrice carré de similarité entre les échantillons
        sample_lables (list): liste de labels correspondant aux échantillons de <sample_ids>
        file (str): path du fichier ou sauvegarder la figure, si None, affiche la figure directement avec pyplot

    Return:
        None : sauvegarde la heatmap dans file ou l'affiche directement
    """
    fig, ax = plt.subplots()
    im = ax.imshow(pairwise_sim)

    ax.set_xticks(range(len(sample_labels)), labels=sample_labels, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(len(sample_labels)), labels=sample_labels)

    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("", rotation=-90, va="bottom")

    # Loop over data dimensions and create text annotations.
    for i in range(len(sample_labels)):
        for j in range(len(sample_labels)):
            text = ax.text(j, i, round(pairwise_sim[i, j], 2), ha="center", va="center", color="w")
    
    ax.set_title("Pairwise similarity between samples")
    fig.tight_layout()
    if file:
        plt.savefig(file)
    else:
        plt.show()


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
