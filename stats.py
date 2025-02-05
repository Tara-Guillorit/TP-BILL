from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

## Collect data ##

def count_types(types):
    """ Compte le nombre de variant structurel par type

    Args:
        types (list): liste des types de chaque variation
    
    Return:
        dict: dictionaire avec les type pour clés et leur nombre d'occurences pour valeur
    """
    types_count = defaultdict(int)
    for t in types:
        types_count[t] += 1
    return types_count

def len_by_type(lengths):
    """ Récupère la distribution de longueur (en valeur absolue) des variants par type 

    Args: 
        lengths (list): liste de tuple sous la forme (type de variation, longueur)
    
    Return:
        dict: dictionaire avec les type pour clés et les distributions de taille pour valeur 
    """
    types_lengths = defaultdict(list)
    for t, l in lengths:
        types_lengths[t].append(abs(l))
    return types_lengths


## Plot data ##

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
        grouped_counts (dict): dictionaire avec les groupes pour clés et une liste de nombre d'occurence pour valeur
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