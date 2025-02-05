from collections import defaultdict
import matplotlib.pyplot as plt

## Collect data ##

def count_by_type(variations):
    """ Compte le nombre de variant structurel par type

    Args:
        variations (list): liste des variations, chacune sous forme de dictionnaire
    
    Return:
        dict: dictionaire avec les type pour clés et leur nombre d'occurences pour valeur
    """
    types_count = defaultdict(int)
    for v in variations:
        types_count[v["svtype"]] += 1
    return types_count

def len_by_type(variations):
    """ Récupère la distribution de longueur (en valeur absolue) des variants par type 

    Args: 
        variations (list): liste des variations, chacune sous forme de dictionnaire
    
    Return:
        dict: dictionaire avec les type pour clés et les distributions de taille pour valeur 
    """
    types_lengths = defaultdict(list)
    for v in variations:
        types_lengths[v["svtype"]].append(abs(v["svlen"]))
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