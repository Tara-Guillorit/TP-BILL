from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

from .similarity import contain_from_sample


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


def pairwise_similarity(grouped_sv, sample_labels):
    """ Calcule le pourcentage de variations partagé entre les échantillons représentés dans <sample_ids>

    Args:
        grouped_sv (list): liste de groupes de variations, chaque groupe est une liste de tuples (sample id, variation)
        sample_labels (list): labels des échantillons représentés dans <grouped_sv>

    Return:
        np.array : matrice carré de la taille de <sample_ids> avec la proportion de variation partagé par pair d'échantillons
    """
    sims = np.zeros(shape=(len(sample_labels), len(sample_labels)), dtype=np.float64)
    np.fill_diagonal(sims, np.nan)
    for i in range(len(sample_labels)):
        #sims[i][i] = 1
        for j in range(0, i):
            from_i_or_j = [g for g in grouped_sv if contain_from_sample(sample_labels[i], g) or contain_from_sample(sample_labels[j], g)]
            from_i_and_j = [g for g in grouped_sv if contain_from_sample(sample_labels[i], g) and contain_from_sample(sample_labels[j], g)]
            sims[i][j] = len(from_i_and_j) / len(from_i_or_j) if len(from_i_or_j) > 0 else 0
            sims[j][i] = sims[i][j]
    return sims
