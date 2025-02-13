from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

from .similarity import contain_from_sample


def count_types(types):
    types_count = defaultdict(int)
    for t in types:
        types_count[t] += 1
    return types_count


def count_types_by_sample(grouped_sv, labels):
    result = {}
    for g in grouped_sv:
        t = g[0][1]['svtype']
        if t not in result:
            result[t] = [0] * len(labels)
        for sv in g:
            result[t][labels.index(sv[0])] += 1
    return result


def count_types_by_group(grouped_sv, grouped_labels):
    result = {}
    for g in grouped_sv:
        t = g[0][1]['svtype']
        if t not in result:
            result[t] = [0] * len(grouped_labels)

        found = [0] * len(grouped_labels)
        for sv in g:
            for i in range(len(grouped_labels)):
                if sv[0] in grouped_labels[i]:
                    found[i] = 1

        result[t] = [result[t][i] + found[i] for i in range(len(result[t]))]

    return result


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


def unique_variants(grouped_variants, group_1, group_2, label_1, label_2):
    uniques = []
    
    for g in grouped_variants:
        present_in = 0
        counts = len(g)
        for v in g:
            if v[0] in group_1:
                present_in = present_in | 1
            elif v[0] in group_2:
                present_in = present_in | 2
            else:
                raise ValueError("Unexpected label")
        
        if present_in == 1 or present_in == 2:
            to_add = g[0][1].copy()
            to_add["n_found"] = counts
            to_add["samples_found"] = [v[0] for v in g]
            to_add["ids"] = [v[1]["id"] for v in g]
            to_add["group"] = label_1 if present_in == 1 else label_2
            uniques.append(to_add)
    
    return uniques
        