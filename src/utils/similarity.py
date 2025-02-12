import numpy as np


def seq_identity(s1, s2):
    identity = 0
    for i in range(len(s1)):
        if s1[i] == s2[i]:
            identity += 1
    return identity / len(s1)


def variant_equal(v1, v2, sim_thresold=1):
    """ Détermine si deux variants sont identique à une proportion superieur ou égale à <sim_thresold>.
    Prend en compte le produit de :
        (i) la proportion de base partagé (même position) entre v1 et v2,
        (ii) le pourcentage d'identidé entre les bases partagés

    Args:
        v1, v2 (dict): dictionnaire avec au minimum les clés "pos", "svlen", "svtype" et "alt" pour les inversions
        sim_thresold (float): proportion de similarité minimum entre v1 et v2
    
    Return:
        boolean : vrai si la similarité est superieur ou égale à <sim_thresold>, faux sinon
    """
    if v1["svtype"] != v2["svtype"]:
        return False

    full_length = max(v1["pos"] + abs(v1["svlen"]), v2["pos"] + abs(v2["svlen"])) - min(v1["pos"], v2["pos"])
    common_length = min(v1["pos"] + abs(v1["svlen"]), v2["pos"] + abs(v2["svlen"])) - max(v1["pos"], v2["pos"])
    common_length = max(common_length, 0)
    shared = common_length / full_length if full_length > 0 else 0

    # Check if variation is an inversion and if sequence are actually given (not "<INS>" instead)
    if v1["svtype"] == "INS" and shared > 0 and v1["alt"] != "<INS>" and v2["alt"] != "<INS>":
        # Gather the shared part of each sequence
        first = v1 if v1["pos"] < v2["pos"] else v2
        second = v1 if v1["pos"] >= v2["pos"] else v2

        common_start = second["pos"] - first["pos"]
        seq1 = first["alt"][common_start:]
        seq2 = second["alt"]

        common_stop = min(len(seq1), len(seq2))
        seq1 = seq1[:common_stop]
        seq2 = seq2[:common_stop]

        # Compute sequence identity
        shared *= seq_identity(seq1, seq2)

    return shared >= sim_thresold


def merge_samples(samples, samples_labels, sim_thresold=1):
    """ Groupe les variants similaires ensemble

    Args:
        samples (list): une liste contenant, pour chaque échantillon (p90-1, p90-2, ...) une liste de variants
        samples_labels (list): le label correspondant à chaque liste de variants dans <samples>
        sim_thresold (float): proportion de similarité minimum pour grouper deux variants

    Return:
        list : la liste des groupes des variants, chaque groupe représenté par une liste de tuple sous la forme (échantillon d'origine, variant)
        exemple => [ [(sample_0, {variant ...}), (sample_3, {variant ...})], [(sample_x, {var ...})], ...]
    """
    # merge toute les sv dans une liste de tuples (échantillon d'origin, variant)
    sv_total = [(samples_labels[i], v) for i in range(len(samples)) for v in samples[i]]
    # tri les sv par position de départ
    sv_total = sorted(sv_total, key=lambda x: x[1]["pos"])

    # groupe les sv similaire ensemble
    sv_merged = []
    while len(sv_total) > 0: # tant qu'il reste des éléments à grouper
        v1 = sv_total[0]
        del sv_total[0]
        group = [v1]

        j = 0
        while j < len(sv_total):
            v2 = sv_total[j]

            if v2[1]["end"] < v1[1]["pos"]: # v2 finit avant le début de v1 (pas encore d'intersection possible)
                j += 1
                pass
            if v2[1]["pos"] > v1[1]["end"]: # v2 et toutes les prochaines sv commence apres la fin de v1 (plus aucune intersection possible)
                break

            if variant_equal(v1[1], v2[1], sim_thresold): # add to group, delete from tab, do not increment
                group.append(v2)
                del sv_total[j]
            else:
                j += 1

        sv_merged.append(group)

    return sv_merged


def contain_from_sample(sample_label, group):
    """ Détermine si un groupe de variations contient une variation de l'échantillon <sample_id>

    Args:
        sample_label (str): le label de l'échantillon dans les groupes
        group (list): liste de variations sous la forme de tuples (sample_id, variation)

    Return:
        boolean : Vrai si le groupe contient une variation de cet échantillon
    """
    for variant in group:
        if variant[0] == sample_label:
            return True
    return False


def find_similar_variant(variant_set, v, similarity=1):
    """ Recherche un variant similaire à v dans variant_set
    """
    result = None
    for i in range(len(variant_set)):
        if variant_equal(v, variant_set[i], similarity):
            if result is None:
                result = variant_set[i]
            else:
                print("Attention plus de un résultat similaire trouvé dans le set de variants")
                
    return result