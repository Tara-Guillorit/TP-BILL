def variant_equal(v1, v2, similarity=0.5):
    """ Détermine si la proportion de sites commun entre v1 et v2
    est superieur à <similarity> (basé sur la plus courte séquence)

    Args:
        v1, v2 (dict) : dictionnaire avec au minimum les clés "start", "end" et "svtype"
        similarity (float) : proportion de similarité entre 0 et 1
    
    
    """
    overlap_start = max(v1["start"], v2["start"])
    overlap_end = min(v1["end"], v2["end"])
    overlap_len = max(overlap_end - overlap_start, 0)

    len_v1 = v1["end"] - v1["start"]
    len_v2 = v2["end"] - v2["start"]
    shortest_v = min(len_v1, len_v2)

    print("sim : ", (overlap_len / shortest_v))
    return overlap_len >= similarity * shortest_v and v1["svtype"] == v2["svtype"]


def group_variants(samples, min_sim=0.5):
    # merge toute les sv dans une liste de tuples (échantillon d'origin, variant)
    sv_total = [(s, v) for v in samples[s] for s in range(len(samples))]
    # tri les sv par position de départ
    sv_total = sorted(sv_total, lambda x: x[1]["start"])

    # groupe les sv similaire ensemble
    sv_grouped = []
    while len(sv_total) > 0: # tant qu'il reste des éléments à grouper
        v1 = sv_total[0]
        del sv_total[0]
        group = [v1]

        j = 0
        while j < len(sv_total):
            v2 = sv_total[j]

            if v2[1]["end"] < v1[1]["start"]: # v2 finit avant le début de v1 (pas encore d'intersection possible)
                j += 1
                pass
            if v2[1]["start"] > v1[1]["end"]: # v2 et toutes les prochaines sv commence apres la fin de v1 (plus aucune intersection possible)
                break

            if variant_equal(v1[1], v2[1], min_sim): # add to group, delete from tab, do not increment
                group.append(v2)
                del sv_total[j]
            else:
                j += 1

        sv_grouped.append(group)

    return sv_grouped
