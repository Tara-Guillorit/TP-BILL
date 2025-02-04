from collecitons import defaultdict

def compute_variation_types(variations):
    types_count = defaultdict(int)
    for v in variations:
        types_count[v["svtype"]] += 1
    return types_count

def compute_variation_len(variations):
    types_lengths = defaultdict(list)
    for v in variations:
        types_lengths[v["svtype"]].append()