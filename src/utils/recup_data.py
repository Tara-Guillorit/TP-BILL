from misc import build_vcf_path
from read_vcf import parse_vcf_noerror
from similarity import merge_samples_nolabel
import pandas as pd
import sys

samples = list(range(1, 11))
iterations = [15, 30, 50 , 65, 90]

all_var = []
for s in samples:
    for i in iterations:
        print("parsing passage", i, " of sample", s)
        data = parse_vcf_noerror(build_vcf_path(s, i))
        for d in data:
            d["sample"] = s
            d["iteration"] = i
        all_var += data
    
print("grouping variants")
grouped = merge_samples_nolabel(all_var, 1)

flat_grouped = []
for i, g in enumerate(grouped):
    for v in g:
        v["group"] = i
        flat_grouped.append(v)

print("converting data to csv in ", sys.argv[1])
df = pd.DataFrame(flat_grouped)
df.to_csv(sys.argv[1])


