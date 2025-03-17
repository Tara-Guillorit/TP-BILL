""" Script for generating variant ORF.
"""
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import ast
from termcolor import colored


def get_ref_seq(pos, end):
    """ Récupere la sous séquence de pos à end (inclues) dans le génome de ref
    """
    for record in SeqIO.parse("data/NC_009127.fa", "fasta"):
        reference = record.seq[pos-1:end]
        return str(reference)

def delete(src, pos, end, offset):
    """ Calcule la séquence src mutée par une délétion de pos à end (inclues)
        pos et end sont les position de la délétion dans le génome
            WARNING : dans les VCF, la position end indique la position qui suit la délétion (il faut enlever 1)
        offset correspond à la position de la séquence src dans le génome (sa première base)

        Retourne la séquence source et la séquence muté (avec des '-' pour l'insertion)
    """
    rel_pos, rel_end = pos - offset, end - offset
    if rel_end < 0 or rel_pos > len(src):
        raise Exception("Deletion not inside the source sequence")
    deleted = ""
    dest = list(src)
    for i in range(max(0, pos-offset), min(end-offset+1, len(dest))):
        deleted += dest[i]
        dest[i] = "-"
    return src, "".join(dest)

def insert(src, pos, ins, offset):
    """ Calcule la séquence src mutée par une insertion
        pos correspond à la position qui précède l'insertion dans le génome de référence
        offset correspond à la position de la séquence src dans le génome (sa première base)

        Retourne la séquence src avec des '-' au niveau de l'insertion et la séquence muté avec l'insertion
    """
    rel_pos = pos - offset
    if rel_pos < -1 or rel_pos > len(src):
        raise Exception("Insertion not inside the source sequence")
    dest = src[:rel_pos + 1] + ins + src[rel_pos + 1:]
    src_mut = src[:rel_pos + 1] + "-" * len(ins) + src[rel_pos + 1:]
    return src_mut, dest

def get_complement(src):
    """ Donne la séquence inverse complément sous forme string
    """
    return str(Seq(src).reverse_complement())


def search_stop_onframe(seq):
    seq_trim = seq.replace("-", "")
    for i in range(0, len(seq_trim) - 3, 3):
        if seq_trim[i:i+3] == "TGA" or seq_trim[i:i+3] == "TAA" or seq_trim[i:i+3] == "TAG":
            print(f"\t{colored(f'BINGO', 'red')} : Codon STOP found on frame at position {i}")


def search_start_stop(seq, pos):
    seq_trim = seq.replace("-", "")
    for i in range(len(seq_trim) - 2):
        cod = seq_trim[i:i+3]
        if cod in ["TGA", "TAA", "TAG"]:
            print(f"\tCodon STOP found in position {i}, frame {colored(f'(+{(i+pos) % 3})', 'red')}")
        if cod == "ATG":
            print(f"\tCodon START found in position {i}, frame {colored(f'(+{(i+pos) % 3})', 'red')}")
        
        if cod in ["TCA", "TTA", "CTA"]:
            print(f"\tCodon STOP found in position {i}, frame {colored(f'(-{(i+pos) % 3})', 'blue')}")
        if cod == "CAT":
            print(f"\tCodon START found in position {i}, frame {colored(f'(-{(i+pos) % 3})', 'blue')}")



orfs = {}
for seq in SeqIO.parse("data/ORF.fasta", "fasta"):
    header = seq.description

    locus_tag = header.split("[")[1][:-2].split("=")[1]

    locations = header.split("[")[5][:-2].split("(")[-1].split("=")[-1].split(")")[0].split(",")
    locations = [inter.split("..") for inter in locations]
    locations = [[int(p) for p in inter] for inter in locations]

    complement = 'complement' in header.split("[")[5]

    start_orf = min(locations, key=lambda x: x[0])[0]
    end_orf = min(locations, key=lambda x: x[1])[1]

    orfs[locus_tag] = {'id': locus_tag, 'location': locations, 'complement': complement, 'start': start_orf, 'end': end_orf}

candidats = pd.read_csv('results/candidates.csv', index_col="index")

for _, row in candidats.iterrows():
    print("==================================================\n")
    print(f"Checking variant {row["group"]} ({row["id"]}) from position {row["pos"]} to {row["end"]} in orfs {row["orfs"]}")
    var_orfs = [orfs[o] for o in ast.literal_eval(row["orfs"])]
    pos, end = int(row["pos"]), int(row["end"])
    end = end if row["svtype"] == "INS" else end-1

    for orf in var_orfs:
        print()
        print(f"Analyse orf {orf["id"]} at positions {orf["location"]}")
        print()
        if len(orf["location"]) > 1:
            print(f"\t{colored('WARNING', 'red')} : the orf has multiple exons")
        if not orf["complement"]:
            print(f"\tORF is on frame {colored(f'+{(orf["start"]-1) % 3}', 'red')}")
        else:
            print(f"\tORF is on frame {colored(f'-{(orf["start"]-1) % 3}', 'blue')}")

        gene_seq = get_ref_seq(orf["start"], orf["end"])
        try:
            if row["svtype"] == "INS":
                ref, alt = insert(gene_seq, pos, row["alt"], orf["start"])
            else:
                ref, alt = delete(gene_seq, pos, end, orf["start"])

            print()
            if (ref.count("-") + alt.count("-")) % 3 != 0:
                print(f"\n\t{colored('BINGO', 'red')}: orf is frameshifted !")

            if orf["complement"]:
                ref, alt = get_complement(ref), get_complement(alt)
            search_stop_onframe(alt)

            print()
            print(f"\tResult of the {row["svtype"]} :")
            print("\tREF :", ref)
            print("\tALT :", alt)
            
            

        except Exception:
            print(f"\t{colored('WARNING', 'red')} : The variant is not included inside the ORF, searching around the variant\n")
            
            if row["svtype"] == "INS":
                src_seq = get_ref_seq(pos - 1, end + 2)
                ref, alt = insert(src_seq, pos, row["alt"], pos-1)
                search_start_stop(alt, pos-2)
            else:
                src_seq = get_ref_seq(pos - 2, end + 2)
                ref, alt = delete(src_seq, pos, end, pos-2)
                search_start_stop(alt, pos-3)

            if orf["complement"]:
                ref, alt = get_complement(ref), get_complement(alt)

            print()
            print(f"\tResult of the {row["svtype"]} :")
            print("\tREF :", ref)
            print("\tALT :", alt)

    
    print("\n==================================================")


#ref = get_ref_seq(262177 ,265177)
#ins = "ACACTTCAAGAA"
#
#ref, alt = insert(ref, 265177, ins, 262177)
#
#
#stop = None
#for i in range(len(alt), 2, -1):
#    cod = alt[i-3:i]
#    if cod in ["TGA", "TAA", "TAG"]:
#        stop = i - 3
#        break
#
#alt = alt[:stop+3]
#print(f"Codon stop {alt[-3:]} found at position {262177 + stop}, on frame {(262177 - 1 + stop) % 3}")
#
#prev_stop = None
#for i in range(len(alt)-3, 2, -3):
#    cod = alt[i-3:i]
#    if cod in ["TGA", "TAA", "TAG"]:
#        prev_stop = i - 3
#        break
#
#alt = alt[prev_stop:]
#print(f"Previous codon stop {alt[:3]} found at position {262177 + prev_stop}, on frame same frame")
#
#start = None
#for i in range(0, len(alt)-2, 3):
#    cod = alt[i:i+3]
#    if cod == "ATG":
#        start = i
#        break
#
#alt = alt[start:]
#print(f"Codon start of ORF found at position {262177 + prev_stop + 3 + start}\n")
#print(alt)