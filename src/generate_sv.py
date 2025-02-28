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
            ATTENTION : dans les VCF, la position end indique la position qui suit la délétion (il faut enlever 1)
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
    return str(Seq(src).reverse_complement())


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
        print(f"\tAnalyse orf {orf["id"]} at positions {orf["location"]}")
        if len(orf["location"]) > 1:
            print(f"\t{colored('ATTENTION', 'red')} : the orf has multiple exons")
        gene_seq = get_ref_seq(orf["start"], orf["end"])
        try:
            if row["svtype"] == "INS":
                ref, alt = insert(gene_seq, pos, row["alt"], orf["start"])
                if orf["complement"]:
                    ref, alt = get_complement(ref), get_complement(alt)
                print("\tResult of the insertion")
                print("\tREF :", ref)
                print("\tALT :", alt)
                if ref.count("-") % 3 != 0:
                    print(f"\n\t{colored('BINGO', 'red')}: cadre de lecture décalé")
            else:
                ref, alt = delete(gene_seq, pos, end, orf["start"])
                if orf["complement"]:
                    ref, alt = get_complement(ref), get_complement(alt)
                print("\tResult of the deletion")
                print("\tREF :", ref)
                print("\tALT :", alt)
                if alt.count("-") % 3 != 0:
                    print(f"\n\t{colored('BINGO', 'red')}: cadre de lecture décalé")
            

        except Exception:
            print("\tThe variant is not included inside the ORF, searching around the variant")
            
            if row["svtype"] == "INS":
                src_seq = get_ref_seq(pos - 1, end + 2)
                ref, alt = insert(src_seq, pos, row["alt"], pos-1)
                if orf["complement"]:
                    ref, alt = get_complement(ref), get_complement(alt)
                print("\tResult of the insertion")
                print("\tREF :", ref)
                print("\tALT :", alt)
            else:
                src_seq = get_ref_seq(pos - 2, end + 2)
                ref, alt = delete(src_seq, pos, end, pos-2)
                if orf["complement"]:
                    ref, alt = get_complement(ref), get_complement(alt)
                print("\tResult of the deletion")
                print("\tREF :", ref)
                print("\tALT :", alt)

            if "TGA" in alt.replace("-", "") or "TAG" in alt.replace("-", "") or "TAA" in alt.replace("-", ""):
                print(f"\t\t{colored('BINGO', 'red')} : un codon stop est dans la séquence alt")
            if "ATG" in alt.replace("-", ""): 
                print(f"\t\t{colored('BINGO', 'red')} : un codon start est dans la séquence alt")
    
    print("\n==================================================")



#orf_tot = get_ref_seq(24962, 29147)
#orf_rf = orf_tot[:25538-24962+1] + orf_tot[25642-24962:]
#
#expect = """
#    ATGGCTTCAACAACGACGCCTTCTGCTTCAACGACACCAACAACAACACCCGCAGCGGTACCCACCAAGACAACCACCAA
#    GACATCATCGGCGACCGAAAATGGCGCAAGCACACGAGAGTTTGAGTGGGTCAGGTTGTGCGGAGGCACCGTCGAGAGGG
#    CACCGTGGACGTGCGTCTTTGTGACTTCCACCGTCGAGAGTCTCGACAGGTTCAATCGCGGCGCCTCTCTCTTCTCCGCC
#    AAGCACACCATCACCAAGAGCGAGTGCAGGCGCCTCAACGCGCTCTGGAGCGCTCTGAGCGACACCGACAAACAGAAGGC
#    GGCGTGGACGATGGCCGGCATCCTCGACAGATGCGACAATCACGCCTACGCGCTCGGCAAGGCCTCGACCATGCCCGAGG
#    TCGTCAAGAGGTGCCTGACGCTCTGCTCCGTGAGGCTGCCCCTCGCTGACGATCAGGACCTGCTCGAGCATCGCAAGTCG
#    ACGCCGCTCTGCGCCAGGGCCGTCTCGCTGCTGGCTCACCCTGCTCCACGCGACCTCTCCGCAGAGGCGCTCACCACCCT
#    CAGGACCGTCGCCTGCTTGGCCCAGCTCAACAGGCTGCTCACCGTCCTCGCACCCGTCAGACCCGACGTCTCCATGCCCG
#    CCAGTCCCTTCGCAGAAGTCACGGACCTCACCAGGGAGTGGATACGATGGGCCCTCGCGCACTACAGGGGAGGCAAGCAG
#    TGGTACGCTCACCAGATCGCCGCCTCGAACGTCGTCGTCGGGTTGGTGCTCTACGACAGAGAGGCCTGCGCCAAGCCGTC
#    GCAGCCCATCGCCTCCCCCATGCCCATCTTCATGCCGCCTTCGCCCCAGCAACAACAGGACCTGCAGGAGCTGGAGCTGG
#    AGTTGGCGCAAGCACACGAGGTCTTTACAGGTGGGTTCGAGGTGATGGACGAGCTGGTAGAGTCCGTGCGCACCTTCACC
#    TCGCAGAGCTCGGCCTACGACGAGATGGCCATGGACCTCGAGGAGCTGCTCTCGGGACCCGTGCCCGAACAAGAAGTCGA
#    CATCGACCAAGCCACACTCCGCGACACAGGCCTCGAAGCCATCCACCCATGCCAAGACCACGACTACCCGCACCACGAGT
#    TTGGGACCGTCACCGAGGCGCAGGCCGACGGACACGCCGAGCCCGTCCAGGGCTCACTCGCCACCCTCTTTGGACCCGTC
#    CAGGTCGGCACCGAGATCGCCGCCGAGACGACAGCGACTGCCTTCATCAAAGAACGAGAGGACGGCACCTACCTCTTTGC
#    CGTGCCCAGCGGCCTCATCGACGTGCACTACTGGCCCACACTCATGCAGGTCCTGCTCCAGCCCAACGTGCCCTCCACCG
#    TCTACCTCGTCCTCAAGAAAACGGACATGCCCATCCGCAGAGGCACCCTCAGAGTCAAGGCCAACGCGTCCACCACGCTA
#    CAAGCCAGAAACAGATGCGCCGAGATGAGGGTGCCCTCGGGATCGGAGGCCAGCGGCATGCTGCAGATGTGGAGGCCAGA
#    GGGACTCTACTTTGAGACGTGCGGAGCCAACAAACCCATCCTCAGGTGCCCTTACCACGCCGACGAGGAGAGCGCCGACG
#    TGCCCTTCGTCCAGATGGTGCGACTCTGCAAGTACCTCAGACACATCAACCACTGGCTCGACGCCGCCGAAGATCACTCC
#    GACCCAGAGATCGCCGCGCACGTCGGGATCAGGCGCAGACAGGGACTCATGGCCGAGATCGTCAAGGTGGCCATGAGCAC
#    CCTCCTAGAGAGGTCTGCCGTCGCCGCTATCCATCACAGCGACACGCCGACTCCTCACCCCGCCAACATGATATCATCGT
#    TCGCCAACAACACCACCGTCGTCTCGGGCGCCTTCAGCACCCCCGACCAGTTCACGCCGTACGGCACCACGCCGCAGTTT
#    ATCCACCCGCACCAGACGCTCGTCAGACAGTCGGTACCAGTGCTCCATCAGCAGCAGCAGATCCCCATCCTCGACTCTGG
#    CGACGCGCTCTCGGCCATCATCGGCCAGACGCTCATCACCGACGGCGAAGAGGCCAACAGAAACATCACGATGGTCGACG
#    TGCCCGTCATCCGCATGGAACATCTGCAGAGGCTCCAGCAGACCTTCCAGATGGTACAAGTGCCGCAAGGACACGAGGGG
#    CAGTTTGTGTGGAACGCCGTCGCACAAGGCCAGACGGTCGCCACGGGAGGAGAGTACGAGGCGGCGCAGATCAGCATCGC
#    AGACACCAACAGGTACTCGGCCGAGGTCCCAGAGTCTGCCAGAGAGATCATGCCTCCCACGCCTTACTGCCTGCACCCCG
#    TGCAGCTCGAGTACTACGATCCGAAACCCGTGGAAGGCGAGAGGGACGAGGAGCCCTACCAGACCACTACCGTCGAGGGC
#    ATCGTGCAAGGCGCGCCGCTCAAGAGGAGCCAGATGCCGCTCTACAACAAGGACCCGTCCTGCGCCTTCTTCATCCCCAT
#    GAAGGGCGTCAGGGACCCATCAGAGATCCCGCAAGAACACTACGAGGACGCGAGCCGCATGGCCTGCGAGTACGAGGAGA
#    GGCTCCGCAGCAAGACCAAGGGCAAGAAGAAGGCATTCAACACTTGCGACGCCATGGCATATGCCCCCACACACTACGCC
#    GCGCAGTCTTACGGCGACGGTCCAGACGGGCCCTTCTACCACGATCTGGGCAACGCGTTCAGGTACCACGAGAGACACTT
#    TCTGCCTCCCCTCAACAACACCTACGACACCATCAAGAGACGAGCGATGGGGCAGGAGAAGGCTCCGCCGGGTGCCAGGT
#    CCATCGACGCCCTCGTCACTTCGTTCTGGGCCACTCACCCCAACACCCGAGTCTTCTACGACGAGCTCGTCACCCTCGTC
#    CGCAATCAGAACTCGCAGCCAGAGTACCTCAAGCAGTACTACGCTCACGCCAACCTCAACCCCTCGCCAGACAGCCCGGC
#    AGCCCTCAAGCTCACCACAAACATCTCGGGAGATCCCATCAAGGACGGCCTCACGCTCTTCTTTCTGGCCGTGGAGTACT
#    TTGGGCACAGGCCGCAGACCGTCGTCTGGTGGGCCCGCACCAAGACGCTCGGAGAGCTCCTCACCGACGACTACTGGGCC
#    AGGTTCATCACCATGTGGGGATGCTGGAAGAAGATGATCGTGAGGCTGGCCGGCGGGGACATGTGCGCCAACGTTCGCAA
#    CGCGGTCAGGGCGCACAGGCACGTTGGCGTGGAAGAATTCCACAAGAAAGTCACAGACGGTCTCCAGGCAGTCACTCGCA
#    TGCACACCGGCGTCTTGCCATTCATGGGAGACAGAGTCGAGCCCCTACCACAGCAGCCCGTCCACCTCAGGCCCTACACC
#    AACAGCCACAACAGGTCCATGGCCCAGTCCCTCGTCGACGCATACAGCACCAACGGCAAGAGGACCTGCAACAACAGGGA
#    GCCCGTCGTGCCGCGCAAGCAGCACGCCGTCCCCTTCTTTCTACCAGACAAGACGCTCATGCTGGCGCTCAAGCAGGGCC
#    ACTTTATGGACGAGTACGCCTTCGTGCCCATCAGAGCCGCGCCGATGAGAGACAAGCAGAGCCTGGGATTCGACCCCAAG
#    TCTCCCTCCAACACCTGCAACCGGGAGGACAACGCCGCCACCATGAACCTGACCAACATCGTCTTCTGCAAAAACACCAA
#    CCAAGAGATCGGCTCCAACAACAACTCTAGACGCATCCTGGGAGCAGAGACGAGGGGACTGGCGCAGTGCGAGGCAGACA
#    ACCTGGAGCCCCGCATCACCAACCTGCCCTGCGAACCTCTGGGCGCCCTCGACCTCAGCCTCAGAGTCAAACAGACGCCG
#    TACCCAGACGCCGTCATCTTCAAGACCATCCAACCCTCTCCCAACCCCACATCAGCCAACGAACATAGACCCGACCCGGC
#    CACCTGCTCTTTCGTGGAGATCACAGCCTACATGGAGCCGCAAAGCACCACTCAGAGTCCTCCCATGGCCCTCTCTACAC
#    CCTCACCCCGCTACTCTCCCTCTTCTTCCAACACGTTACCCTCGCTGTCTTCTGAACCCTCACGAAAACGCCGCAGAACA
#    TGA
#""".replace("\n", "").replace("\t", "").replace(" ", "")
#
#print(orf_rf == expect)