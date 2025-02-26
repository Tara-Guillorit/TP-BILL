#il faut intaller la librairie suivante
#pip install pyvcf3
#pip intall Biopython
#pip intall matplotlib
#pip install numpy

#faire un fonction qui retourne le AF selon un seuil

#lien du package : https://github.com/dridk/PyVCF3
#https://github.com/Tara-Guillorit/TP-BILL

import vcf #package pour extraire info des vcf
from Bio import SeqIO #package pour extraire info du FASTA contenant les ORF
import matplotlib.pyplot as plt #pour faire representation graphique
import sys #pour pouvoir metre des parametre lors de l'execution de la fontion

#from analyse_vcf import seuil_de_AF
from read_vcf import parse_vcf #on importe la fontion pour cree la liste de dico
from read_ORF import list_interval_with_dico



#On choisit une valuer minimale pour la freqeunce allelique AF:
doc_vcf = sys.argv[1]
seuil_de_AF = sys.argv[2]
seuil_de_cv = sys.argv[3]

#doc_vcf = 'data-p90/P90-1.trimed1000.sv_sniffles.vcf'
#seuil_de_AF = "0.1"
#seuil_de_cv = "0"

dico_vcf = parse_vcf(doc_vcf) #la variable contien la liste de dico

#cette fonction retoune une liste avec les dico qui on une AF >= au seuil
def list_vcf_with_dico (seuil_AF, seuil_cv):
    list_vfc_with_seuil = []
    for ligne in dico_vcf:
        if ligne['af']>= float(seuil_AF) and max(ligne['depth'])>= float(seuil_cv):
            list_vfc_with_seuil.append(ligne)
    return list_vfc_with_seuil

list_vcf=list_vcf_with_dico(seuil_de_AF,seuil_de_cv)
#for dico in list_cvf :
#    print (dico)



list_ORF = list_interval_with_dico()

#for line in list_ORF:
#    print(line)

def extract_info (list_vcf , list_orf):
    info = ""
    for dico_vcf in list_vcf :
        for dico_orf in list_orf:
            if ((int(dico_orf['location'][0][0]) <= dico_vcf['pos'] <= int(dico_orf['location'][-1][-1]))
                    or (int(dico_orf['location'][0][0]) <= dico_vcf['end'] <= int(dico_orf['location'][-1][-1]))
                    or (int(dico_orf['location'][0][0]) >= dico_vcf['pos'] and int(dico_orf['location'][-1][-1]) <= dico_vcf['end'])):
                #info_not_str = dico_orf['locus_tag'],dico_orf['location'], dico_vcf['pos'], dico_vcf['end'], dico_vcf['svlen'], dico_vcf['svtype']
                info += (str(dico_orf['locus_tag'])+"\t"+str(dico_orf['location'])+
                         "\t"+str(dico_vcf['pos'])+"\t"+str(dico_vcf['end'])+"\t"+str(dico_vcf['svlen'])+"\t"+str(dico_vcf['svtype'])+"\n")
    return info

print(extract_info(list_vcf,list_ORF))

#creation d'un doc text avec les info de clean extract
with open("extract.txt", "w", encoding="utf-8") as fichier:
    fichier.write(extract_info(list_vcf, list_ORF))
