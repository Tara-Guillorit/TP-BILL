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

from read_vcf import parse_vcf #on importe la fontion pour cree la liste de dico
from read_ORF import list_interval_with_dico



#On choisit une valuer minimale pour la freqeunce allelique AF:
doc_vcf = sys.argv[1]
seuil_de_AF = sys.argv[2]

dico_vcf = parse_vcf(doc_vcf) #la variable contien la liste de dico

#cette fonction retoune une liste avec les dico qui on une AF >= au seuil
def list_vcf_with_dico (seuil):
    list_vfc_with_seuil = []
    for ligne in dico_vcf:
        if ligne['af']>= float(seuil):
            list_vfc_with_seuil.append(ligne)
    return list_vfc_with_seuil

#print (list_vcf_with_dico(seuil_de_AF))




#on verifie si non mutation sont dans un ORF avec un double boucle for et retourne une liste des dico vfc et ORF (mis dans une liste) qui respecte cette condition
#ex list_pos_in_interval = [[{dico_vcf},{dico_ORF}],[{dico_vcf},{dico_ORF}],[{dico_vcf},{dico_ORF}],[{dico_vcf},{dico_ORF}],....]
#ex dico de vcf :  { 'pos': , 'id': , 'svtype': , 'svlen': , 'end': , 'af': }
#ex de dico de ORF : { 'locus_tag': , 'protein_id' : , 'location' : [[start, end],[start,end]] , 'direct' : True or False , 'complement' : True or False , 'join' : True or False }
def list_pos_in_interval_with_dico (list_vcf, list_intervals):
    list_pos_in_interval = []
    for line_vcf in list_vcf:
        for line_interval in list_intervals:
            #on doit considerer qui peut avoir plusieur intervales pour un même ORF
            for interval in  line_interval['location'] :
                if ((line_vcf['pos']>=int(interval[0]) and line_vcf['pos']<=int(interval[1])) or(line_vcf['end']>=int(interval[0]) and line_vcf['end']<=int(interval[1]))): #chancher le if pour que on considere toute la longueur de la mutation , im peut avoir des indel qui sont dans plusieur ORF
                    list_pos_in_interval.append([line_vcf,line_interval])
    return list_pos_in_interval

#print(list_pos_in_interval_with_dico(list_vcf_with_dico(seuil_de_AF),list_interval_with_dico()))

#dans cette fonction on retourne que une liste de dico, et ke dico aura les info suivantes:
#list_pos_in_interval2 = [ { 'pos': , 'id': , 'svtype': , 'svlen': , 'end': , 'af': , 'locus_tag': , 'location' : [[start, end],[start,end]] }]
def list_pos_in_interval_with_dico_2 (list_vcf, list_intervals):
    list_pos_in_interval2 = []
    for line_vcf in list_vcf:
        for line_interval in list_intervals:
            #on doit considerer qui peut avoir plusieur intervales pour un même ORF
            for interval in  line_interval['location'] :
                if ((line_vcf['pos']>=int(interval[0]) and line_vcf['pos']<=int(interval[1])) or(line_vcf['end']>=int(interval[0]) and line_vcf['end']<=int(interval[1]))): #chancher le if pour que on considere toute la longueur de la mutation , im peut avoir des indel qui sont dans plusieur ORF
                    line_vcf['locus_tag'] = line_interval['locus_tag']
                    line_vcf['location'] = line_interval['location']
                    list_pos_in_interval2.append(line_vcf)
    return list_pos_in_interval2

#print(list_pos_in_interval_with_dico_2(list_vcf_with_dico(seuil_de_AF),list_interval_with_dico()))


# fonction qui extrait les informations essentielle pour pouvoir faire une representation graphique
def extract_info ():
    list_all = list_pos_in_interval_with_dico(list_vcf_with_dico(seuil_de_AF),list_interval_with_dico())
    list_filter = []
    for line in list_all:
        # on extrait la position de la mutation : line[0]['pos']
        #on extrait la postition de fin de mutation : line[0]['end']
        # on extrait le type de la mutation : line[0]['svtype']
        # on extrait la longueur de la mutation : line[0]['svlen']
        # on extrait le nom de l'ORF : line[1]['locus_tag']
        # on extrait la ou les intervales de l'ORF : line[1]['location']
        list_filter.append([ line[0]['pos'],line[0]['end'],line[0]['svtype'],line[0]['svlen'],line[1]['locus_tag'],line[1]['location']])
    return list_filter

#print(extract_info())


# fonction qui extrait les informations essentielle pour pouvoir faire une representation graphique
def extract_info2 ():
    list_all = list_pos_in_interval_with_dico(list_vcf_with_dico(seuil_de_AF),list_interval_with_dico())
    list_filter = []
    for line in list_all:
        # on extrait la position de la mutation : line[0]['pos']
        #on extrait la postition de fin de mutation : line[0]['end']
        # on extrait le type de la mutation : line[0]['svtype']
        # on extrait la longueur de la mutation : line[0]['svlen']
        # on extrait le nom de l'ORF : line[1]['locus_tag']
        # on extrait la ou les intervales de l'ORF : line[1]['location']
        list_filter.append([ line[0]['pos'],line[0]['svtype'],line[0]['svlen'],line[1]['location']])
    return list_filter

#print(extract_info2())

def clean_extract ():
    doc = ""
    list = extract_info2()
    for line in list :
        for info in line :
            doc+=str(info)+" "
        doc+= "\n"
    return doc

#print(clean_extract())

#creation d'un doc text avec les info de clean extract
with open("extract.txt", "w", encoding="utf-8") as fichier:
    fichier.write(str(clean_extract()))

#representation graphique
#initialisation du graph
plt.figure(figsize=(40, 6))

# Parcourir les données et tracer les ORF et les mutations
list_for_graph = extract_info()


for idx, (pos,end, mut_type, length, name_mut , orf_interval) in enumerate(list_for_graph):
    # Tracer l'ORF (ligne horizontale) pour chaque position il peut avoir 1 ou plusieur ORF
    for ORF in orf_interval:
        plt.hlines(y=idx, xmin=int(ORF[0]), xmax=int(ORF[1]), colors='blue', label='ORF' if idx == 0 else "")

    # Tracer le point de début de la mutation
    color = 'green' if mut_type == "INS" else 'red'
    plt.scatter(pos, idx, color=color, label=f'Début {mut_type}' if idx == 0 else "")
    plt.scatter(end, idx, color=color, label=f'Fin {mut_type}' if idx == 0 else "")


    #ecrire la longueur de la mutation
    plt.annotate(f"{mut_type} {length} {name_mut}", (pos, idx), textcoords="offset points", xytext=(0, 5), ha='center')

# Ajouter des labels et une légende
plt.xlabel("Position sur le génome")
plt.ylabel("Index des ORF/mutations")
plt.title("Représentation des ORF et des mutations filtré selon le AF de : "+seuil_de_AF)
plt.legend()

# Afficher le graphique
plt.tight_layout()
plt.show()


