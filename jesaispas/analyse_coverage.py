
import vcf #package pour extraire info des vcf
from Bio import SeqIO #package pour extraire info du FASTA contenant les ORF
import matplotlib.pyplot as plt #pour faire representation graphique
import sys #pour pouvoir metre des parametre lors de l'execution de la fontion

from read_vcf import parse_vcf #on importe la fontion pour cree la liste de dico

from read_ORF import list_interval_with_dico #on importe la fontion pour cree la liste de dico


#On choisit une valuer minimale pour la freqeunce allelique AF:
doc_vcf = sys.argv[1]
seuil_de_cv = sys.argv[2]

dico_vcf = parse_vcf(doc_vcf) #la variable contien la liste de dico

#cette fonction retoune une liste avec les dico qui on une cv max >= au seuil
def list_vcf_with_dico (seuil):
    list_vfc_with_seuil = []
    for ligne in dico_vcf:
        if max(ligne['depth'])>= float(seuil):
            list_vfc_with_seuil.append(ligne)
    return list_vfc_with_seuil

#print (list_vcf_with_dico(seuil_de_cv))

list_ORF = list_interval_with_dico()

#print(list_ORF)

#on verifie si non mutation sont dans un ORF avec un double boucle for et retourne une liste des dico vfc et ORF (mis dans une liste) qui respecte cette condition
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

#print(list_pos_in_interval_with_dico_2(list_vcf_with_dico(seuil_de_cv),list_ORF))


# fonction qui extrait les informations essentielle pour pouvoir faire une representation graphique
def extract_info_cv ():
    list_all = list_pos_in_interval_with_dico_2(list_vcf_with_dico(seuil_de_cv),list_ORF)
    list_filter = []
    for line in list_all:
        # on extrait la position de la mutation : line['pos']
        #on extrait la postition de fin de mutation : line['end']
        # on extrait le type de la mutation : line['svtype']
        # on extrait la longueur de la mutation : line['svlen']
        # on extrait le nom de l'ORF : line['locus_tag']
        # on extrait la ou les intervales de l'ORF : line['location']
        list_filter.append([ line['pos'],line['end'],line['svtype'],line['svlen'],line['locus_tag'],line['location']])
    return list_filter

#print(extract_info())

#representation graphique
#initialisation du graph
plt.figure(figsize=(40, 6))

# Parcourir les données et tracer les ORF et les mutations
list_for_graph = extract_info_cv()


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
plt.title("Représentation des ORF et des mutations filtré selon la profondeur de : "+ seuil_de_cv)
plt.legend()

# Afficher le graphique
plt.tight_layout()
plt.show()
