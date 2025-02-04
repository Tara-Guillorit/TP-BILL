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
import sys





#On choisit une valuer minimale pour la freqeunce allelique AF:
doc_vcf = sys.argv[1]
seuil_de_AF = sys.argv[2]


def list_AF (seuil):
    list_of_info = []
    #utilisation de vcf pour lire le fichier vcf
    vcf_reader = vcf.Reader(open(doc_vcf,'r'))
    #print("AF ID POS TYPE LEN ")
    for record in vcf_reader:
        if (record.INFO['AF']>=float(seuil)) : #verifie si le AF est superieur ou egal a ce seuil
            #on ajoute a la liste list_of_data les info pertinentes
            list_of_info.append([record.INFO['AF'],record.ID,record.POS,record.INFO['SVTYPE'] ,record.INFO['SVLEN']])
    liste_triee = sorted(list_of_info, key=lambda x: x[0], reverse = True) #trier la liste en ordre croissant des AF
    return liste_triee

#pour voir le contenu de la liste list_of_info decommenter ceci :
#for i in range (len(liste_triee)):
#    print(liste_triee[i])


#on va extraire les position de debut et de fin de ORFS pour pouvoir verifier si nos mutations sont incluse dans ces ORF
#ensuite on verifas si les mutation decale de cadre de lecture
def list_interval ():
    list_of_intervals = []
    for ORF in SeqIO.parse("ORF.fasta", "fasta"):
        description = ORF.description
        #description = lcl|NC_009127.1_cds_YP_001096040.1_1 [locus_tag=CyHV3_ORF1_1] [db_xref=GeneID:11266495] [protein=protein ORF1] [protein_id=YP_001096040.1] [location=426..1199] [gbkey=CDS]

        #print (description.split("[")[4:6])
        location = description.split("[")[5][:-2].split("(")[-1].split("=")[-1].split(")")[0].split(",")

    for interval in location:
        list_of_intervals.append(interval.split(".."))

    return list_of_intervals


#for i in range (len(list_of_intervals)):
    #print(list_of_intervals[i])
    # explication des splits: à faire



#on verifie si non mutation sont dans un ORF avec un double boucle for
def list_pos_in_interval(liste_des_AF, liste_des_intervales):
    list_of_POS_and_ORF_interval = []
    print ("liste de mutation qui sont dans des ORF et l'intervale de l'ORF  : ")
    for line_AF in liste_des_AF :
    #line_AF[2] correspont a la POS
    # line_AF[3] correspont a la LEN
    #line_AF[4] correspont a la LEN
        for line_interval in liste_des_intervales:

            if (line_AF[2]>=int(line_interval[0]) and line_AF[2]<=int(line_interval[1])):
                list_of_POS_and_ORF_interval.append([line_AF[2],line_AF[3],line_AF[4], [int(line_interval[0]),int(line_interval[1] )]])
                print (line_AF[2],line_AF[3],line_AF[4], line_interval) #renvoit la position de la mutation et l'intervale ou elle a ete touver
    return list_of_POS_and_ORF_interval




# Initialisation du graphique
plt.figure(figsize=(40, 6))

# Parcourir les données et tracer les ORF et les mutations
for idx, (pos, mut_type, length, orf_interval) in enumerate(list_of_POS_and_ORF_interval):
    # Tracer l'ORF (ligne horizontale)
    plt.hlines(y=idx, xmin=orf_interval[0], xmax=orf_interval[1], colors='blue', label='ORF' if idx == 0 else "")

    # Tracer le point de début de la mutation
    color = 'green' if mut_type == "INS" else 'red'
    plt.scatter(pos, idx, color=color, label=f'Début {mut_type}' if idx == 0 else "")

    #ecrire la longueur de la mutation
    plt.annotate(f"{mut_type} {length}", (pos, idx), textcoords="offset points", xytext=(0, 5), ha='center')

# Ajouter des labels et une légende
plt.xlabel("Position sur le génome")
plt.ylabel("Index des ORF/mutations")
plt.title("Représentation des ORF et des mutations")
plt.legend()

# Afficher le graphique
plt.tight_layout()
plt.show()
