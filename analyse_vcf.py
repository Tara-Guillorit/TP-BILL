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


#on va extraire les position de debut et de fin de ORFS pour pouvoir verifier si nos mutations sont incluse dans ces ORF
#on adapte ce code pour utiliser des listes de dico aussi qui contien les valeur:
# { 'locus_tag': , 'protein_id' : , 'location' : [[start, end],[start,end]] , 'direct' : True or False , 'complement' : True or False , 'join' : True or False }
def list_interval_with_dico ():
    list_ORF = []
    for ORF in SeqIO.parse("ORF.fasta", "fasta"):
        header = ORF.description
        #exemple d'un header :
        # header = lcl|NC_009127.1_cds_YP_001096040.1_1 [locus_tag=CyHV3_ORF1_1] [db_xref=GeneID:11266495] [protein=protein ORF1] [protein_id=YP_001096040.1] [location=426..1199] [gbkey=CDS]

        dico = {}
        locus_tag = header.split("[")[1][:-2].split("=")[1]

        protein_id = header.split("[")[4][:-2].split("=")[1]

        location = header.split("[")[5][:-2].split("(")[-1].split("=")[-1].split(")")[0].split(",") #je decide de conserver les diferents intervales pour les join
        #que peut contenir un seul ORF (ex: location = ['55128..56513', '83162..83365', '83445..84221'])
        location_list = [] #contiendra une liste de 1 ou plusieurs intevales (ex [['264821', '265162']] ou [['265486', '265825'], ['265938', '270469']])
        for intervales in location:
            location_list.append(intervales.split(".."))

        #si la partie [location = ...] du header contien le mot "complement" alors direct = False sinon True
        #si la partie [location = ...] du header contien le mot "complement" alors complement = True sinon False
        if 'complement' in header.split("[")[5]:
            direct = False
            complement = True
        else :
            direct = True
            complement = False

        # si la partie [location = ...] du header contien le mot "join" alors join = True sinon False
        if 'join' in header.split("[")[5]:
            join = True
        else:
            join = False

        #on ajoute tout ca a au dico qui sera ensuite ajouter a la liste
        dico.update({'locus_tag':locus_tag , 'protein_id':protein_id , 'location':location_list, 'direct':direct, 'complement':complement,'join':join})
        list_ORF.append(dico)


    return list_ORF

#print(list_interval_with_dico())


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
                if (line_vcf['pos']>=int(interval[0]) and line_vcf['pos']<=int(interval[1])):
                    list_pos_in_interval.append([line_vcf,line_interval])
    return list_pos_in_interval

#print(list_pos_in_interval_with_dico(list_vcf_with_dico(seuil_de_AF),list_interval_with_dico()))


# fonction qui extrait les informations essentielle pour pouvoir faire une representation graphique
def extract_info ():
    list_all = list_pos_in_interval_with_dico(list_vcf_with_dico(seuil_de_AF),list_interval_with_dico())
    list_filter = []
    for line in list_all:
        # on extrait la position de la mutation : line[0]['pos']
        # on extrait le type de la mutation : line[0]['svtype']
        # on extrait la longueur de la mutation : line[0]['svlen']
        # on extrait le nom de l'ORF : line[1]['locus_tag']
        # on extrait la ou les intervales de l'ORF : line[1]['location']
        list_filter.append([ line[0]['pos'],line[0]['svtype'],line[0]['svlen'],line[1]['locus_tag'],line[1]['location']])
    return list_filter

print(extract_info())



#representation graphique
#initialisation du graph
plt.figure(figsize=(40, 6))

# Parcourir les données et tracer les ORF et les mutations
list_for_graph = extract_info()


for idx, (pos, mut_type, length, name_mut , orf_interval) in enumerate(list_for_graph):
    # Tracer l'ORF (ligne horizontale) pour chaque position il peut avoir 1 ou plusieur ORF
    for ORF in orf_interval:
        plt.hlines(y=idx, xmin=int(ORF[0]), xmax=int(ORF[1]), colors='blue', label='ORF' if idx == 0 else "")

    # Tracer le point de début de la mutation
    color = 'green' if mut_type == "INS" else 'red'
    plt.scatter(pos, idx, color=color, label=f'Début {mut_type}' if idx == 0 else "")

    #ecrire la longueur de la mutation
    plt.annotate(f"{mut_type} {length} {name_mut}", (pos, idx), textcoords="offset points", xytext=(0, 5), ha='center')

# Ajouter des labels et une légende
plt.xlabel("Position sur le génome")
plt.ylabel("Index des ORF/mutations")
plt.title("Représentation des ORF et des mutations")
plt.legend()

# Afficher le graphique
plt.tight_layout()
plt.show()


