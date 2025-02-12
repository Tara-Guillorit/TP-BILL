from Bio import SeqIO #package pour extraire info du FASTA contenant les ORF

#on va extraire les position de debut et de fin de ORFS pour pouvoir verifier si nos mutations sont incluse dans ces ORF
#on adapte ce code pour utiliser des listes de dico aussi qui contien les valeur:
# { 'locus_tag': , 'protein_id' : , 'location' : [[start, end],[start,end]] , 'direct' : True or False , 'complement' : True or False , 'join' : True or False }
def list_interval_with_dico (filename):
    list_ORF = []
    for ORF in SeqIO.parse(filename, "fasta"):
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
