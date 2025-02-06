
import vcf #package pour extraire info des vcf
#cette fonction prend en entrée le nom de du fichier vcf et donne en sorti une liste de dico qui contiene les valeur:
# { 'pos': , 'id': , 'svtype': , 'svlen': , 'end': , 'af': }

def parse_vcf (nom_vcf):

        list_vcf = [] #cette liste contiendra des dico
        # utilisation de vcf pour lire le fichier vcf
        vcf_reader = vcf.Reader(open(nom_vcf, 'r'))

        for record in vcf_reader:
            dico = {}
            dico.update({'pos':record.POS , 'id':record.ID , 'svtype':record.INFO['SVTYPE'], 'svlen':record.INFO['SVLEN'], 'end':record.INFO['END'],'af': record.INFO['AF'], 'depth': record.INFO['COVERAGE']})
            if dico['svtype'] == 'INS':
                dico['alt'] = record.ALT
    
            list_vcf.append(dico)
        return list_vcf
