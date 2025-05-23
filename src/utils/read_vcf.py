
import vcf #package pour extraire info des vcf
import tempfile
from contextlib import contextmanager
import re
import os
import numpy as np
#cette fonction prend en entrée le nom de du fichier vcf et donne en sorti une liste de dico qui contiene les valeur:
# { 'pos': , 'id': , 'svtype': , 'svlen': , 'end': , 'af': }


@contextmanager
def temp_file_fixed(nom_vcf):
    """ Créé un fichier temporaire avec le vcf corrigé (coverage None remplacés par 0) et retourne son chemin
    """
    temp_file = tempfile.NamedTemporaryFile(delete=False, mode='w')
    try:
        with open(nom_vcf, 'r') as infile:
            file_contents = infile.read()
        modified_contents = re.sub(r'(None)([;,])', r'0\g<2>', file_contents)

        temp_file.write(modified_contents)
        temp_file.close()

        yield temp_file.name
    finally:
        # Ensure the file is closed and deleted after use
        os.remove(temp_file.name)


def get_dr_dv(record, af):
    called = record.get_unknowns() + record.get_hom_refs() + record.get_hom_alts() + record.get_hets()
    if len(called) == 0:
        print("None found")
        return 0, 0
    elif len(called) > 1:
        print(f"{len(called)} founds")
    return called[0].data.DR, called[0].data.DV

def parse_vcf (nom_vcf):
        """ Retourne la liste des variants structurelles contenus dans <nom_vcf> sous la forme d'une liste de dictionnaire
            Chaque variation est sous la forme {'id': str, 'pos': int, 'end': int, 'svlen': int, 'svtype': str, 'depth': list[int], 'af': float}
        """
        list_vcf = [] #cette liste contiendra des dico
        # utilisation de vcf pour lire le fichier vcf
        with open(nom_vcf, 'r') as f:
            vcf_reader = vcf.Reader(f)

            for record in vcf_reader:
                dico = {}
                af = record.INFO['AF'] if 'AF' in record.INFO else 0
                dr, dv = get_dr_dv(record, af)
                
                dico.update({'pos':record.POS , 'id':record.ID , 'svtype':record.INFO['SVTYPE'], 'svlen':record.INFO['SVLEN'], 'end':record.INFO['END'],'af': af, 'dv': dv, 'dr': dr, 'depth': dr+dv})
                dico['alt'] = str(record.ALT[0])
                dico['ref'] = str(record.REF)
        
                list_vcf.append(dico)
            return list_vcf


def parse_vcf_noerror(nom_vcf):
    """ Parse un fichier vcf en gérant les erreurs causé par les coverage None
    """
    try:
        return parse_vcf(nom_vcf)
    except:
        with temp_file_fixed(nom_vcf) as tp_name:
            return parse_vcf(tp_name)

