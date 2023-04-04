##     _____  ____  _    _ _____   _____  ____  _____  
##   / ____|/ __ \| |  | |  __ \ / ____|/ __ \|  __ \ 
##  | (___ | |  | | |  | | |__) | (___ | |  | | |__) |
##   \___ \| |  | | |  | |  _  / \___ \| |  | |  ___/ 
##   ____) | |__| | |__| | | \ \ ____) | |__| | |     
##  |_____/ \____/ \____/|_|  \_\_____/ \____/|_|     

## Alex Holehouse (Pappu Lab and Holehouse Lab) and Jared Lalmansing (Pappu lab)
## Simulation analysis package
## Copyright 2014 - 2022
##


"""
ssdata contains all specific data used, and is disinct from configs, which are 
package specific global configurations, while ssdata contains information that
is unambigious as well-defined.
"""



THREE_TO_ONE = {'ALA':'A', 
                'CYS':'C',
                'ASP':'D',
                'GLU':'E',
                'PHE':'F',
                'GLY':'G',
                'HIS':'H', 
                'ILE':'I',
                'LYS':'K',
                'LEU':'L',
                'MET':'M',
                'ASN':'N',
                'PRO':'P',
                'GLN':'Q',
                'ARG':'R',
                'SER':'S',
                'THR':'T',
                'VAL':'V',
                'TRP':'W',
                'TYR':'Y',
                'ACE':'<',
                'NME':'>',
                'NAC':'>'}

ONE_TO_THREE = {'A':'ALA', 
                'C':'CYS',
                'D':'ASP',
                'E':'GLU',
                'F':'PHE',
                'G':'GLY',
                'H':'HIS', 
                'I':'ILE',
                'K':'LYS',
                'L':'LEU',
                'M':'MET',
                'N':'ASN',
                'P':'PRO',
                'Q':'GLN',
                'R':'ARG',
                'S':'SER',
                'T':'THR',
                'V':'VAL',
                'W':'TRP',
                'Y':'TYR',
                '<':'ACE',
                '>':'NME'}

DEFAULT_SIDECHAIN_VECTOR_ATOMS = {'ALA': 'CB', 
                                  'CYS': 'SG',
                                  'ASP': 'CG',
                                  'ASH': 'CG',
                                  'GLU': 'CD',
                                  'GLH': 'CD',
                                  'PHE': 'CZ',
                                  'GLY': 'ERROR',
                                  'HIS': 'NE2',
                                  'HID': 'NE2',
                                  'HIE': 'NE2',
                                  'HIP': 'NE2', 
                                  'ILE': 'CD1',
                                  'LYS': 'NZ',
                                  'LYD': 'NZ',
                                  'KAC': 'NZ',
                                  'KM1': 'NZ',
                                  'KM2': 'NZ',
                                  'KM3': 'NZ',
                                  'LEU': 'CG',
                                  'MET': 'CE',
                                  'ASN': 'CG',
                                  'PRO': 'CG',
                                  'GLN': 'CD',
                                  'ARG': 'CZ',
                                  'SER': 'OG',
                                  'SEP': 'OG',
                                  'THR': 'CB',
                                  'TPO': 'CB',
                                  'VAL': 'CB',
                                  'TRP': 'CG',
                                  'TYR': 'CZ',
                                  'PTR': 'CZ',
                                  'ACE': 'ERROR',
                                  'NME': 'ERROR'}

# list of valid residue names as supported by SOURSOP. These are the resnames compared againt when assessing for if a molecule is a protein or not. 
ALL_VALID_RESIDUE_NAMES = ['ALA','CYS','ASP','ASH','GLU','GLH','PHE','GLY','HIE','HIS','HID','HIP','ILE','LEU', 'LYS','LYD','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR','AIB', 'ABA','NVA','NLE', 'ORN', 'DAB','PTR','TPO','SEP', 'KAC', 'KM1', 'KM2' 'KM3', 'ACE','NME', 'FOR', 'NH2']
