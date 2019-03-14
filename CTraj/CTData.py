##
################################################
##  ,-----.,--------.                 ,--.    ##
## '  .--./'--.  .--',--.--. ,--,--.  `--'    ##
## |  |       |  |   |  .--'' ,-.  |  ,--.    ##
## '  '--'\   |  |   |  |   \ '-'  |  |  |    ##
##  `-----'   `--'   `--'    `--`--'.-'  /    ##
##                                  '---'     ##
################################################
##
## Alex Holehouse (Pappu Lab)
## Simulation analysis package
## Copyright 2014 - 2018
##

# updated in 0.2.25 -> added NAC as a '>'
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
                                  'GLU': 'CD',
                                  'PHE': 'CZ',
                                  'GLY': 'ERROR',
                                  'HIS': 'NE2', 
                                  'ILE': 'CD1',
                                  'LYS': 'NZ',
                                  'LEU': 'CG',
                                  'MET': 'CE',
                                  'ASN': 'CG',
                                  'PRO': 'CG',
                                  'GLN': 'CD',
                                  'ARG': 'CZ',
                                  'SER': 'OG',
                                  'THR': 'CB',
                                  'VAL': 'CB',
                                  'TRP': 'CG',
                                  'TYR': 'CZ',
                                  'ACE': 'ERROR',
                                  'NME': 'ERROR'}
