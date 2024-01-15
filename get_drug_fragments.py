######################################
# Author(s): Nicole Zatorski         #
# Date last updated: 10/27/23        #
# Specs: molecule fragmentation code #
######################################

from rdkit import Chem
from rdkit.Chem import AllChem
from  rdkit.DataStructs import FingerprintSimilarity
from rdkit.Chem import FragmentCatalog
from rdkit import RDConfig
import os
import numpy as np

def get_molecules(filename_in):
    '''
    returns a list of rdkit molecules from drug SMILES
    '''
    out = []
    with open(filename_in) as fo:
        for line in fo:
            split_line = line[:-1].split(',')
            m1 = Chem.MolFromSmiles(split_line[1])
            out.append(m1)
    return out

def get_fps(filename_in):
    '''
    returns a list of fingerprints from drug SMILES
    '''
    out = []
    with open(filename_in) as fo:
        for line in fo:
            split_line = line[:-1].split(',')
            m1 = Chem.MolFromSmiles(split_line[1])
            fp = str(list(AllChem.GetMorganFingerprintAsBitVect(m1, radius=2, nBits = 512))).replace(' ','')[1:-1]
            out.append(fp)
    return out


def drug_fragments(m1):
    '''
    returns drug fragments for a drug 
    '''
    out_dict = {}
    fName=os.path.join(RDConfig.RDDataDir,'FunctionalGroups.txt')
    fcgen=FragmentCatalog.FragCatGenerator()
    fparams = FragmentCatalog.FragCatParams(1,6,fName)
    fcat=FragmentCatalog.FragCatalog(fparams)
    fcgen=FragmentCatalog.FragCatGenerator()
    num_fragments = fcgen.AddFragsFromMol(m1,fcat)
    for n in range(num_fragments):
        fragment = fcat.GetEntryDescription(n)
        if fragment not in out_dict:
            out_dict[fragment] = 1
        else:
            out_dict[fragment] = out_dict[fragment] + 1
    return out_dict
        
def iterates_over_drug_files_to_get_features(file_in):
    '''
    returns a dictionary of molecules (keys) and their fragments (dictonary containing frequencies)
    '''
    output = {}
    drugs = get_molecules(file_in)
    for drug in drugs:
        fragment_dictionary = drug_fragments(drug)
        output[drug] = fragment_dictionary
    return output
