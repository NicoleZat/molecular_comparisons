######################################
# Author(s): Nicole Zatorski         #
# Date last updated: 1/15/24         #
# Specs: molecule fragmentation code #
######################################

from rdkit import Chem
from rdkit import DataStructs

def drug_sim(input_csv_file_name, output_csv_file_name):
    '''
    reads a csv file of comma sepparated drug names 
    writes the tanimoto, dice, cosine, and average of those three similarities to a new csv file
    '''
    drug_molecules = {}
    with open(input_csv_file_name) as fo:
        for line in fo:
            split_line = line.replace('\n','').split(',')
            drug_molecules[split_line[0]] = Chem.RDKFingerprint(Chem.MolFromSmiles(split_line[2]))

    keys = list(drug_molecules.keys())
    for key_index in range(len(keys)):
        k1 = keys[key_index]
        m1 = drug_molecules[k1]
        for key_index2 in range(key_index+1,len(keys)):
            if key_index+1 == len(keys):
                continue
            k2 = keys[key_index2]
            m2 = drug_molecules[k2]
            tanimoto = DataStructs.FingerprintSimilarity(m1,m2)
            dice = DataStructs.FingerprintSimilarity(m1,m2, metric=DataStructs.DiceSimilarity)
            cos = DataStructs.FingerprintSimilarity(m1,m2, metric=DataStructs.CosineSimilarity)
            average = (tanimoto+dice+cos)/3
            f = open(output_csv_file_name, 'a+')
            f.write(k1 +',' + k2 +',' + str(tanimoto)+ ',' + str(dice) + ',' + str(cos)+ ',' + str(average)+ '\n')
            f.close()

