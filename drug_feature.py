######################################
# Author(s): Nicole Zatorski         #
# Date last updated: 1/15/24         #
# Specs: molecule features code      #
######################################

from rdkit import Chem
from rdkit.Chem import AllChem
from  rdkit.DataStructs import FingerprintSimilarity
from rdkit.Chem import FragmentCatalog
from rdkit import RDConfig
import os
import numpy as np

def make_m_dict():
    c2f_dict = {}
    c3f_dict = {}
    target_dict = {}
    atc_dict = {}
    sages_dict = {}

    with open('files_in_use/chem_to_fda.csv') as fo:
        for line in fo:
            split_line = line[:-1].split(',')
            c2f_dict[split_line[1]] = split_line[0]
            c3f_dict[split_line[0]] = split_line[1]

            drug_targets = {}
            try:
                with open('all_drugs_weighted/'+ split_line[0] + '.csv') as fo2:
                    for line in fo2:
                        split_line2 = line[:-1].split(',')
                        if split_line2[0] not in drug_targets:
                            drug_targets[split_line2[0]] = [float(split_line2[1])]
                        else:
                            drug_targets[split_line2[0]] = drug_targets[split_line2[0]] + [float(split_line2[1])]
                out_targets = {}
                for t in drug_targets:
                    if len(drug_targets[t]) > 1:
                        valu = sum(drug_targets[t])/len(drug_targets[t])
                    else:
                        valu = drug_targets[t][0]
                    # if valu <= 1000:
                    out_targets[t] = valu
                target_dict[split_line[1]] = out_targets

                drug_sages = {}
                with open('sages_all_drugs_out/frequency_'+ split_line[0] +'.csv') as fo3:
                    lineindex = 0
                    for line in fo3:
                        if lineindex >0:
                            split_line3 = line[:-1].split(',')
                            if float(split_line3[3]) != 0:
                                drug_sages[split_line3[1]] = float(split_line3[3])
                        lineindex += 1
                with open('sages_all_drugs_out/average_'+ split_line[0] +'.csv') as fo3:
                    lineindex = 0
                    for line in fo3:
                        if lineindex > 0:
                            split_line3 = line[:-1].split(',')
                            if float(split_line3[1]) != 0:
                                drug_sages[split_line3[0]] = float(split_line3[1])
                        lineindex += 1
                sages_dict[split_line[1]] = drug_sages
            except:
                pass

    smiles_dict = {}
    with open('files_in_use/drug_info.csv') as fo:
        for line in fo:
            split_line = line[:-1].split(',')
            fpout = []
            smiles_dict[split_line[0]] = split_line[30]
            atcs = split_line[1].split('_')
            for elt in atcs:
                if len(elt) > 0:
                    if split_line[0] in c3f_dict:
                        if c3f_dict[split_line[0]] not in atc_dict:
                            atc_dict[c3f_dict[split_line[0]]] = [elt[:3]]
                        else:
                            atc_dict[c3f_dict[split_line[0]]] = atc_dict[c3f_dict[split_line[0]]] + [elt[:3]]
    fp_dict = {}
    with open('files_in_use/clean_drug_fps.csv') as fo:
        for line in fo:
            split_line = line[:-1].split(',')
            fpout = []
            for elt in split_line[1:]:
                fpout.append(int(elt))
            fp_dict[split_line[0]] = smiles_dict[c2f_dict[split_line[0]]]

    return fp_dict, target_dict, atc_dict, sages_dict

def add_withdrawn():
    fp_dict, target_dict, atc_dict, sages_dict = make_m_dict()

    with open('files_in_use/withdrawn.csv') as fo:
        for line in fo:
            split_line = line[:-1].split(',')
            fp_dict[split_line[0]] = split_line[31]
            
            atc_out = []
            atcs = split_line[2].split('_')
            for atc in atcs:
                if len(atc) > 1:
                    atc_out.append(atc[:3])
            atc_dict[split_line[0]] = atc_out
            
            drug_targets = {}
            try:
                with open('all_drugs_weighted/'+ split_line[0] + '.csv') as fo2:
                    for line in fo2:
                        split_line2 = line[:-1].split(',')
                        if split_line2[0] not in drug_targets:
                            drug_targets[split_line2[0]] = [float(split_line2[1])]
                        else:
                            drug_targets[split_line2[0]] = drug_targets[split_line2[0]] + [float(split_line2[1])]
                out_targets = {}
                for t in drug_targets:
                    if len(drug_targets[t]) > 1:
                        valu = sum(drug_targets[t])/len(drug_targets[t])
                    else:
                        valu = drug_targets[t][0]
                    out_targets[t] = valu

                target_dict[split_line[0]] = out_targets

                drug_sages = {}
                with open('sages_all_drugs_out/frequency_'+ split_line[0] +'.csv') as fo3:
                    lineindex = 0
                    for line in fo3:
                        if lineindex >0:
                            split_line3 = line[:-1].split(',')
                            if float(split_line3[3]) != 0:
                                drug_sages[split_line3[1]] = float(split_line3[3])
                        lineindex += 1
                with open('sages_all_drugs_out/average_'+ split_line[0] +'.csv') as fo3:
                    lineindex = 0
                    for line in fo3:
                        if lineindex > 0:
                            split_line3 = line[:-1].split(',')
                            if float(split_line3[1]) != 0:
                                drug_sages[split_line3[0]] = float(split_line3[1])
                        lineindex += 1
                sages_dict[split_line[1]] = drug_sages

            except:
                pass
    return fp_dict, target_dict, atc_dict, sages_dict

def get_header(target_dict):
    all_targets = []
    for key in target_dict:
        all_targets = all_targets + list(target_dict[key].keys())
    all_targets = list(set(all_targets))
    all_targets.sort()
    return all_targets

def write_file(filename, str_to_write, write=True):
    if write:
        smiles_file = open(filename, 'a+')
        smiles_file.write(str_to_write)
        smiles_file.close()

def make_data_str_from_dict(drug, header_list, drug_dict):
    temp_targets = drug
    for header in header_list:
        if header in drug_dict:
            temp_targets = temp_targets + ',' + str(drug_dict[header])
        else:
            temp_targets = temp_targets + ',0'
    return temp_targets + '\n'

def dict_check(key_name, dict_name):
    if key_name not in dict_name:
        return {}
    else:
        return dict_name[key_name]


def get_feature_files(write=True):
    fp_dict, target_dict, atc_dict, sages_dict = add_withdrawn()
    indir = 'files_in_use/'
    
    all_targets = get_header(target_dict)
    all_sages = get_header(sages_dict)
    all_atc = ['A01', 'A02', 'A03', 'A04', 'A05', 'A06', 'A07', 'A08', 'A10', 'A11', 'A12', 'A14', 'A16', 'B01', 'B02', 'B05', 'C01', 'C02', 'C03', 'C04', 'C05', 'C07', 'C08', 'C09', 'C10', 'D01', 'D04', 'D05', 'D06', 'D07', 'D08', 'D09', 'D10', 'D11', 'G01', 'G02', 'G03', 'G04', 'H01', 'H02', 'H03', 'H04', 'H05', 'J01', 'J02', 'J04', 'J05', 'L01', 'L02', 'L03', 'L04', 'M01', 'M02', 'M03', 'M04', 'M05', 'M09', 'N01', 'N02', 'N03', 'N04', 'N05', 'N06', 'N07', 'P01', 'P02', 'P03', 'R01', 'R02', 'R03', 'R05', 'R06', 'R07', 'S01', 'S02', 'S03', 'V03', 'V04']

    write_file(indir + 'targets.csv','drug,' + ','.join(all_targets) + '\n', write=write)
    write_file(indir + 'sages.csv','drug,' + ','.join(all_sages)+ '\n', write=write)
    write_file(indir + 'atcs.csv','drug,' + ','.join(all_atc)+ '\n', write=write)

    for fp in fp_dict:
        write_file(indir + 'smiles.csv',fp + ',' + fp_dict[fp] + '\n', write=write)
        m1 = Chem.MolFromSmiles(fp_dict[fp])
        fpstring = str(np.array(AllChem.GetMorganFingerprintAsBitVect(m1, radius=2, nBits = 512))).replace(' ',',')[1:-1]
        write_file(indir + 'fp.csv', fp + ',' + fpstring + '\n', write=write)

        temp_targets = make_data_str_from_dict(fp, all_targets, dict_check(fp, target_dict))
        write_file(indir + 'targets.csv', temp_targets, write=write)

        temp_sages = make_data_str_from_dict(fp, all_sages, dict_check(fp, sages_dict))
        write_file(indir + 'sages.csv', temp_sages, write=write)

        write_file(indir + 'atcs.csv', fp, write=write)
        temp_atc_dict = dict_check(fp, atc_dict)
        for atc in all_atc:
            if atc in temp_atc_dict:
                write_file(indir + 'atcs.csv', ',1', write=write)
            else:
                write_file(indir + 'atcs.csv', ',0', write=write)
        write_file(indir + 'atcs.csv', '\n', write=write)
        

# get_feature_files(write=True)

with open('files_in_use/smiles.csv') as fo:
    for line in fo:
        split_line = line[:-1].split(',')
        m1 = Chem.MolFromSmiles(split_line[1])
        fp = str(list(AllChem.GetMorganFingerprintAsBitVect(m1, radius=2, nBits = 512))).replace(' ','')[1:-1]
        out = split_line[0] + ',' + fp + '\n'
        f=open('files_in_use/fp.csv','a+')
        f.write(out)
        f.close()








def make_cluster_dict(label_file):
    label_dict = {}
    out_dict = {}
    with open(label_file) as fo:
        for line in fo:
            split_line = line[:-1].split(',')
            if split_line[0] not in label_dict:
                label_dict[split_line[0]] = [int(split_line[1])]
            else:
                if label_dict[split_line[0]] != split_line[1]:
                    label_dict[split_line[0]] = label_dict[split_line[0]] + [int(split_line[1])]
    for l in label_dict:
        if len(label_dict[l]) > 1:
            valu = max(set(label_dict[l]), key=label_dict[l].count)
        else:
            valu = label_dict[l][0]
        if valu not in out_dict:
            out_dict[valu] = [l]
        else:
            out_dict[valu] = out_dict[valu] + [l]
    return(out_dict)



def cluster_drug_stats(cluster_file):
    lens = []
    all_drug_fps = []
    background_structs = {}
    by_cluster_structs = {}
    background_targets = {}
    outdir = 'cluster_frequencies_agglo/'

    fName=os.path.join(RDConfig.RDDataDir,'FunctionalGroups.txt')
    drug_cluster_dict = make_cluster_dict(cluster_file)
    m_dict,target_dict, atc_dict = make_m_dict()

    for k in drug_cluster_dict:
        cluster_fps = []
        cluster_structs = {}
        cluster_targets = {}
        cluster_atcs = {}

        for drug in drug_cluster_dict[k]:
            m1 = Chem.MolFromSmiles(m_dict[drug])
            fp = Chem.RDKFingerprint(m1)
            cluster_fps.append(fp)
            all_drug_fps.append(fp)

            if drug in target_dict:
                if len(target_dict) > 0:
                    for t in target_dict[drug]:
                        if t not in background_targets:
                            background_targets[t] = 1
                        else:
                            background_targets[t] = background_targets[t] + 1
                        if t not in cluster_targets:
                            cluster_targets[t] = 1
                        else:
                            cluster_targets[t] = cluster_targets[t] + 1
            
            if drug in atc_dict:
                atcs = atc_dict[drug]
                for a in atcs:
                    if a not in cluster_atcs:
                        cluster_atcs[a] = 1
                    else:
                        cluster_atcs[a] = cluster_atcs[a] + 1


            fcgen=FragmentCatalog.FragCatGenerator()
            fparams = FragmentCatalog.FragCatParams(1,6,fName)
            fcat=FragmentCatalog.FragCatalog(fparams)
            fcgen=FragmentCatalog.FragCatGenerator()
            num_fragments = fcgen.AddFragsFromMol(m1,fcat)
            for n in range(num_fragments):
                fragment = fcat.GetEntryDescription(n)
                if fragment not in cluster_structs:
                    cluster_structs[fragment] = 1
                else:
                    cluster_structs[fragment] = cluster_structs[fragment] + 1
                if fragment not in background_structs:
                    background_structs[fragment] = 1
                else:
                    background_structs[fragment] = background_structs[fragment] + 1 


        tanimoto = []
        for i in range(len(cluster_fps)):
            fp1 = cluster_fps[i]
            for index in range(i+1, len(cluster_fps)):
                fp2 = cluster_fps[index]
                tanimoto.append(FingerprintSimilarity(fp1, fp2))
        tanimoto = sum(tanimoto)/len(tanimoto)

        lens.append(len(drug_cluster_dict[k]))
        by_cluster_structs[k] = cluster_structs

        # f = open(outdir + str(k) + '_num_drugs_tanimoto.csv', 'a+')
        # f.write(str(len(drug_cluster_dict[k])) + ',' + str(tanimoto))
        # f.close()

        f1 = open(outdir + str(k) + '_atc.csv', 'a+')
        for a in cluster_atcs:
            f1.write(a + ',' + str(cluster_atcs[a])+ '\n')
        f1.close()

        # f2 = open(outdir + str(k) + '_cluster_targets.csv', 'a+')
        # for target in cluster_targets:
        #     f2.write(target+ ',' +str(cluster_targets[target]) + '\n')
        # f2.close()


    all_tan = []
    for i in range(len(all_drug_fps)):
        fp1 = all_drug_fps[i]
        for index in range(i+1, len(all_drug_fps)):
            fp2 = all_drug_fps[index]
            all_tan.append(FingerprintSimilarity(fp1, fp2))

    # f3 = open(outdir + 'background_num_drugs_tanimoto.csv', 'a+')
    # f3.write(str(sum(all_tan)/len(all_tan)) + ',' + str(min(lens))+ ',' + str(max(lens))+ ',' +str(sum(lens)/len(lens)))
    # f3.close()

    # for k in by_cluster_structs:
    #     f4 = open(outdir+str(k)+'_structs_background_removed.csv', 'a+')
    #     for struct in by_cluster_structs[k]:
    #         if struct in background_structs:
    #             f4.write(struct + ',' + str(by_cluster_structs[k][struct]/len(drug_cluster_dict[k])- background_structs[struct]/sum(lens)) + '\n')
    #         else:
    #             f4.write(struct + ',' + str(by_cluster_structs[k][struct]) + '\n')
    #     f4.close()
    
    # f5 = open(outdir + '/background_targets.csv', 'a+')
    # for target in background_targets:
    #     f5.write(target+ ',' + str(background_targets[target]) + '\n')
    # f5.close()

# cluster_drug_stats('files_in_use/nodupes_cluster20Agglo.csv')





















def get_tox_labels(map_file, input_file):
    map_list = []
    with open(map_file) as fo:
        for line in fo:
            split_line = line[:-1].split(',')
            map_list.append(split_line[1])

    all_words = set([])
    with open(input_file) as fo:
        for line in fo:
            split_line = line[:-1].split(',')
            if split_line[0].lower() in map_list:
                all_words = all_words.union(set(split_line))
    f = open('all_words.csv', 'a+')
    for word in all_words:
        try:
            float(word)
        except:
            if len(word) > 2:
                f.write(word + '\n')
    f.close()
# get_tox_labels('chem_to_fda.csv', 'squeaky_clean_labels.csv')



def map_id_to_fda():
    map_dict = {}
    with open('drug_info.csv') as fo:
        for line in fo:
            split_line = line[:-1].split(',')
            drug_id = split_line[0]
            names = split_line[-1].split('_')
            for name in names:
                if name != '':
                    if name.lower() not in map_dict:
                        map_dict[name.lower()] = [drug_id]
                    else:
                        map_dict[name.lower()] = map_dict[name.lower()] + [drug_id]
    with open('8505_drug_label_words.csv') as fo:
        for line in fo:
            drug_name = line[:-1].split(',')[0].lower()
            if drug_name in map_dict:
                all_names = map_dict[drug_name]
                all_names.sort()
                f = open('chem_to_fda.csv', 'a+')
                f.write(all_names[-1] + ',' + drug_name + '\n')
                f.close()
            else:
                if len(drug_name.split(' ')) > 1 and 'and' not in drug_name:
                    all_names = []
                    for elt in drug_name.split(' '):
                        if elt in map_dict:
                            all_names = all_names + map_dict[elt]
                    if len(all_names) >= 1:
                        all_names.sort()
                        f = open('chem_to_fda.csv', 'a+')
                        f.write(all_names[-1] + ',' + drug_name + '\n')
                        f.close()


def get_protein_features_and_labels(label_file,feature_file, map_file):
    map_dict = {}
    with open(map_file) as fo:
        for line in fo:
            split_line = line[:-1].split(',')
            map_dict[split_line[0]] = split_line[1]

    feature_dict = {}
    with open(feature_file) as fo:
        for line in fo:
            split_line = line[:-1].split(',')
            chemid = split_line[0]
            if chemid in map_dict:
                feats = ','.join(split_line[2:30]).lower().replace('none','0').replace('true','1').replace('false','0').replace('neutral', '0,0,1,0').replace('acid', '1,0,0,0').replace('base', '0,1,0,0').replace('zwitterion', '0,0,0,1').replace('na', '0,0,0,0').replace('n','0').replace('y','1')
                smiles = split_line[30]
                m1 = Chem.MolFromSmiles(smiles)
                fp = str(list(AllChem.GetMorganFingerprintAsBitVect(m1, radius=2, nBits = 512))).replace(' ','')[1:-1]
                out = feats + ',' + fp
                feature_dict[map_dict[chemid]] = map_dict[chemid] + ',' + feats + ',' + fp
    
    label_dict = {}
    with open(label_file) as fo:
        for line in fo:
            split_line = line[:-1].split(',')
            if split_line[0] in feature_dict:
                if split_line[0] not in label_dict:
                    label_dict[split_line[0]] = [int(split_line[1])]
                else:
                    if label_dict[split_line[0]] != split_line[1]:
                        label_dict[split_line[0]] = label_dict[split_line[0]] + [int(split_line[1])]
    for l in label_dict:
        if len(label_dict[l]) > 1:
            valu = max(set(label_dict[l]), key=label_dict[l].count)
        else:
            valu = label_dict[l][0]
        f = open('drugfeatures_and_clusterlabels.csv', '+a')
        f.write(str(valu) + ',' + feature_dict[l] + '\n')
        f.close()

# get_protein_features_and_labels('cluster20warn2clusterAgglo.csv','drug_info.csv', 'chem_to_fda.csv')

def clean_drug_features(feature_file, map_file):
    map_dict = {}
    with open(map_file) as fo:
        for line in fo:
            split_line = line[:-1].split(',')
            if split_line[0] not in map_dict:
                map_dict[split_line[0]] = [split_line[1]]
            elif split_line[1] not in map_dict[split_line[0]]:
                map_dict[split_line[0]] = map_dict[split_line[0]] + [split_line[1]]

    feature_dict = {}
    with open(feature_file) as fo:
        for line in fo:
            split_line = line[:-1].split(',')
            chemid = split_line[0]
            if chemid in map_dict:
                feats = ','.join(split_line[2:30]).lower().replace('none','0').replace('true','1').replace('false','0').replace('neutral', '0,0,1,0').replace('acid', '1,0,0,0').replace('base', '0,1,0,0').replace('zwitterion', '0,0,0,1').replace('na', '0,0,0,0').replace('n','0').replace('y','1')
                smiles = split_line[30]
                m1 = Chem.MolFromSmiles(smiles)
                fp = str(list(AllChem.GetMorganFingerprintAsBitVect(m1, radius=2, nBits = 512))).replace(' ','')[1:-1]
                for elt in set(map_dict[chemid]):
                    f = open('clean_drug_info.csv', '+a')
                    f.write(elt + ',' + feats + ',' + fp + '\n')
                    f.close()
                    f = open('clean_drug_fps.csv', '+a')
                    f.write(elt + ',' + fp + '\n')
                    f.close()

# clean_drug_features('drug_info.csv', 'chem_to_fda.csv')

# def clean_drug_features(feature_file, map_file):
#     map_chemid = []
#     map_names = []
#     with open(map_file) as fo:
#         for line in fo:
#             split_line = line[:-1].split(',')
#             map_chemid.append(split_line[0])
#             map_names.append(split_line[1])
#     print(len(map_chemid))
#     print(len(set(map_chemid)))
#     print(len(map_names))
#     print(len(set(map_names)))

#     label_names = []
#     with open(feature_file) as fo:
#         for line in fo:
#             split_line = line[:-1].split(',')
#             chemid = split_line[0]
#             label_names.append(chemid)
#     print(len(label_names))
#     print(len(set(label_names)))

# clean_drug_features('8505_drug_label_vectors.csv', 'chem_to_fda.csv')