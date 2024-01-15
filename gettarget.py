######################################
# Author(s): Nicole Zatorski         #
# Date last updated: 1/15/24         #
# Specs: molecular target gen        #
######################################

from chembl_webresource_client.new_client import new_client
import pandas as pd

activities = new_client.activity
drugdf = pd.read_csv('atc_chemID_map.csv')

target_id_name_df = pd.DataFrame(columns=['ChEBML ID', 'Name'])
target_id_name_dict = {}
# target_id_set = set([])
drug_dict = {}

print('searching targets')
for index, row in drugdf.iterrows():
    chem_id = row['ChEMBL ID']
    print('Index: ' + str(index)+ ' Compound: ' + chem_id)
    res = activities.filter(molecule_chembl_id=chem_id)
    target_list = []
    for element in res:
        if element['type'] == 'IC50':
            if element['standard_value']:
                if float(element['standard_value']) <= 1000:
                    target_id = element['target_chembl_id']
                    target_name = element['target_pref_name']
                    target_list.append(target_id)
                    # target_id_set.add(target_id)
                    if target_id not in target_id_name_dict.keys():
                        target_id_name_dict[target_id] = target_name
                        temp = pd.DataFrame([[target_id, target_name]],columns=['ChEBML ID', 'Name'])
                        target_id_name_df = target_id_name_df.append(temp)
    drug_dict[chem_id] = target_list

target_id_name_df.to_csv('target_id_name_map.csv', index=False)
target_ids_list = list(target_id_name_dict.keys())
# target_ids_list = list(target_id_set)
out_cols = ['ChEMBL ID'] + target_ids_list
out = pd.DataFrame(columns=out_cols)

print('making output file')
for chem_id in drug_dict:
    print(chem_id)
    target_list = drug_dict[chem_id]
    target_for_drug = []
    for target in target_ids_list:
        if target in target_list:
            target_for_drug.append(1)
        else: 
            target_for_drug.append(0)
    tempdf = pd.DataFrame([[chem_id]+target_for_drug], columns=out_cols)
    out = out.append(tempdf)

out.to_csv('drug_targets.csv', index=False)