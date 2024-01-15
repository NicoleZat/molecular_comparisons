########################################
# Author(s): Nicole Zatorski           #
# Date last updated: 1/15/24           #
# Specs: molecule target quality check #
########################################

from chembl_webresource_client.new_client import new_client
import pandas as pd
target = new_client.target

# whole_file = pd.read_csv('./TargetInfo/drug_targets.csv')
whole_file = pd.read_csv('drug_targets.csv')
targets_list = whole_file.columns.tolist()

kept_targ = []
rid_targ = []
for targ in targets_list[1:]:
    print('~~~~~~~~~~')
    print(targ)
    next_target = new_client.target.filter(target_chembl_id = targ)
    if next_target[0]['organism'] == 'Homo sapiens':
    # if next_target[0]['target_type'] != 'CELL-LINE' and next_target[0]['target_type'] != 'TISSUE' and next_target[0]['target_type'] != 'UNCHECKED' and next_target[0]['organism'] == 'Homo sapiens':
        kept_targ.append(targ)
    else:
        whole_file = whole_file.drop([targ], axis=1)
        rid_targ.append(targ)

whole_file.to_csv('targets_human.csv', index=False)
# whole_file.to_csv('targets_human_non_cell_tissue_unchecked.csv', index=False)
print(rid_targ)
print('*********')
print(kept_targ)
print('@@@@@@@') 
print(len(targets_list))
print(len(kept_targ))

