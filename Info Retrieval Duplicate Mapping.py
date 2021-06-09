'''
Input: {
    qid1:int,
    qid2:int,
    question1:str, # corresponding to qid1
    question2:str, # corresponding to qid2
    is_duplicate:bool
}

Output: [
    set(qids that are duplicates)
]

If is_duplicate=0 but qids are in the same set,
e.g. (1, 2) are duplicates, (2, 3) are duplicates, but (1, 3) are not,
(1, 2, 3) are kept in the same set, but (1, 3) is logged in contradictions.txt
'''

import pandas as pd
from tqdm import tqdm

data_source = 'train.csv'
output_file = 'duplicates.txt'
contradiction_file = 'contradictions.txt'

df = pd.read_csv(data_source)

set_num = 0
qid_dict = dict() # key: qid, value: set_num
dup_dict = dict() # key: set_num, value: [qids]
contradiction_count = 0 # cases when is_duplicate=0 but qids are in the same set

for row in tqdm(df.itertuples(index=False)):
    qid1, qid2, str1, str2 = row[1], row[2], row[3], row[4]
    is_duplicate = row[5] == 1
    if is_duplicate:
        if qid1 not in qid_dict and qid2 not in qid_dict:
            qid_dict[qid1] = set_num
            qid_dict[qid2] = set_num
            dup_dict[set_num] = [qid1, qid2]
            set_num += 1
        elif qid1 not in qid_dict: # implied: qid2 in qid_dict
            set_num_qid2 = qid_dict[qid2]
            qid_dict[qid1] = qid_dict[qid2]
            if set_num_qid2 in dup_dict.keys():
                dup_dict[set_num_qid2].append(qid1)
            else:
                dup_dict[set_num_qid2] = [qid1]
        elif qid2 not in qid_dict: # implied: qid1 in qid_dict
            set_num_qid1 = qid_dict[qid1]
            qid_dict[qid2] = qid_dict[qid1]
            if set_num_qid1 in dup_dict.keys():
                dup_dict[set_num_qid1].append(qid2)
            else:
                dup_dict[set_num_qid1] = [qid2]
        elif qid1 in qid_dict and qid2 in qid_dict:
            # contradiction: is_duplicate=1 but qids are in different sets
            # fix: merge the set qid2 is in into the set qid1 is in
            set_num_qid1 = qid_dict[qid1]
            set_num_qid2 = qid_dict[qid2]
            if set_num_qid1 != set_num_qid2:
                for _qid in dup_dict[set_num_qid2]:
                    dup_dict[set_num_qid1].append(_qid)
                    qid_dict[_qid] = set_num_qid1
                del dup_dict[set_num_qid2]
        else:
            print(f'Unhandled case: qid1={qid1}, qid2={qid2}')
    else:
        for qid in [qid1,qid2]:
            if qid not in qid_dict:
                qid_dict[qid] = set_num
                dup_dict[set_num] = [qid]
                set_num += 1
        if qid_dict[qid1] == qid_dict[qid2]:
            # contradiction: is_duplicate=0 but qids are in the same set
            # not fixable: likely oversights in dataset
            contradiction_count += 1
            f_contra = open(contradiction_file, 'a')
            f_contra.write(f'qid1={qid1}, qid2={qid2}, question1={str1}, question2={str2}\n')
            f_contra.close()

f_output = open(output_file, 'w+')
for dup_list in dup_dict.values():
    f_output.write(f'{set(sorted(dup_list))}\n')
f_output.close()
print(f'Output saved to {output_file}')

print(f'{contradiction_count} contradictions found where is_duplicate=0 but qids are in the same set.')
print(f'Contradictions saved to {contradiction_file}')
