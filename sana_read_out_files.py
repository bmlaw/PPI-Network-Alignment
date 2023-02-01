import os
from matplotlib import pyplot as plt
import sys


def read_out(filename):
    ofile = open(filename, 'r')

    ec = ""
    s3 = ""
    sequence = ""
    weight = ""

    for line in ofile:
        if line[0:3] == "ec:":
            ec = line[4:]
            
        if line[0:3] == "s3:":
            s3 = line[4:]

        if line[0:9] == "sequence:":
            sequence = line[10:]

        if line[0:16] == "weight sequence:":
            weight = line[17:]

    ofile.close()

    return [ec, s3, sequence, weight]


def read_directory(directory):

    ec_list = []
    s3_list = []
    sequence_list = []
    weight_list = []

    for filename in os.listdir(directory):
        out_list = read_out(directory + '/' + filename)
        ec_list.append(float(out_list[0].strip()))
        s3_list.append(float(out_list[1].strip()))
        sequence_list.append(float(out_list[2].strip()))
        weight_list.append(float(out_list[3].strip()))

    return [ec_list, s3_list, sequence_list, weight_list]

'''
def get_all_seq_scores(directory):

    # weight seq -> [unweighted_score, weighted_score]
    score_dict = dict()

    for filename in os.listdir(directory):
        unweighted = float(filename[1].split(':').strip()[1:])
        weighted = float(filename[2].split(':').strip()[1:])
        score_dict[float(filename.split('_')[3])] = [unweighted, weighted]

    return score_dict
'''

# read files names
def get_ec(filename):

    ec = filename.split('_')[2]


    return ec


'''
Below are for finding high seq matches!
'''


def read_align(filename):
    align = open(filename, 'r')

    # [left, right]
    match_list = []

    for line in align:
        by_tab = line.strip().split('\t')
        left = by_tab[0]
        right = by_tab[1]

        match_list.append([left, right])

    align.close()

    return match_list

'''
def read_blast_score(filename):
    ofile = open(filename, 'r')

    # left$right -> score
    match_dict = dict()

    # left -> highest score
    highest_dict = dict()

    for line in ofile:
        by_tab = line.strip().split('\t')

        if len(by_tab) == 12:
            match_dict[by_tab[0] + '$' + by_tab[1]] = float(by_tab[11])

            if by_tab[0] not in highest_dict or float(by_tab[11]) > highest_dict[by_tab[0]]:
                highest_dict[by_tab[0]] = float(by_tab[11])

            if by_tab[1] not in highest_dict or float(by_tab[11]) > highest_dict[by_tab[1]]:
                highest_dict[by_tab[1]] = float(by_tab[11])
                

    ofile.close()

    return [match_dict, highest_dict]
'''

def read_blast_score(filename):
    ofile = open(filename, 'r')

    # left$right -> score
    match_dict = dict()

    # left -> highest score
    highest_dict = dict()

    # {node1 -> [[node2, score1],[...]], node3 ->[...],...}
    rbh_dict = dict()

    for line in ofile:
        by_tab = line.strip().split('\t')


        if len(by_tab) == 12:
            match_dict[by_tab[0] + '$' + by_tab[1]] = float(by_tab[11])

            if by_tab[0] not in highest_dict or float(by_tab[11]) > highest_dict[by_tab[0]]:
                highest_dict[by_tab[0]] = float(by_tab[11])

            if by_tab[1] not in highest_dict or float(by_tab[11]) > highest_dict[by_tab[1]]:
                highest_dict[by_tab[1]] = float(by_tab[11])

            left = by_tab[0]
            right = by_tab[1]
            score = float(by_tab[11])

            if (left in rbh_dict) and (right not in rbh_dict):
                # deal with equal highest scores
                for sub in rbh_dict[left]:
                    if score > sub[1]:
                        # remove the low score keys
                        rbh_dict[left].remove(sub)
                        rbh_dict[left].append([right, score])
                        
                    elif score == sub[1]:
                        # coexist of highest score pairs
                        rbh_dict[left].append([right, score])


            elif (left not in rbh_dict) and (right not in rbh_dict):
                rbh_inner = []
                rbh_inner.append([right, score])
                rbh_dict[left] = rbh_inner

                rbh_inner = []
                rbh_inner.append([left, score])
                rbh_dict[right] = rbh_inner
                            

            elif (right in rbh_dict) and (left not in rbh_dict):
                # deal with equal highest scores
                for sub in rbh_dict[right]:
                    if score > sub[1]:
                        # remove the low score keys
                        rbh_dict[right].remove(sub)
                        rbh_dict[right].append([left, score])
                        
                    elif score == sub[1]:
                        # coexist of highest score pairs
                        rbh_dict[right].append([left, score])                
                
            

    rbh_list = []

    for node in rbh_dict:
        for sub in rbh_dict[node]:
            if sub[0] in rbh_dict:
                # [[score, node0, node1], [...], ...]
                rbh_list.append([sub[1], node, sub[0]])

                # remove pairs from both sides
                rbh_dict[node].remove(sub)
                for sub_list in rbh_dict[sub[0]]:
                    if sub_list[0] == node:
                        rbh_dict[sub[0]].remove(sub_list)

                '''
                if len(rbh_dict[node]) == 0:
                    rbh_dict.pop(node)

                if len(rbh_dict[sub[0]]) == 0:
                    rbh_dict.pop(sub[0])
                '''

    rbh_list.sort(reverse = True)
                         
                

    ofile.close()

    return [match_dict, highest_dict, rbh_list]


def find_high_score_match(align_list, blast_score_dict):

    seq_match_list = []

    for match in align_list:
        #if blast_score_dict.in(match[0] + '$' + match[1]):
        if match[0] + '$' + match[1] in blast_score_dict:
            seq_match_list.append([match[0], match[1], blast_score_dict[match[0] + '$' + match[1]]])

    return seq_match_list


# [unweighted_avg, weight_avg]
def calculate_seq_scores(align_list, score_dict, highest_dict):
    total = 0
    total_highest = 0
    unweighted_avg = 0
    
    for match in align_list:
        #print('TEST:', match[0] + '$' + match[1])
        if match[0] + '$' + match[1] in score_dict:
            print('TEST2')
            '''
            if match[0] in highest_dict:
                max_score = highest_dict[match[0]]
            else:
                max_score = 1
            '''
            max_score = highest_dict[match[0]]
            unweighted_avg += float(score_dict[match[0] + '$' + match[1]]) / max_score
            
            total += float(score_dict[match[0] + '$' + match[1]])
            
        if match[0] in highest_dict:
            total_highest += highest_dict[match[0]]

    weighted_avg = total / max(total_highest, 1)
    unweighted_avg = unweighted_avg / len(align_list)

    print(unweighted_avg, weighted_avg)

    return [unweighted_avg, weighted_avg]


# reads alignment directory
def get_correct_seq_scores(directory, score_dict, highest_dict):
    # seq_weight -> [unweighted, weighted]
    return_me = dict()

    for filename in os.listdir(directory):
        seq = 1 - float(get_ec(filename))

        align_list = read_align(directory + '/' + filename)

        seq_scores = calculate_seq_scores(align_list, score_dict, highest_dict)

        return_me[seq] = seq_scores

    return return_me


        
            
        
'''
'''

    
'''
files = sys.argv[1] # <-- directory contains sana out files
align_files = sys.argv[2] # <-- directory of align files
blast_file = sys.argv[3] # <-- blast out file

blast = read_blast_score(blast_file)
match_dict = blast[0]
highest_dict = blast[1]

print(match_dict)
print(highest_dict)

score_dict = get_correct_seq_scores(align_files, match_dict, highest_dict)

print(score_dict)

plot_list = read_directory(files)

ec_list = plot_list[0]
s3_list = plot_list[1]
sequence_list = plot_list[2]
weight_list = plot_list[3]

unweighted_seq_list = []
weighted_seq_list = []


for seq in score_dict:
    unweighted_seq_list.append(score_dict[seq][0])
    weighted_seq_list.append(score_dict[seq][1])

plt.plot(weight_list, ec_list, label = 'ec', linestyle='--', marker='o')
plt.plot(weight_list, s3_list, label = 's3', linestyle='--', marker='o')
plt.plot(weight_list, sequence_list, label = 'sequence similarity score', linestyle='--', marker='o')
plt.plot(weight_list, unweighted_seq_list, label = 'unweighted sequence similarity score', linestyle='--', marker='o')
plt.plot(weight_list, weighted_seq_list, label = 'weighted sequence similarity score', linestyle='--', marker='o')

plt.xlabel('sequence weight')
plt.ylabel('score')
plt.title('SANA Scores vs. Sequence Weight')

plt.legend()

plt.show()
'''

blast_file = "RNorvegicus_AThaliana_blast.out"
print(read_blast_score(blast_file)[2])

#print(plot_list)



