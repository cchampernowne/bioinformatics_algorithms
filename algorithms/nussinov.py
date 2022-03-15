import numpy as np
import sys
import pandas as pd


np.set_printoptions(threshold=np.inf)
struct_list = []


def main():
    seq = get_sequence(sys.argv[1])
    smatrix, pmatrix = create_matrices(seq)

    dot_bracket = []
    dot_bracket = ['-' for i in range(len(seq))]

    traceback(len(seq)-1,0,dot_bracket,smatrix,pmatrix,0,seq)

    global struct_list
    string_list = all_dot_brackets(struct_list)

    templ = []
    for ele in string_list:
        if ele[1] + round(struct_energy(ele[0]),1) < 0:
            templ.append([ele[0],ele[1] + round(struct_energy(ele[0]),1)])

    top_structs = sorted(templ, key=lambda ele: ele[1])[:10]

    same_scores = []
    for struct in top_structs:
        if struct[1] == top_structs[0][1]:
            same_scores.append(struct[0])

    # Q 5.a
    f = open("5.o1", "w+")
    f.write(seq + '\n' + str(top_structs[0][1]) + '\n' + top_structs[0][0])
    f.close()

    # Q 5.b
    f = open("5.o2", "w+")
    if len(same_scores) > 1:
        f.write('YES')
    else:
        f.write('NO')
    f.close()

     # Q 5.c
    to_write = str(len(same_scores))
    for struct in same_scores: to_write+= '\n' + struct
    f = open("5.o3", "w+")
    f.write(to_write)
    f.close()

    # Q 5.d
    to_write = seq
    for struct in top_structs: to_write+= '\n' + struct[0] + ' ' + str(struct[1])
    f = open("5.o4", "w+")
    f.write(to_write)
    f.close()


def struct_energy(dot_bracket):  # predict the loop penalty for a given structure
    score = 0
    for ele in dot_bracket:
        if ele == '(':
            break
        dot_bracket = dot_bracket[1:]

    for ele in dot_bracket[::-1]:
        if ele == ')':
            break
        dot_bracket = dot_bracket[:-1]
    
    score+=.1*dot_bracket.count('.')
    return score


def loop_count(dot_bracket):  # check if structure has viable counts of bases between loops
    loop_list = []
    if (dot_bracket.count('(') != dot_bracket.count(')')) or '(' not in dot_bracket or ')' not in dot_bracket:
        return False

    start = -1
    for i in range(len(dot_bracket)):
        if (dot_bracket[i] == '(' or dot_bracket[i] == ')') and start == -1:
            start = i
        elif (dot_bracket[i] == '(' or dot_bracket[i] == ')') and start != -1:
            if i - start < 4 and i - start > 0:
                loop_list.append(False)
            else:
                loop_list.append(True)
            start = i

    if False in loop_list:
        return False
    return True


def all_dot_brackets(struct_list):  # change to string format, remove duplicate dot brackets, remove inelligible dot brackets using loop_count()
    string_list = []
    for struct in struct_list:
        dot_bracket = ''
        for ele in struct[0]:
            dot_bracket+=ele
        string_list.append([dot_bracket,struct[1]])
    temp1_list = []
    temp2_list = []

    for ele in string_list:
        if ele[0] not in temp1_list:
            temp1_list.append(ele[0])
            temp2_list.append(ele)

    final_list = []
    for ele in temp2_list:
        if loop_count(ele[0]) == True:
            final_list.append(ele)
    
    return final_list


def get_sequence(file_name):  # retrieve sequences from file and return a list of all given sequences
    f = open(file_name)
    raw = f.read()
    if '>' in raw:
        raw_list = raw.split('\n')
        f.close()
        raw_list.pop(0)
        sequence = ''.join(raw_list)
        return sequence
    return raw


def min_adj_val(i,j,smatrix,seq):  # return min of left, bottom, and if i and j are a base pair, bottom-left cell, as well as whether the current cell represent a pair
    pair = False
    if (seq[i] == 'A' and seq[j] == 'U') or (seq[i] == 'U' and seq[j] == 'A') or (seq[i] == 'G' and seq[j] == 'U') or (seq[i] == 'U' and seq[j] == 'G'):
        pair_score = -2
        pair = True
    elif (seq[i] == 'G' and seq[j] == 'C') or (seq[i] == 'C' and seq[j] == 'G'):
        pair_score = -3
        pair = True
    else:
        pair_score = np.inf

    if pair_score == min(smatrix.loc[j+1,i], smatrix.loc[j,i-1], pair_score):
        min_score = smatrix.loc[j+1,i-1] + pair_score
    else:
        min_score = min(smatrix.loc[j+1,i], smatrix.loc[j,i-1])

    return min_score, pair


def traceback(i,j,dot_bracket,smatrix,pmatrix,score,seq):  # traceback all possible cells and append to list of total structures
    if '-' not in dot_bracket:
        global struct_list
        struct_list.append([dot_bracket,score])
        return
    
    if pmatrix.loc[j,i] == True:
        dot_bracket[min(i,j)] = '('
        dot_bracket[max(i,j)] = ')'  
        if (seq[i] == 'A' and seq[j] == 'U') or (seq[i] == 'U' and seq[j] == 'A') or (seq[i] == 'G' and seq[j] == 'U') or (seq[i] == 'U' and seq[j] == 'G'):
            traceback(i-1,j+1,dot_bracket,smatrix,pmatrix,score-2,seq)
        else:
            traceback(i-1,j+1,dot_bracket,smatrix,pmatrix,score-3,seq)

    i_dot = dot_bracket.copy()
    j_dot = dot_bracket.copy()
    i_dot[i] = '.'
    j_dot[j] = '.'
    traceback(i-1,j,i_dot,smatrix,pmatrix,score,seq)
    traceback(i,j+1,j_dot,smatrix,pmatrix,score,seq)


def create_matrices(seq):  # create matrix using nussinov algorithm
    smatrix = pd.DataFrame(np.full((len(seq),len(seq)), -np.inf), index=range(len(seq)), columns=range(len(seq)))
    pmatrix = pd.DataFrame(np.full((len(seq),len(seq)), -np.inf), index=range(len(seq)), columns=range(len(seq)))

    i = 0
    j = 0
    while j < (len(seq)):
        smatrix.loc[j,i] = 0
        j+=1
        if j < len(seq):
            smatrix.loc[j,i] = 0
        i+=1

    col_list = [*range(1,len(seq))]
    row_list = [*range(0,len(seq) - 1)]

    while row_list:    
        for i in range(len(col_list)):
            score,pair = min_adj_val(col_list[i],row_list[i],smatrix,seq)
            smatrix.loc[row_list[i],col_list[i]] = score
            pmatrix.loc[row_list[i],col_list[i]] = pair
        del row_list[-1]
        del col_list[0]

    return smatrix, pmatrix


if __name__ == "__main__":
    main()