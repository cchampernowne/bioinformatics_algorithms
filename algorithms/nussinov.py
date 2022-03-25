import numpy as np
import pandas as pd
import re


np.set_printoptions(threshold=np.inf)
ERROR_MSG = '\nERROR!\nInput not understood.\n - Sequence lengths must be between 5 - 200\n - Sequences must not have spaces, or any other characters between nucleotides\n - Nucleotides must be represented using only the characters A, T, C, G, U, a, t, c, g, and u.\n\nPlease try again'


def main():
    seq = interperate_input(input('Enter sequence,\nor enter Q to quit: '))
    if not seq:
        print('GOODBYE :)')
        return
    smatrix, pmatrix = create_matrices(seq)
    dot_bracket = []
    dot_bracket = ['-' for i in range(len(seq))]
    struct_list = []
    traceback(len(seq)-1, 0, dot_bracket, smatrix, pmatrix, 0, seq, struct_list)
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
    print('\nRESULTS:\n\nSTRUCTURES:')
    for struct in top_structs:
        print(struct[0] + ' ' + str(struct[1]))


def server_results(seq):
    valid_sequence = re.compile(r'^ *(A|T|C|G|U|a|t|c|g|u)* *$')
    if re.fullmatch(valid_sequence, seq):
        smatrix, pmatrix = create_matrices(seq)
        dot_bracket = []
        dot_bracket = ['-' for i in range(len(seq))]
        struct_list = []
        traceback(len(seq)-1, 0, dot_bracket, smatrix, pmatrix, 0, seq, struct_list)
        string_list = all_dot_brackets(struct_list)
        templ = []
        for ele in string_list:
            if ele[1] + round(struct_energy(ele[0]),1) < 0:
                templ.append([ele[0],ele[1] + round(struct_energy(ele[0]),1)])
        top_structs = sorted(templ, key=lambda ele: ele[1])[:10]
        to_write = ''
        for struct in top_structs:
            to_write+= '\n' + struct[0] + ' ' + str(struct[1])
        return to_write
    return False


def interperate_input(usr_input):
    re_quit = re.compile(r'^ *(Q|q) *$')
    re_valid_seq = re.compile(r'^ *(A|T|C|G|U|a|t|c|g|u)* *$')
    if re.fullmatch(re_quit, usr_input):
        return False
    if re.fullmatch(re_valid_seq, usr_input) and 4 < len(usr_input.strip()) < 201:
        return usr_input.strip()
    return interperate_input(input(ERROR_MSG + ',\nor enter Q to quit:'))


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
    if (seq[i] in ('A','a') and seq[j] in ('U','u','T','t')) or (seq[i] in ('U','u','T','t') and seq[j] in ('A','a')) or (seq[i] in ('G','g') and seq[j] in ('U','u','T','t')) or (seq[i] in ('U','u','T','t') and seq[j] in ('G','g')):
        pair_score = -2
        pair = True
    elif (seq[i] in ('G','g') and seq[j] in ('C','c')) or (seq[i] in ('C','c') and seq[j] in ('G','g')):
        pair_score = -3
        pair = True
    else:
        pair_score = np.inf
    if pair_score == min(smatrix.loc[j+1,i], smatrix.loc[j,i-1], pair_score):
        min_score = smatrix.loc[j+1,i-1] + pair_score
    else:
        min_score = min(smatrix.loc[j+1,i], smatrix.loc[j,i-1])
    return min_score, pair


def traceback(i, j, dot_bracket, smatrix, pmatrix, score, seq, struct_list):  # traceback all possible cells and append to list of total structures
    if '-' not in dot_bracket:
        struct_list.append([dot_bracket,score])
        return
    if pmatrix.loc[j,i] == True:
        dot_bracket[min(i,j)] = '('
        dot_bracket[max(i,j)] = ')'  
        if (seq[i] == 'A' and seq[j] == 'U') or (seq[i] == 'U' and seq[j] == 'A') or (seq[i] == 'G' and seq[j] == 'U') or (seq[i] == 'U' and seq[j] == 'G'):
            traceback(i-1, j+1, dot_bracket, smatrix, pmatrix, score-2, seq, struct_list)
        else:
            traceback(i-1, j+1, dot_bracket, smatrix, pmatrix, score-3, seq, struct_list)
    i_dot = dot_bracket.copy()
    j_dot = dot_bracket.copy()
    i_dot[i] = '.'
    j_dot[j] = '.'
    traceback(i-1, j, i_dot, smatrix, pmatrix, score, seq, struct_list)
    traceback(i, j+1, j_dot, smatrix, pmatrix, score, seq, struct_list)


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


if __name__ == '__main__':
    main()