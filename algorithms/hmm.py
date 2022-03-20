import numpy as np
import sys
import pandas as pd
import re


np.set_printoptions(threshold=np.inf)


def main():
    seq = get_sequence(sys.argv[1])
    fprob, hpost, lpost = forward_algorithm(seq)  # 'GGCA' for provided output files
    lvmatrix = log_viterbi_algorithm_matrix(seq)
    paths, mul = most_prob_path(lvmatrix,seq)
    vprob = viterbi_algorithm(seq)

    # Q 1)a
    f = open("4.o1", "w+")
    f.write("{:.2e}".format(fprob).replace("e-0", "e-"))
    f.close()

    # Q 1)b
    f = open("4.o2", "w+")
    f.write(print_matrix(lvmatrix,seq))
    f.close()

    # Q 1)c
    f = open("4.o3", "w+")
    f.write(paths[0])
    f.close()

    # Q 1)d
    f = open("4.o4", "w+")
    f.write("{:.2e}".format(vprob).replace("e-0", "e-"))
    f.close()

    # Q 1)e
    f = open("4.o5", "w+")
    f.write(mul)
    f.close()

    # Q 1)f
    f = open("4.o6", "w+")
    f.write("{:.2e}".format(hpost).replace("e-0", "e-") + '\n' + "{:.2e}".format(lpost).replace("e-0", "e-"))
    f.close()

    # Q 1)g
    f = open("4.o7", "w+")
    f.write(str(len(paths)) + '\n')
    for path in paths:
        f.write(path + '\n')
    f.close()


def server_results(seq):
    valid_sequence = re.compile(r"^(A|T|C|G|U)*$")
    if re.fullmatch(valid_sequence, seq):
        fprob, hpost, lpost = forward_algorithm(seq)  # 'GGCA' for provided output files
        lvmatrix = log_viterbi_algorithm_matrix(seq)
        paths, mul = most_prob_path(lvmatrix,seq)
        vprob = viterbi_algorithm(seq)
        x_probability = "{:.2e}".format(fprob).replace("e-0", "e-")
        viterbi = print_matrix(lvmatrix,seq)
        most_prob_paths = ""
        for path in paths:
            most_prob_paths = most_prob_paths + path + '\n\n'
        most_prob_probability = "{:.2e}".format(vprob).replace("e-0", "e-")
        return [x_probability, viterbi, most_prob_paths, most_prob_probability]
    return False


def print_matrix(matrix, seq):  # return the smatrix as a string for q 1)b
    total_matrix = '- 0 '
    for char in seq: total_matrix+= char + ' '
    total_matrix = total_matrix[:-1]
    total_matrix+= '\n' + '0 0.0 '
    for i in range(len(seq)): total_matrix+= '-INF '
    total_matrix = total_matrix[:-1]
    total_matrix+='\nH -INF '
    for column in range(len(seq)): total_matrix+= str(round(matrix.loc['H',column],1)) + ' '
    total_matrix = total_matrix[:-1]
    total_matrix+='\nL -INF '
    for column in range(len(seq)): total_matrix+= str(round(matrix.loc['L',column],1)) + ' '
    total_matrix = total_matrix[:-1]
    return total_matrix


def get_sequence(file_name):  # retrieve sequence from file and return
    f = open(file_name)
    raw_list = f.read().split('\n')
    f.close()
    raw_list.pop(0)
    sequence = ''.join(raw_list)
    return sequence


def most_prob_path(vmatrix,seq):  # calculate the most probable path as a result of viterbi algorithm. if final column in matrix has same vals, multiple paths can occur
    mpp = ''
    ms = 'NO'
    if vmatrix.loc['L',len(seq) - 1] == vmatrix.loc['H',len(seq) - 1]:
        ms = 'YES'
        for pos in vmatrix.columns:
            if vmatrix.loc['H',pos] > vmatrix.loc['L',pos]:
                mpp+='H'
            else:
                mpp+='L'
        mpp1 = mpp[:-1]
        mpp2 = mpp1 + 'H'
        mpp1+= 'L'

        return [mpp1,mpp2],ms

    for pos in vmatrix.columns:
        if vmatrix.loc['H',pos] > vmatrix.loc['L',pos]:
            mpp+='H'
        else:
            mpp+='L'
    return [mpp],ms


def hmm_probability(init_state, end_state, residue):  # hardcoded probabilities as given in assignment
    prob_score = 0
    if init_state == 'init':

        if end_state == 'H':
            if residue == 'A':
                prob_score = 0.5*0.2
            if residue == 'C':
                prob_score = 0.5*0.3
            if residue == 'G':
                prob_score = 0.5*0.3
            if residue == 'T':
                prob_score = 0.5*0.2

        if end_state == 'L':
            if residue == 'A':
                prob_score = 0.5*0.3
            if residue == 'C':
                prob_score = 0.5*0.2
            if residue == 'G':
                prob_score = 0.5*0.2
            if residue == 'T':
                prob_score = 0.5*0.3

    if init_state == 'H':

        if end_state == 'H':
            if residue == 'A':
                prob_score = 0.5*0.2
            if residue == 'C':
                prob_score = 0.5*0.3
            if residue == 'G':
                prob_score = 0.5*0.3
            if residue == 'T':
                prob_score = 0.5*0.2

        if end_state == 'L':
            if residue == 'A':
                prob_score = 0.5*0.3
            if residue == 'C':
                prob_score = 0.5*0.2
            if residue == 'G':
                prob_score = 0.5*0.2
            if residue == 'T':
                prob_score = 0.5*0.3

    if init_state == 'L':

        if end_state == 'H':
            if residue == 'A':
                prob_score = 0.4*0.2
            if residue == 'C':
                prob_score = 0.4*0.3
            if residue == 'G':
                prob_score = 0.4*0.3
            if residue == 'T':
                prob_score = 0.4*0.2

        if end_state == 'L':
            if residue == 'A':
                prob_score = 0.6*0.3
            if residue == 'C':
                prob_score = 0.6*0.2
            if residue == 'G':
                prob_score = 0.6*0.2
            if residue == 'T':
                prob_score = 0.6*0.3

    return prob_score


def forward_algorithm(seq):  # performs forward algorithm
    states = ['H','L']
    fmatrix = pd.DataFrame(np.full((2,len(seq)), -np.inf), index=states, columns=range(len(seq)))
    for i in range(len(seq)):
        for state in states:
            if i == 0:
                fmatrix.loc[state,i] = hmm_probability('init', state, seq[i])
            else:
                fmatrix.loc[state,i] = hmm_probability('H', state, seq[i])*fmatrix.loc['H',i-1] + hmm_probability('L', state, seq[i])*fmatrix.loc['L',i-1]
    return fmatrix.loc['L',len(seq)-1] + fmatrix.loc['H',len(seq)-1], fmatrix.loc['H',3] / (fmatrix.loc['L',3] + fmatrix.loc['H',3]), fmatrix.loc['L',3] / (fmatrix.loc['L',3] + fmatrix.loc['H',3])


def viterbi_algorithm(seq):  # perform viterbi alg with normal outputs
    states = ['H','L']
    vmatrix = pd.DataFrame(np.full((2,len(seq)), -np.inf), index=states, columns=range(len(seq)))
    for i in range(len(seq)):
        for state in states:
            if i == 0:
                vmatrix.loc[state,i] = hmm_probability('init', state, seq[i])
                continue
            if vmatrix.loc['L',i-1] < vmatrix.loc['H',i-1]:
                vmatrix.loc[state,i] = hmm_probability('H', state, seq[i]) * vmatrix.loc['H',i-1]
            else:
                vmatrix.loc[state,i] = hmm_probability('L', state, seq[i]) * vmatrix.loc['L',i-1]

    return max(vmatrix.loc['L',len(seq)-1],vmatrix.loc['H',len(seq)-1])


def log_viterbi_algorithm_matrix(seq):  # perform viterbi alg with log outputs to produce matrix
    states = ['H','L']
    vmatrix = pd.DataFrame(np.full((2,len(seq)), -np.inf), index=states, columns=range(len(seq)))
    for i in range(len(seq)):
        for state in states:
            if i == 0:
                vmatrix.loc[state,i] = np.log2(hmm_probability('init', state, seq[i]))
                continue
            if vmatrix.loc['L',i-1] < vmatrix.loc['H',i-1]:
                vmatrix.loc[state,i] = np.log2(hmm_probability('H', state, seq[i])) + vmatrix.loc['H',i-1]
            else:
                vmatrix.loc[state,i] = np.log2(hmm_probability('L', state, seq[i])) + vmatrix.loc['L',i-1]
    return vmatrix


if __name__ == "__main__":
    main()