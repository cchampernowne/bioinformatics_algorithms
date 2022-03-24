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
    fprob, hpost, lpost = forward_algorithm(seq)  # 'GGCA' for provided output files
    lvmatrix = log_viterbi_algorithm_matrix(seq)
    paths = most_prob_path(lvmatrix,seq)
    vprob = viterbi_algorithm(seq)
    print('\nRESULTS:\n\nPROBABILITY X:\n' + '{:.2e}'.format(fprob).replace('e-0', 'e-') + '\n\nVITERBI TABLE:')
    for path in paths:
        print(path + '\n')
    print('MOST PROBABLE PROBABILITY:\n' + '{:.2e}'.format(vprob).replace('e-0', 'e-'))


def server_results(seq):
    valid_sequence = re.compile(r'^ *(A|T|C|G|U|a|t|c|g|u)* *$')
    if re.fullmatch(valid_sequence, seq):
        seq = seq.strip()
        fprob, hpost, lpost = forward_algorithm(seq)  # 'GGCA' for provided output files
        lvmatrix = log_viterbi_algorithm_matrix(seq)
        paths = most_prob_path(lvmatrix,seq)
        vprob = viterbi_algorithm(seq)
        x_probability = '{:.2e}'.format(fprob).replace('e-0', 'e-')
        viterbi = print_matrix(lvmatrix,seq)
        most_prob_paths = ''
        for path in paths:
            most_prob_paths = most_prob_paths + path + '\n\n'
        most_prob_probability = '{:.2e}'.format(vprob).replace('e-0', 'e-')
        return [x_probability, viterbi, most_prob_paths, most_prob_probability]
    return False


def interperate_input(usr_input):
    re_quit = re.compile(r'^ *(Q|q) *$')
    re_valid_seq = re.compile(r'^ *(A|T|C|G|U|a|t|c|g|u)* *$')
    if re.fullmatch(re_quit, usr_input):
        return False
    if re.fullmatch(re_valid_seq, usr_input) and 4 < len(usr_input.strip()) < 201:
        return usr_input.strip()
    return interperate_input(input(ERROR_MSG + ',\nor enter Q to quit:'))


def print_matrix(matrix, seq):
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
    if vmatrix.loc['L',len(seq) - 1] == vmatrix.loc['H',len(seq) - 1]:
        for pos in vmatrix.columns:
            if vmatrix.loc['H',pos] > vmatrix.loc['L',pos]:
                mpp+='H'
            else:
                mpp+='L'
        mpp1 = mpp[:-1]
        mpp2 = mpp1 + 'H'
        mpp1+= 'L'

        return [mpp1,mpp2]

    for pos in vmatrix.columns:
        if vmatrix.loc['H',pos] > vmatrix.loc['L',pos]:
            mpp+='H'
        else:
            mpp+='L'
    return [mpp]


def hmm_probability(init_state, end_state, residue):  # hardcoded probabilities as given in assignment
    prob_score = 0
    if init_state == 'init':
        if end_state == 'H':
            if residue in ('A', 'a'):
                prob_score = 0.5*0.2
            if residue in ('C', 'c'):
                prob_score = 0.5*0.3
            if residue in ('G', 'g'):
                prob_score = 0.5*0.3
            if residue in ('T', 't'):
                prob_score = 0.5*0.2
        if end_state == 'L':
            if residue in ('A', 'a'):
                prob_score = 0.5*0.3
            if residue in ('C', 'c'):
                prob_score = 0.5*0.2
            if residue in ('G', 'g'):
                prob_score = 0.5*0.2
            if residue in ('T', 't'):
                prob_score = 0.5*0.3
    if init_state == 'H':
        if end_state == 'H':
            if residue in ('A', 'a'):
                prob_score = 0.5*0.2
            if residue in ('C', 'c'):
                prob_score = 0.5*0.3
            if residue in ('G', 'g'):
                prob_score = 0.5*0.3
            if residue in ('T', 't'):
                prob_score = 0.5*0.2
        if end_state == 'L':
            if residue in ('A', 'a'):
                prob_score = 0.5*0.3
            if residue in ('C', 'c'):
                prob_score = 0.5*0.2
            if residue in ('G', 'g'):
                prob_score = 0.5*0.2
            if residue in ('T', 't'):
                prob_score = 0.5*0.3
    if init_state == 'L':
        if end_state == 'H':
            if residue in ('A', 'a'):
                prob_score = 0.4*0.2
            if residue in ('C', 'c'):
                prob_score = 0.4*0.3
            if residue in ('G', 'g'):
                prob_score = 0.4*0.3
            if residue in ('T', 't'):
                prob_score = 0.4*0.2
        if end_state == 'L':
            if residue in ('A', 'a'):
                prob_score = 0.6*0.3
            if residue in ('C', 'c'):
                prob_score = 0.6*0.2
            if residue in ('G', 'g'):
                prob_score = 0.6*0.2
            if residue in ('T', 't'):
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


if __name__ == '__main__':
    main()