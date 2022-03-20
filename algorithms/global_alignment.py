import numpy as np
import re
from algorithms.local_alignment import ERROR_MSG


np.set_printoptions(threshold=np.inf)


def main():
    s1 = input('Enter the first sequence: ')
    s2 = input('Enter the second sequence: ')
    if validity_check(s1, s2):
        smatrix = create_matrix(s1, s2)
        sequence_list = []
        seq1 = ''
        seq2 = ''
        j = len(s1)
        i = len(s2)
        max_score = smatrix[i,j][0]
        trace_paths(i, j, s1, s2, seq1, seq2, smatrix, sequence_list)
        print('\nRESULTS:\n')
        print('OPTIMAL SCORE:')
        print(str(max_score) + '\n')
        print('DYNAMIC PROGRAMMING MATRIX:')
        print(print_matrix(smatrix) + '\n')
        print('OPTIMAL ALLIGNMENTS:')
        for pair in sequence_list:
            print(pair[0] + '\n' + pair[1] + '\n')
    else:
        print("\nERROR!\n" + ERROR_MSG)


def server_results(s1, s2):
    if validity_check(s1, s2):
        sequence_list = []
        smatrix = create_matrix(s1, s2)
        seq1 = ''
        seq2 = ''
        j = len(s1)
        i = len(s2)
        max_score = smatrix[i,j][0]
        trace_paths(i, j, s1, s2, seq1, seq2, smatrix, sequence_list)
        optimal_alignments = ""
        for pair in sequence_list:
            optimal_alignments = optimal_alignments + pair[0] + "\n" + pair[1] + "\n\n"
        return [str(max_score), print_matrix(smatrix), optimal_alignments]
    return False


def validity_check(s1, s2):
    valid_sequence = re.compile(r"^(A|T|C|G|U|a|t|c|g|u)*$")
    if re.fullmatch(valid_sequence, s1) and re.fullmatch(valid_sequence, s2) and 4 < len(s1.strip()) < 51 and 4 < len(s2.strip()) < 51:
        return True
    return False


def print_matrix(smatrix):  # return the values of each cell in the matrix as a string
    i,j,z = smatrix.shape
    total_matrix = ''
    for row in range(j):
        row_string = ''
        for ele in smatrix[:,row]:
            row_string = row_string + str(ele[0]) + ' '
        row_string = row_string[:-1]
        total_matrix = total_matrix + row_string + '\n'
    return total_matrix


def create_matrix(s1,s2):  # convert given sequences to a matrix using the Needleman-Wunsch algorithm
    s1 = ' ' + s1
    s2 = ' ' + s2
    m = len(s1)
    n = len(s2)
    smatrix = np.array([[[-np.inf, []] for j in range(m)] for i in range(n)], dtype=object)
    for j in range(m):  # iterate through each cell of the matrix going row by row and starting at (0,0)
        for i in range(n):
            maxval = []
            if i == 0 and j == 0:  # set pos (0,0) to a score of 0
                smatrix[0,0][0] = 0
                continue
            if s1[j] == s2[i]:  # if there is a match, top-left cell +2 is added to list of possible maximum scores for the cell
                maxval.append([smatrix[i-1,j-1][0] + 2, (i-1,j-1)])
            elif i-1 >=0 and j-1 >= 0:  # if there is no match, top-left cell -1 is added to list of possible maximum scores for the cell
                maxval.append([smatrix[i-1,j-1][0] - 1, (i-1,j-1)])
            if i-1 >=0 and j >= 0:  # if the cell directly to the left exists, left cell - 2 is added to list of possible maximum scores for the cell
                maxval.append([smatrix[i-1,j][0] - 2, (i-1,j)])
            if i >=0 and j-1 >= 0:  # if the cell directly above exists, above cell - 2 is added to list of possible maximum scores for the cell
                maxval.append([smatrix[i,j-1][0] - 2, (i,j-1)])
            temp_np = np.array(maxval)
            for ele in temp_np:  # append coordinates that result in max score of cell
                if int(ele[0]) == max(temp_np[:,0]):
                    smatrix[i,j][1].append(ele[1])
            smatrix[i,j][0] = max(temp_np[:,0])  # assign max score to cell in smatrix
    return smatrix


def trace_paths(i, j, s1, s2, seq1, seq2, smatrix, sequence_list):  # get all possible maximum scoring local alignments from matrix
    if not smatrix[i,j][1]:  # if the cell coordinates are (0,0), this is the end of the local alignment
        sequence_list.append([seq1, seq2])
        return
    for tup in smatrix[i,j][1]:   # for each possible cell this score in the cell came from, call trace_paths again
        if tup == (i-1,j-1):
            trace_paths(i-1, j-1, s1[:-1], s2[:-1], s1[-1] + seq1, s2[-1] + seq2, smatrix, sequence_list)
        if tup == (i-1,j):
            trace_paths(i-1, j, s1, s2[:-1], '-' + seq1, s2[-1] + seq2, smatrix, sequence_list)
        if tup == (i,j-1):
            trace_paths(i, j-1, s1[:-1], s2, s1[-1] + seq1, '-' + seq2, smatrix, sequence_list)


if __name__ == "__main__":
    main()