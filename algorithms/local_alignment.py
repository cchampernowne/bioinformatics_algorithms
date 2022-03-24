import numpy as np
import re


np.set_printoptions(threshold=np.inf)
ERROR_MSG = 'One or both of your sequences is incorrectly formatted:\n - Sequence lengths must be between 5 - 50\n - Sequences must not have spaces, or any other characters between nucleotides\n - Nucleotides must be represented using only the characters A, T, C, G, U, a, t, c, g, and u'


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
        s1 = '-' + s1
        s2 = '-' + s2
        max_coords, max_score = max_coordinants(i,j,smatrix)
        for coord in max_coords:
            trace_paths(coord, s1, s2, seq1, seq2, smatrix, sequence_list)
        print('\nRESULTS:\n\nOPTIMAL SCORE:\n' + str(max_score) + '\n\nDYNAMIC PROGRAMMING MATRIX:\n' + print_matrix(smatrix) + '\n\nOPTIMAL ALLIGNMENTS:')
        for pair in sequence_list:
            print(pair[0] + '\n' + pair[1] + '\n')
    else:
        print('\nERROR!\n' + ERROR_MSG)


def server_results(s1, s2):
    if validity_check(s1, s2):
        sequence_list = []
        smatrix = create_matrix(s1, s2)
        seq1 = ''
        seq2 = ''
        j = len(s1)
        i = len(s2)
        s1 = '-' + s1
        s2 = '-' + s2
        max_coords, max_score = max_coordinants(i,j,smatrix)
        for coord in max_coords:
            trace_paths(coord, s1, s2, seq1, seq2, smatrix, sequence_list)
        optimal_alignments = ''
        for pair in sequence_list:
            optimal_alignments = optimal_alignments + pair[0] + '\n' + pair[1] + '\n\n'
        return [str(max_score), print_matrix(smatrix), optimal_alignments]
    return False


def validity_check(s1, s2):
    valid_sequence = re.compile(r'^ *(A|T|C|G|U|a|t|c|g|u)* *$')
    if re.fullmatch(valid_sequence, s1) and re.fullmatch(valid_sequence, s2) and 4 < len(s1.strip()) < 51 and 4 < len(s2.strip()) < 51:
        return True
    return False


def max_coordinants(i,j, smatrix):  # get the max score in the matrix, and the coordinates of everywhere in the matrix this max occurs
    max_coords = []
    max_val = 0
    for x in range(i + 1):
        for y in range(j + 1):
            if smatrix[x,y][0] > max_val:
                max_coords = [(x,y)]
                max_val = smatrix[x,y][0]
            elif smatrix[x,y][0] == max_val:
                max_coords.append((x,y))
    return max_coords, max_val


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


def create_matrix(s1,s2):  # convert given sequences to a matrix using the Smith-Waterman algorithm
    s1 = ' ' + s1
    s2 = ' ' + s2
    m = len(s1)
    n = len(s2)
    smatrix = np.array([[[-np.inf, []] for j in range(m)] for i in range(n)], dtype=object)  # create empty array of lists containing a -inf int and an empty list. The score of the cell, and the cells the score came from are held in each respectively
    for j in range(m):  # iterate through each cell of the matrix going row by row and starting at (0,0)
        for i in range(n):
            maxval = [[0, (-1,-1)]]  # initialize list of possible max values to include 0
            if i == 0 and j == 0:  # set pos (0,0) to a score of 0
                smatrix[0,0][0] = 0
                continue
            if s1[j] == s2[i]:  # if there is a match, top-left cell +2 is added to list of possible maximum scores for the cell
                maxval.append([smatrix[i-1,j-1][0] + 2, (i-1,j-1)])
            elif i-1 >=0 and j-1 >= 0:  # if there is no match, top-left cell -1 is added to list of possible maximum scores for the cell
                maxval.append([smatrix[i-1,j-1][0] - 1, (i-1,j-1)])
            if i-1 >=0 and j >= 0: # if the cell directly to the left exists, left cell - 2 is added to list of possible maximum scores for the cell
                maxval.append([smatrix[i-1,j][0] - 2, (i-1,j)])
            if i >=0 and j-1 >= 0:  # if the cell directly above exists, above cell - 2 is added to list of possible maximum scores for the cell
                maxval.append([smatrix[i,j-1][0] - 2, (i,j-1)])
            temp_np = np.array(maxval)
            for ele in temp_np:  # append coordinates that result in max score of cell
                if int(ele[0]) == max(temp_np[:,0]):
                    smatrix[i,j][1].append(ele[1])
            smatrix[i,j][0] = max(temp_np[:,0])  # assign max score to cell in smatrix
    return smatrix


def trace_paths(coord, s1, s2, seq1, seq2, smatrix, sequence_list):  # get all possible maximum scoring local alignments from matrix
    i, j = coord
    if smatrix[i,j][0] == 0:  # if the score in the cell is 0, this is the end of the local alignment
        sequence_list.append([seq1, seq2])
        return
    for tup in smatrix[i,j][1]:  # for each possible cell this score in the cell came from, call trace_paths again
        if tup == (i-1,j-1):
            trace_paths(tup, s1, s2, s1[j] + seq1, s2[i] + seq2, smatrix, sequence_list)
        if tup == (i-1,j):
            trace_paths(tup, s1, s2, '-' + seq1, s2[i] + seq2, smatrix, sequence_list)
        if tup == (i,j-1):
            trace_paths(tup, s1, s2, s1[j] + seq1, '-' + seq2, smatrix, sequence_list)


if __name__ == '__main__':
    main()