import numpy as np
import pandas as pd
import re


np.set_printoptions(threshold=np.inf)
ERROR_MSG = '\nERROR!\nInput not understood.\n - Sequence lengths must be between 5 - 200\n - Sequences must not have spaces, or any other characters between nucleotides\n - Nucleotides must be represented using only the characters A, T, C, G, U, a, t, c, g, and u.\n\nPlease try again,\n'
SERVER_ERROR_MSG = 'CHANGE LATER!! not proper FASTA format'

class Node():
    def __init__(self, name):
        self.name = name
        self.distance = -np.inf
        self.left = None
        self.right = None

    def add_left_child(self, left):
        self.left = left

    def add_right_child(self, right):
        self.right = right

    def has_children(self):
        if self.left is not None and self.right is not None:
            return True
        else: 
            return False

    def left_child_dist(self):
        total_dist = 0
        if self.has_children():
            total_dist += self.left.distance
            total_dist += self.left.left_child_dist()
        return total_dist
    
    def how_many_children(self):
        count = 0
        if self.left is not None:
            count = count + self.left.how_many_children()
        if self.right is not None:
            count = count + self.right.how_many_children()
        else:
            count = 1
        return count

    def children_list(self):
        children_list = []
        if self.left is not None:
            children_list.append(self.left)
        if self.right is not None:
            children_list.append(self.right)
        return children_list

    def tree_string(self,first):
        if first:
            tree_string = 'S' + self.name
        else:
            tree_string = 'S' + self.name + ':' + str(truncate(self.distance,1))
        if self.left is not None:
            tree_string += '(' + self.left.tree_string(False) + ')'
        if self.right is not None:
            tree_string += '(' + self.right.tree_string(False) + ')'
        return tree_string
            
    def deep_copy(self):
        copy = Node(self.name)
        copy.distance = self.distance
        if self.left is not None: 
            copy.left = self.left.deep_copy()
        if self.right is not None: 
            copy.right = self.right.deep_copy()
        return copy


def main():
    tree_list = []
    sequence_dict = get_sequences()
    if not sequence_dict:
        print('GOODBYE :)')
        return
    smatrix, node_dict = create_matrix(sequence_dict)
    upgma(smatrix, node_dict, tree_list)
    print('\nRESULTS:\nPAIRWISE DISTANCE MATRIX:\n' + print_matrix(smatrix,list(sequence_dict.keys())) + '\nPHYLOGENETIC TREE(S):')
    for tree in tree_list:
        print(tree.tree_string(True) + '\n')


def server_results(fasta_input):
    tree_list = []
    sequence_dict = fasta_to_sequences(fasta_input)
    if sequence_dict is not False:
        smatrix, node_dict = create_matrix(sequence_dict)
        upgma(smatrix, node_dict, tree_list)
        phylogenetic_trees = ''
        for tree in tree_list: phylogenetic_trees = tree.tree_string(True) + '\n\n'
        return [print_matrix(smatrix,list(sequence_dict.keys())), phylogenetic_trees]
    return False


def get_sequences():
    n = 1
    sequence_dict = {}
    usr_input = interperate_input(input('Enter the first sequence,\nor enter Q to quit:'), False)
    if not usr_input:
        return False
    sequence_dict[str(n)] = usr_input.strip()
    while True:
        n += 1
        if 10 < n:
            return sequence_dict
        if 3 < n:
            usr_input = interperate_input(input('Enter the next sequence,\nor enter Q to quit,\nor enter C to continue:'), True)
            if not usr_input:
                return False
            if isinstance(usr_input, str):
                sequence_dict[str(n)] = usr_input.strip()
                continue
            return sequence_dict
        else:
            usr_input = interperate_input(input('Enter the next sequence,\nor enter Q to quit:'), False)
            if not usr_input:
                return False
            if isinstance(usr_input, str):
                sequence_dict[str(n)] = usr_input.strip()
                continue


def interperate_input(usr_input, can_continue):
    re_quit = re.compile(r'^ *(Q|q) *$')
    re_continue = re.compile(r'^ *(C|c) *$')
    re_valid_seq = re.compile(r'^ *(A|T|C|G|U|a|t|c|g|u)* *$')
    if re.fullmatch(re_quit, usr_input):
        return False
    if can_continue and re.fullmatch(re_continue, usr_input):
        return True
    if re.fullmatch(re_valid_seq, usr_input) and 4 < len(usr_input.strip()) < 201:
        return usr_input
    if can_continue:
        return interperate_input(input(ERROR_MSG + 'or enter Q to quit,\n or enter C to continue:'), True)
    else:
        return interperate_input(input(ERROR_MSG + 'or enter Q to quit:'), False)


def truncate(num, digits):
    split = str(num).split('.', 1)
    try:
        split[1] = split[1][:digits]
    except:
        True
    return float(str(split[0]) + '.' + str(split[1]))


def fasta_to_sequences(fasta_input):  # retrieve sequences from file and return a list of all given sequences
    valid_sequence = re.compile(r'^(\r|\n| )*(>.*(\r|\n)(A|T|C|G|U|a|t|c|g|u)*(\r|\n| )*)*$')
    if re.fullmatch(valid_sequence, str(fasta_input)):
        raw_list = fasta_input.split('>')
        raw_list.pop(0)
        sequence_dict = {}
        n = 1
        for ele in raw_list:
            temp_list = ele.split('\n', 1)
            sequence_dict[str(n)] = temp_list[1].replace('\n', '').replace(' ', '')
            n = n + 1
        if not 2 < len(sequence_dict) < 11:
            return False
        for seq in sequence_dict.items():
            if 4 < len(seq) < 201:
                return False
        return sequence_dict
    return False


def min_diff(smatrix):  # get val and list of labels of all locations of min val
    min_val = np.inf
    min_label_list = []
    for index in smatrix.index:
        for column in smatrix.columns:
            if smatrix.loc[index,column] < min_val and smatrix.loc[index,column] > 0:
                min_val = smatrix.loc[index,column]
                min_label_list = [[index,column]]
            elif smatrix.loc[index, column] == min_val:
                if [column,index] not in min_label_list:
                    min_label_list.append([index,column])
    return min_val, min_label_list


def pos_diff(s1, s2):  # return the number of base positions in which s1 and s2 differ
    score = 0
    for i in range(max(len(s1), len(s2))):
        try:
            if s1[i] != s2[i]:
                score = score + 1
        except:
            score = score + 1
    return score


def print_matrix(smatrix, names):  # return the smatrix as a string for q 3.a
    total_matrix = '- '
    for name in names:
        total_matrix = total_matrix + 'S' + name + ' '
    total_matrix = total_matrix[:-1]
    total_matrix = total_matrix + '\n'

    for index in names:
        row_string = 'S' + index + ' '
        for column in names:
            row_string = row_string + str(int(smatrix.loc[index, column])) + ' '
        row_string = row_string[:-1]
        total_matrix = total_matrix + row_string + '\n'
    return total_matrix


def create_matrix(sequence_dict):  # create pairwise distance matrix of given sequences
    names = list(sequence_dict.keys())
    n = len(names)
    node_dict = {}

    for name in names:
        node_dict[name] = create_leaf(name, 0)

    smatrix = pd.DataFrame(np.full((n,n), -1), index=names, columns=names)

    for index in names:
        for column in names:
            smatrix.loc[index,column] = pos_diff(sequence_dict[index],sequence_dict[column])

    return smatrix, node_dict


def create_leaf(name, distance):
    node = Node(name)
    node.distance = distance
    return node


def join_nodes(a, b):
    ab = Node(a.name + b.name)
    ab.add_left_child(a)
    ab.add_right_child(b)
    return ab


def new_matrix(smatrix,node_dict):
    names = sorted(node_dict.keys(),key=len,reverse=True)
    m = len(node_dict)
    nmatrix = pd.DataFrame(np.full((m,m), 0), index=names, columns=names)

    for index in names:
        for column in names:
            if index == column:
                nmatrix.loc[index,column] = 0
                continue
            div = 0
            temp_val = 0
            if index not in smatrix.index:
                for child in node_dict[index].children_list():
                    if child.has_children():
                        div += child.how_many_children()
                        temp_val += smatrix.loc[child.name,column]*(child.how_many_children())
                    else:
                        div += 1
                        temp_val += smatrix.loc[child.name,column]
                
                nmatrix.loc[index,column] = temp_val / div
            elif column not in smatrix.columns:
                for child in node_dict[column].children_list():
                    if child.has_children():
                        div += child.how_many_children()
                        temp_val += smatrix.loc[index,child.name]*(child.how_many_children())
                    else:
                        div += 1
                        temp_val += smatrix.loc[index,child.name]

                nmatrix.loc[index,column] = temp_val / div
            else:
                nmatrix.loc[index,column] = smatrix.loc[index,column]
    return nmatrix


def dict_copy(node_dict):
    copy_dict = {}
    for key in node_dict.keys():
        copy_dict[key] = node_dict[key].deep_copy()
    return copy_dict


def upgma(smatrix, node_dict, tree_list):
    min_val, min_label_list = min_diff(smatrix)
    for index, column in min_label_list:
        copy_dict = dict_copy(node_dict)
        left_child = copy_dict[index]
        left_child.distance = min_val/2 - left_child.left_child_dist()
        del copy_dict[index]

        right_child = copy_dict[column]
        right_child.distance = min_val/2 - right_child.left_child_dist()
        del copy_dict[column]

        ab = join_nodes(left_child, right_child)
        copy_dict[ab.name] = ab
        nmatrix = new_matrix(smatrix,copy_dict)

        if len(copy_dict) <= 1:
            tree_list.append(list(copy_dict.values())[0])
            return
        upgma(nmatrix, copy_dict, tree_list)


if __name__ == '__main__':
    main()