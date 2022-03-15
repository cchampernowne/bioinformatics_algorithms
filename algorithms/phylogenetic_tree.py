import numpy as np
import sys
import pandas as pd


np.set_printoptions(threshold=np.inf)
tree_list = []


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
    sequence_dict = get_sequences(sys.argv[1])
    smatrix, node_dict = create_matrix(sequence_dict)
    upgma(smatrix, node_dict)

    # pairwise distance matrix
    f = open("3.o1", "w+")
    f.write(print_matrix(smatrix,list(sequence_dict.keys())))
    f.close()

    # phylogenetic trees
    f = open("3.o4", "w+")
    f.write(str(len(tree_list)) + '\n')
    for tree in tree_list:
        f.write(tree.tree_string(True) + '\n\n')
    f.close()


def truncate(num, digits):
    split = str(num).split('.', 1)
    try:
        split[1] = split[1][:digits]
    except:
        True
    return float(str(split[0]) + '.' + str(split[1]))


def get_sequences(file_name):  # retrieve sequences from file and return a list of all given sequences
    f = open(file_name)
    raw_list = f.read().split('>')
    f.close()

    raw_list.pop(0)
    sequence_dict = {}

    for ele in raw_list:
        temp_list = ele.split('\n', 1)
        sequence_dict[temp_list[0].replace('S', '')] = temp_list[1].replace('\n', '')

    return sequence_dict


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


def upgma(smatrix,node_dict):
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
            global tree_list
            tree_list.append(list(copy_dict.values())[0])
            return
        upgma(nmatrix, copy_dict)


if __name__ == "__main__":
    main()