from flask import Flask, request, render_template
import logging
import algorithms.local_alignment
import algorithms.global_alignment
import algorithms.phylogenetic_tree
import algorithms.hmm
logging.basicConfig(level=logging.INFO)

app = Flask(__name__)


@app.route('/')
def home():
    return render_template('home.html')


@app.route('/smith-waterman')
def smith_waterman():
    return render_template('smith_waterman.html')


@app.route('/smith-waterman/results', methods=['POST'])
def smith_results():
    s1, s2 = request.form['sequence1'], request.form['sequence2']
    smatrix = algorithms.local_alignment.create_matrix(s1, s2)

    seq1 = ''
    seq2 = ''
    j = len(s1)
    i = len(s2)
    s1 = '-' + s1
    s2 = '-' + s2

    max_coords, max_score = algorithms.local_alignment.max_coordinants(i,j,smatrix)
    for coord in max_coords:
        algorithms.local_alignment.trace_paths(coord, s1, s2, seq1, seq2, smatrix)
    optimal_alignments = ''
    for pair in algorithms.local_alignment.sequence_list:
        optimal_alignments = optimal_alignments + pair[0] + "\n" + pair[1] + "\n\n"
    return render_template('smith_waterman_results.html', sequence1 = request.form['sequence1'], sequence2 = request.form['sequence2'], optimal_score = str(max_score), dynamic_matrix = algorithms.local_alignment.print_matrix(smatrix), optimal_alignments = optimal_alignments)


@app.route('/needleman-wunsch')
def needleman_wunsch():
    return render_template('needleman_wunsch.html')


@app.route('/needleman-wunsch/results', methods=['POST'])
def needleman_results():
    s1, s2 = request.form['sequence1'], request.form['sequence2']
    smatrix = algorithms.global_alignment.create_matrix(s1, s2)
    seq1 = ''
    seq2 = ''
    j = len(s1)
    i = len(s2)
    max_score = smatrix[i,j][0]
    algorithms.global_alignment.trace_paths(i, j, s1, s2, seq1, seq2, smatrix)
    optimal_alignments = ""
    for pair in algorithms.global_alignment.sequence_list:
        optimal_alignments = optimal_alignments + pair[0] + "\n" + pair[1] + "\n\n"
    return render_template('needleman_wunsch_results.html', sequence1 = request.form['sequence1'], sequence2 = request.form['sequence2'], optimal_score = str(max_score), dynamic_matrix = algorithms.global_alignment.print_matrix(smatrix), optimal_alignments = optimal_alignments)


@app.route('/phylogenetic-tree')
def phylogenetic_tree():
    return render_template('phylogenetic_tree.html')


@app.route('/phylogenetic-tree/results', methods=['POST'])
def phylogenetic_tree_results():
    sequence_dict = algorithms.phylogenetic_tree.get_sequences(request.form["sequence_file"])
    smatrix, node_dict = algorithms.phylogenetic_tree.create_matrix(sequence_dict)
    algorithms.phylogenetic_tree.upgma(smatrix, node_dict)
    phylogenetic_trees = ""
    for tree in algorithms.phylogenetic_tree.tree_list:
        phylogenetic_trees = tree.tree_string(True) + "\n\n"
    return render_template('phylogenetic_tree_results.html', distance_matrix = algorithms.phylogenetic_tree.print_matrix(smatrix,list(sequence_dict.keys())), phylogenetic_trees = phylogenetic_trees)


@app.route('/hmm')
def hmm():
    return render_template('hmm.html')


@app.route('/hmm/results', methods=['POST'])
def hmm_results():
    seq = request.form["sequence_input"]
    fprob, hpost, lpost = algorithms.hmm.forward_algorithm(seq)  # 'GGCA' for provided output files
    lvmatrix = algorithms.hmm.log_viterbi_algorithm_matrix(seq)
    paths, mul = algorithms.hmm.most_prob_path(lvmatrix,seq)
    vprob = algorithms.hmm.viterbi_algorithm(seq)

    x_probability = "{:.2e}".format(fprob).replace("e-0", "e-")
    viterbi = algorithms.hmm.print_matrix(lvmatrix,seq)
    most_prob_paths = ""
    for path in paths:
        most_prob_paths = most_prob_paths + path + '\n\n'
    most_prob_probability = "{:.2e}".format(vprob).replace("e-0", "e-")
    return render_template('hmm_results.html', x_probability = x_probability, viterbi = viterbi, most_prob_paths = most_prob_paths, most_prob_probability = most_prob_probability)


@app.route('/nussinov')
def nussinov():
    return render_template('nussinov.html')


@app.route('/nussinov/results', methods=['POST'])
def nussinov_results():
    return render_template('nussinov_results.html')


# 
# 
# @app.route('/pick-top')
# def pick_top():
#     global middle
#     if middle.top:  # if a top sentence exists, set it as current middle sentence
#         middle = middle.top
#     top, left, right, bottom = curr_choices()
#     return render_template('template.html', middle=middle.sentence, top=top, left=left, right=right, bottom=bottom)
# 
# 
# @app.route('/pick-left')
# def pick_left():
#     global middle
#     if middle.left:  # if a left sentence exists, set it as current middle sentence
#         middle = middle.left
#     top, left, right, bottom = curr_choices()
#     return render_template('template.html', middle=middle.sentence, top=top, left=left, right=right, bottom=bottom)
# 
# 
# @app.route('/pick-right')
# def pick_right():
#     global middle
#     if middle.right:  # if a right sentence exists, set it as current middle sentence
#         middle = middle.right
#     top, left, right, bottom = curr_choices()
#     return render_template('template.html', middle=middle.sentence, top=top, left=left, right=right, bottom=bottom)
# 
# 
# @app.route('/pick-bottom')
# def pick_bottom():
#     global middle
#     if middle.bottom:  # if a bottom sentence exists, set it as current middle sentence
#         middle = middle.bottom
#     top, left, right, bottom = curr_choices()
#     return render_template('template.html', middle=middle.sentence, top=top, left=left, right=right, bottom=bottom)
# 
# 
# @app.route('/top', methods=['POST'])
# def top():  # add the posted top sentence to the current middle node's top node
#     global middle
#     middle.add_top(Sentence_Node(x))
#     top, left, right, bottom = curr_choices()
#     return render_template('template.html', middle=middle.sentence, top=top, left=left, right=right, bottom=bottom)
# 
# 
# @app.route('/left', methods=['POST'])
# def left():  # add the posted left sentence to the current middle node's left node
#     global middle
#     middle.add_left(Sentence_Node(request.form['left']))
#     top, left, right, bottom = curr_choices()
#     return render_template('template.html', middle=middle.sentence, top=top, left=left, right=right, bottom=bottom)
# 
# 
# @app.route('/right', methods=['POST'])
# def right():  # add the posted right sentence to the current middle node's right node
#     global middle
#     middle.add_right(Sentence_Node(request.form['right']))
#     top, left, right, bottom = curr_choices()
#     return render_template('template.html', middle=middle.sentence, top=top, left=left, right=right, bottom=bottom)
# 
# 
# @app.route('/bottom', methods=['POST'])
# def bottom():  # add the posted bottom sentence to the current middle node's bottom node
#     global middle
#     middle.add_bottom(Sentence_Node(request.form['bottom']))
#     top, left, right, bottom = curr_choices()
#     return render_template('template.html', middle=middle.sentence, top=top, left=left, right=right, bottom=bottom)
# 
# 
# @app.route('/top')
# def top_reload():  # triggers on reload of localhost:5000/top
#     global middle
#     top, left, right, bottom = curr_choices()
#     return render_template('template.html', middle=middle.sentence, top=top, left=left, right=right, bottom=bottom)
# 
# 
# @app.route('/left')
# def left_reload():  # triggers on reload of localhost:5000/left
#     global middle
#     top, left, right, bottom = curr_choices()
#     return render_template('template.html', middle=middle.sentence, top=top, left=left, right=right, bottom=bottom)
# 
# 
# @app.route('/right')
# def right_reload():  # triggers on reload of localhost:5000/right
#     global middle
#     top, left, right, bottom = curr_choices()
#     return render_template('template.html', middle=middle.sentence, top=top, left=left, right=right, bottom=bottom)
# 
# 
# @app.route('/bottom')
# def bottom_reload():  # triggers on reload of localhost:5000/bottom
#     global middle
#     top, left, right, bottom = curr_choices()
#     return render_template('template.html', middle=middle.sentence, top=top, left=left, right=right, bottom=bottom)


if __name__ == "__main__":
    app.run(debug=True)