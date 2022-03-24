from flask import Flask, request, render_template
from algorithms.local_alignment import server_results as local_results, ERROR_MSG as ERROR_ALIGNMENT
from algorithms.global_alignment import server_results as global_results
from algorithms.phylogenetic_tree import server_results as upgma_results, SERVER_ERROR_MSG as ERROR_UPGMA
from algorithms.hmm import server_results as hidden_mm_results, ERROR_MSG as ERROR_HMM
from algorithms.nussinov import server_results as nuss_results


app = Flask(__name__)


@app.route('/')
def home():
    return render_template('home.html')


@app.route('/smith-waterman')
def smith_waterman():
    return render_template('smith_waterman.html')


@app.route('/smith-waterman/results', methods=['POST'])
def smith_results():
    results = local_results(request.form['sequence1'], request.form['sequence2'])
    if results == False:
        return render_template('smith_waterman_error.html', sequence1 = request.form['sequence1'], sequence2 = request.form['sequence2'], error_msg = ERROR_ALIGNMENT)
    return render_template('smith_waterman_results.html', sequence1 = request.form['sequence1'], sequence2 = request.form['sequence2'], optimal_score = results[0], dynamic_matrix = results[1], optimal_alignments = results[2])


@app.route('/needleman-wunsch')
def needleman_wunsch():
    return render_template('needleman_wunsch.html')


@app.route('/needleman-wunsch/results', methods=['POST'])
def needleman_results():
    results = global_results(request.form['sequence1'], request.form['sequence2'])
    if results == False:
        return render_template('needleman_wunsch_error.html', sequence1 = request.form['sequence1'], sequence2 = request.form['sequence2'], error_msg = ERROR_ALIGNMENT)
    return render_template('needleman_wunsch_results.html', sequence1 = request.form['sequence1'], sequence2 = request.form['sequence2'], optimal_score = results[0], dynamic_matrix = results[1], optimal_alignments = results[2])


@app.route('/upgma')
def phylogenetic_tree():
    return render_template('phylogenetic_tree.html')


@app.route('/upgma/results', methods=['POST'])
def phylogenetic_tree_results():
    results = upgma_results(request.form['fasta_input'])
    if results == False:
        return render_template('phylogenetic_tree_error.html', fasta_input=request.form['fasta_input'], error_msg = ERROR_UPGMA)
    return render_template('phylogenetic_tree_results.html', fasta_input=request.form['fasta_input'], distance_matrix = results[0], phylogenetic_trees = results[1])


@app.route('/hmm')
def hmm():
    return render_template('hmm.html')


@app.route('/hmm/results', methods=['POST'])
def hmm_results():
    results = hidden_mm_results(request.form['sequence_input'])
    if results == False:
        return render_template('hmm_error.html', sequence_input = request.form['sequence_input'], error_msg = ERROR_HMM)
    return render_template('hmm_results.html', sequence_input = request.form['sequence_input'], x_probability = results[0], viterbi = results[1], most_prob_paths = results[2], most_prob_probability = results[3])


@app.route('/nussinov')
def nussinov():
    return render_template('nussinov.html')


@app.route('/nussinov/results', methods=['POST'])
def nussinov_results():
    results = nuss_results(request.form['sequence_input'])
    if results == False:
        return render_template('nussinov_error.html', sequence_input = request.form['sequence_input'])
    return render_template('nussinov_results.html', sequence_input = request.form['sequence_input'], structures_ranked = results)


if __name__ == '__main__':
    app.run(debug=True)