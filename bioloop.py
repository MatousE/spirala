""""This is the Flask Python program that runs the HTML
for the bioloop.life website
"""
import os
from flask import Flask, render_template, request
import bioloopcalc as bio
app = Flask(__name__)


@app.route('/')
@app.route('/home')
def home():
    """"This function renders the home.html file
        and states that the websites route '/' and '/home'
        are both used for the home.html
    """
    return render_template('home.html')


@app.route('/SeqGen')
def seq_gen():
    """"This function renders the SeqGen.html file
        and states that the websites route '/SeqGen'
        is used for the SeqGen.html
    """
    return render_template('seqgen.html')


@app.route('/SeqGen', methods=['POST'])
def input_seq_gen():
    """"This function gets the user input from the
        SeqGen webpage and parses that input to the
        bio.sequencegen() so that it can return a
        random DNA sequence of the user specified length
    """
    seq_len = request.form['text']
    if len(seq_len) != 0:
        seq_len = int(seq_len)
        if isinstance(seq_len, int) and seq_len <= 10000000:
            new_seq = bio.sequencegen(seq_len)
            return render_template('output.html', value=new_seq)

    return render_template('output.html',
                           value='ERROR: INVALID INPUT MUST BE AN INT AND LESS THAN 10000000 BASES')


@app.route('/GCcontent')
def gc_content():
    """"This function renders the GCcontent.html file
        and states that the websites route '/GCcontent'
        is used for the GCcontent.html
    """
    return render_template('GCcontent.html')


@app.route('/GCcontent', methods=['POST'])
def input_gc_content():
    """"This function gets the user input from the
        GCcontent webpage and passes that input to
        bio.calculate() so that it can return a
        random the GC content for the user specified
        DNA sequence
    """
    seq = request.form['text']
    seq = seq.upper()
    if seq.count('A') > 0 or seq.count('C') > 0 or seq.count('G') > 0 or seq.count('T') > 0:
        gc_val = bio.calculatecg(seq)
        return render_template('output.html', value=gc_val)
    return render_template('output.html', value='ERROR: INVALID INPUT MUST CONTAIN A,C,G AND T')


@app.route('/GCcontent/fileinput')
def fileinput():
    """"Loading the fileinput.html template
    """
    return render_template('fileinput.html')


@app.route('/GCcontent/fileinput', methods=['POST'])
def get_cg_file():
    """"Getting file input"""
    f = request.files['file']
    f.save(f.filename)
    file = open(f.filename, 'r')

    seq = file.read()
    file.close()
    os.remove(f.filename)

    if seq.count('A') > 0 or seq.count('C') > 0 or seq.count('G') > 0 or seq.count('T') > 0:
        cg_val = bio.calculatecg(seq)
        return render_template('output.html', value=cg_val)
    return render_template('output.html', value='ERROR: INVALID INPUT MUST CONTAIN A,C,G AND T')


@app.route('/DNAmotif')
def dna_motif():
    """"This function renders the dnamorif.html file
        and states that the websites route '/DNAmotif'
        is used for the DNAmotif.html
    """
    return render_template('DNAmotif.html')


@app.route('/DNAmotif', methods=['POST'])
def input_motif():
    """Getting the motif and the FASTA file as well as
        outputting the locations in the FASTA file
    """
    file = request.files['file']
    motif = request.form.get('text')
    motif_occurrences = bio.dnamotifsearch(motif, file)
    return render_template('output.html', value=motif_occurrences)


if __name__ == '__main__':
    app.run(debug=True)
