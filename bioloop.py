""""This is the Flask Python program that runs the HTML
for the bioloop.life website
"""
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
    seq_len = int(seq_len)
    new_seq = bio.sequencegen(seq_len)
    return new_seq


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
    upper_seq = seq.upper()
    gc_val = bio.calculatecg(upper_seq)
    return gc_val


if __name__ == '__main__':
    app.run(debug=True)
