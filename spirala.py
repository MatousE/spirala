from flask import Flask, render_template

app = Flask(__name__)


@app.route('/')
@app.route('/home')
def index():
    return render_template('home.html')


@app.route('/SeqGen')
def SeqGen():
    return render_template('SeqGen.html')


if __name__ == '__main__':
    app.run(debug=True)
