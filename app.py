from flask import Flask, jsonify
from routes import route
import os

current_directory = os.path.dirname(__file__)
UPLOAD_FOLDER = os.path.join(current_directory, 'core/input')

app = Flask(__name__)
app.register_blueprint(route.routes)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

@app.route('/', methods=['GET'])
def hello():
    return jsonify(about='Hello from Bayespairing, WP!')

# it appears that it is unnecessary to specify the port/localhost
if __name__ == '__main__':
    app.run(host="localhost", port=5001, debug=True)