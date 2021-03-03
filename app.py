from flask import Flask, jsonify
from routes import route

app = Flask(__name__)
app.register_blueprint(route.routes)

@app.route('/', methods=['GET'])
def hello():
    return jsonify(about='Hello from Bayespairing, WP!')

# it appears that it is unnecessary to specify the port/localhost
if __name__ == '__main__':
    app.run(host="localhost", port=5001, debug=True)