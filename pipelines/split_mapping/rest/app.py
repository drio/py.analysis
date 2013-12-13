#!flask/bin/python

from tornado.wsgi import WSGIContainer
from tornado.httpserver import HTTPServer
from tornado.ioloop import IOLoop
from flask import Flask, jsonify, abort, make_response, request, redirect, url_for
from flask.ext.httpauth import HTTPBasicAuth
import json
import redis
import time

app = Flask(__name__)


""" All the interactions with redis are coded here.
    We serialize and deserialize the samples against
    the samples key in redis.
    There is also a key to keep the next id available.
    And the credentials for the webserver are stored
    here also.
"""
def setup_redis():
    e, g = app.redis.exists, app.redis.get
    if not e("user") or not e("pwd"):
        raise RuntimeError("Please set user and pwd in redis server.")
    f = open("static/credentials.json", "w")
    f.write("{'user': '%s', 'pwd': '%s'}" % (g("user"), g("pwd")))
    f.close()
    if not e("samples"):
        app.redis.set("samples", jsonify([]))
    if not e("current_id"):
        app.redis.set("current_id", 0)


def get_samples():
    return json.loads(app.redis.get("samples"))


def get_id():
    return app.redis.incr("last_id")


def save(samples):
    app.redis.set("samples", json.dumps(samples))


app.redis = redis.StrictRedis(host='localhost', port=6379, db=0)
app.setup_redis = setup_redis
app.get_samples = get_samples
app.get_id = get_id
app.save = save

"""
samples = [
    {
        'id': 1,
        'name': u'dummy',
        'steps': {}
    },
]
"""

auth = HTTPBasicAuth()


@auth.get_password
def get_password(username):
    if username == app.redis.get("user"):
        return app.redis.get("pwd")
    return None


@auth.error_handler
def unauthorized():
    return make_response(jsonify({'error': 'Unauthorized access'}), 401)


@app.errorhandler(404)
def not_found(error):
    return make_response(jsonify({'error': 'Not found'}), 404)


@app.route('/sapi/api/v1.0/samples', methods=['GET'])
@auth.login_required
def get_samples():
    return jsonify({'samples': app.get_samples()})


@app.route('/sapi/api/v1.0/samples/<int:sample_id>', methods=['GET'])
def get_sample(sample_id):
    sample = filter(lambda s: s['id'] == sample_id, app.get_samples())
    if len(sample) == 0:
        abort(404)
    return jsonify({'sample': sample[0]})


@app.route('/sapi/api/v1.0/samples', methods=['POST'])
def create_sample():
    if not request.json or not 'name' in request.json:
        abort(400)
    sample = {
        'id': app.get_id(),
        'name': request.json['name'],
        'steps': {},
    }
    samples = app.get_samples()
    samples.append(sample)
    app.save(samples)
    return jsonify({'sample': sample}), 201


@app.route('/sapi/api/v1.0/samples/<int:sample_id>', methods=['PUT'])
def update_sample(sample_id):
    samples = app.get_samples()
    sample = filter(lambda s: s['id'] == sample_id, samples)
    if len(sample) == 0:
        abort(404)
    if not request.json:
        abort(400)
    if 'name' in request.json and type(request.json['name']) != unicode:
        abort(400)
    if 'step' in request.json and type(request.json['step']) is not unicode:
        abort(400)
    sample[0]['name'] = request.json.get('name', sample[0]['name'])
    _step = request.json.get('step')
    sample[0]['steps'][_step] = time.time()
    app.save(samples)

    return jsonify({'sample': sample[0]})


@app.route('/sapi/api/v1.0/samples/<int:sample_id>', methods=['DELETE'])
def delete_sample(sample_id):
    samples = app.get_samples()
    sample = filter(lambda s: s['id'] == sample_id, samples)
    if len(sample) == 0:
        abort(404)
    samples.remove(sample[0])
    app.save(samples)
    return jsonify({'result': True})


@app.route('/frontend')
def home():
    return redirect(url_for('static', filename='frontend.html'))


if __name__ == '__main__':
    #app.run(debug=True)
    app.setup_redis()
    http_server = HTTPServer(WSGIContainer(app))
    http_server.listen(5000)
    IOLoop.instance().start()
    url_for('static', filename='frontend.html')
