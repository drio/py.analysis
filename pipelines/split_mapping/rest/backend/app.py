#!flask/bin/python

from flask import Flask, jsonify, abort, make_response, request
import time

app = Flask(__name__)

samples = [
    {
        'id': 1,
        'name': u'dummy',
        'steps': {}
    },
]


@app.errorhandler(404)
def not_found(error):
    return make_response(jsonify({'error': 'Not found'}), 404)


@app.route('/todo/api/v1.0/samples', methods=['GET'])
def get_samples():
    return jsonify({'samples': samples})


@app.route('/todo/api/v1.0/samples/<int:sample_id>', methods=['GET'])
def get_sample(sample_id):
    sample = filter(lambda s: s['id'] == sample_id, samples)
    if len(sample) == 0:
        abort(404)
    return jsonify({'sample': sample[0]})


@app.route('/todo/api/v1.0/samples', methods=['POST'])
def create_sample():
    if not request.json or not 'name' in request.json:
        abort(400)
    sample = {
        'id': samples[-1]['id'] + 1,
        'name': request.json['name'],
        'steps': {},
    }
    samples.append(sample)
    return jsonify({'sample': sample}), 201


@app.route('/todo/api/v1.0/samples/<int:sample_id>', methods=['PUT'])
def update_sample(sample_id):
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

    return jsonify({'sample': sample[0]})


@app.route('/todo/api/v1.0/samples/<int:sample_id>', methods=['DELETE'])
def delete_sample(sample_id):
    sample = filter(lambda s: s['id'] == sample_id, samples)
    if len(sample) == 0:
        abort(404)
    samples.remove(sample[0])
    return jsonify({'result': True})


if __name__ == '__main__':
    app.run(debug=True)
