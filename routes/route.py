import pickle
import os
from flask import jsonify, abort, Blueprint
from core.src import parse_sequences
import sys
import os

routes = Blueprint('routes', __name__)
current_directory = os.path.dirname(__file__)

@routes.errorhandler(400)
def resource_not_found(e):
    return jsonify(error=str(e)), 400

@routes.route('/<string:sequence>/<string:secondary_structure>/<int:score_threshold>', methods=['GET'])
def bayespairing(sequence, secondary_structure, score_threshold):
    try:
        arguments = {}

        # missing: verbose, pretrained, delta, samplesize, theta, lambda, w,s,aln,mod,inter,constraints,init,seq, mod
        threshold = (int) (score_threshold)
        if (threshold != -1):
            arguments["t"] = threshold
            
        sequence = str(sequence)
        # no ss should be ""
        secondary_structure = str(secondary_structure)
        dataset = "3dMotifAtlas_RELIABLE"
        # the default dataset is rna3dmotif

        # we load the modules from the dataset to get the number of modules available.
        #graphs = pickle.load(open("../models/" + dataset + "_one_of_each_graph.cPickle", "rb"))
        file_path = os.path.join(current_directory, "../core/models/" + dataset + "_one_of_each_graph.cPickle")
        if os.path.exists(file_path):
            while True:
                try:
                    nets = pickle.load(open(file_path, "rb"))
                    break
                except:
                    abort(400, description='Could not process dataset.')
            graphs = nets
        else:
            abort(400, description='The provided dataset does not appear to exist.')

        # needs to handle user input
        modules_to_check = range(len(graphs))

        # executes BayesPairing on the sequence
        print(parse_sequences.run_fasta(sequence, modules_to_check, dataset, secondary_structure, arguments)[0])

    except ValueError:
        abort(400, description='Received an invalid argument.')

    outputLocation = os.path.join(current_directory, "../core/output/output.pickle")
    output = pickle.load(open(outputLocation, "rb"))
    os.remove(outputLocation)

    return jsonify({'result': output})