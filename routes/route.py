import pickle
import os
from flask import jsonify, abort, Blueprint, request
from core.src import parse_sequences
import sys
import os

routes = Blueprint('routes', __name__)
current_directory = os.path.dirname(__file__)

@routes.errorhandler(400)
def resource_not_found(e):
    return jsonify(error=str(e)), 400

@routes.route('/', methods=['GET'])
def bayespairing():
    try:
        arguments = {}

        # mandatory args
        sequence = request.args.get("sequence", default="", type = str)
        if (sequence == ""):
            abort(400, description='Did not receive a sequence as an argument.')
  
        secondary_structure = request.args.get("secondary_structure", default="", type = str)
        secondary_structure_infile = request.args.get("secondary_structure_infile", default = 0, type=int)

        # a secondary structure can be provided via a fasta file. This means a string cannot be provided, so verify that only one was received
        if (secondary_structure_infile and secondary_structure != ""):
            abort(400, description='Indicated that the secondary structure is provided in the fasta file, but also provided a separate structure.')
        elif (secondary_structure_infile):
            secondary_structure = "infile"

        # the default dataset is always used for now
        dataset = request.args.get("dataset", default="3dMotifAtlas_RELIABLE", type = str)

        arguments["t"] = request.args.get("score_threshold", default=-2.3, type = float)
        arguments["delta"] = request.args.get("delta", default=-1, type = int)
        arguments["samplesize"] = request.args.get("sample_size", default=100, type = int)
        arguments["theta"] = request.args.get("theta", default=1, type = int)
        arguments["lambda"] = request.args.get("lambda", default=0.2, type = float)
        arguments["w"] = request.args.get("window_length", default=200, type = int)
        arguments["s"] = request.args.get("step_size", default=100, type = int)
        arguments["mod"] = request.args.get("modules", default="", type = str)
        arguments["constraints"] = request.args.get("constraints", default="", type = str)
        arguments["aln"] = request.args.get("alignment_provided", default=0, type = int)

        # we load the modules from the dataset to get the number of modules available.
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

        if (arguments["mod"] != ""):
            modules_to_check = [int(input_number) for input_number in arguments["mod"]]          
        else:
            modules_to_check = range(len(graphs))

        # executes BayesPairing on the sequence
        parse_sequences.run_fasta(sequence, modules_to_check, dataset, secondary_structure, arguments)

    except ValueError:
        abort(400, description='Received an invalid argument.')

    outputLocation = os.path.join(current_directory, "../core/output/output.pickle")
    # check whether BayesPairing successfully produced output
    if not os.path.exists(outputLocation):
        abort(400, description='BayesPairing failed to produce output. Please make sure all input variables are valid.')

    output = pickle.load(open(outputLocation, "rb"))
    os.remove(outputLocation)

    return jsonify({'result': output})