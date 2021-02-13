from flask import jsonify, abort, Blueprint
from core.src import parse_sequences
import pickle
import os

@routes.errorhandler(400)
def resource_not_found(e):
    return jsonify(error=str(e)), 400

@app.route('/<int:year>/<int:month>/<title>')

@routes.route('/<string:sequence>/<string:secondary_structure>/<int:score_threshold>', methods=['GET'])
def one_adder(sequence, secondary_structure, score_threshold):
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
        graphs = parse_sequences.unpick(dataset,"models","one_of_each_graph.cPickle")
        # needs to handle user input
        modules_to_check = range(len(graphs))

        # executes BayesPairing on the sequence
        print(parse_sequences.run_fasta(sequence, modules_to_check, dataset, secondary_structure, arguments)[0])

    except ValueError:
        abort(400, description='Received an invalid argument.')

    output = pickle.load( open( "./core/output/output.pickle", "rb"))
    os.remove("./core/output/output.pickle")

    return jsonify({'result': output})