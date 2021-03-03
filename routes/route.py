import pickle
from flask import jsonify, abort, Blueprint, request, flash, request, redirect, url_for, current_app, send_from_directory
from werkzeug.utils import secure_filename
from werkzeug.datastructures import ImmutableMultiDict
from core.src import pipeline_bp as bp
import networkx as nx
import json
import sys
import os
import threading
import atexit
import tempfile

routes = Blueprint("routes", __name__)
#TODO: exit handling
#atexit.register(exit)

CURRENT_DIRECTORY = os.path.dirname(__file__)
ALLOWED_EXTENSIONS = {"fa", "fasta", "sl"}
DEBUG_MODE = False
BAYESPAIRING_VERBOSE = False

@routes.errorhandler(400)
def resource_not_found(e):
    return jsonify(error=str(e)), 400

@routes.route("/graph", methods=["GET"])
def get_graph_per_module():
    '''
    Given a list of modules, or a single module, returns graphs for each from the default dataset.

    :returns: a mapping of module IDs to their representative graphs
    '''
    try:
        modules = request.args.getlist("modules", type = int)
        print(modules)
        dataset_path = os.path.join(CURRENT_DIRECTORY, "../core/models/3dMotifAtlas_ALL_one_of_each_graph.cPickle")
        dataset = pickle.load(open(dataset_path, "rb"))
        module_graph_mapping = {}
        if (not modules):
            abort(400, description="Received an empty list.")
        
        for module in modules:
            # TODO, add master graph representation
            graph = dataset[module][0]
            graph_dict = {}
            graph_dict["nodes"] = list(graph.nodes(data=True))
            graph_dict["edges"] = list(graph.edges(data=True))
            module_graph_mapping[module] = graph_dict
        
        return jsonify({"graphs": module_graph_mapping})
    except Exception as e:
        abort(400, "Received an invalid module: " + str(e))

@routes.route("/string", methods=["POST"])
def bayespairing_string():
    '''
    Represents the primary endpoint for calling BayesPairing with string input.

    Inputs are provided in the form of query strings and are not explicitly defined to allow for flexibility.

    :returns: a jsonified dictionary of modules and their hits
    '''
    if "arguments" not in request.form:
        abort(400, description="Did not receive any arguments.")
    return bayespairing(eval(request.form.get("arguments")))

@routes.route("/file", methods=["POST"])
def bayespairing_file():
    '''
    Represents the primary endpoint for calling BayesPairing with file input.

    Inputs are provided in the request body. The file is provided as an individual argument with key input, 
    with arguments provided as a JSON dictionary with key "arguments". Arguments are the same as in the string input case.

    :returns: a jsonified dictionary of modules and their hits
    '''
    if "input" not in request.files:
        abort(400, description="Did not receive an input file.")
    if "arguments" not in request.form:
        abort(400, description="Did not receive any arguments.")

    file = request.files["input"]

    if file.filename == "":
        abort(400, description="Received an invalid file.")
    if file and allowed_file(file.filename):
        # we store the file as a temporary file
        # this has inherent benefits, it handles uniqueness of naming and also destroys the file once we are done
        file.seek(0)
        temp = tempfile.NamedTemporaryFile()
        temp.write(file.read())
        temp.seek(0)

        if (DEBUG_MODE):
            print("BayesPairing File: Stored temporary file: %s\n" %(temp.name))

        # for consistency, we convert from a string to an ImmutableMultiDict (and a dictionary in-between) 
        # so we can reuse the BaysesPairing method
        return bayespairing(eval(request.form.get("arguments")), file.filename.rsplit(".", 1)[1], temp)
        
    abort(400, description="BayesPairing failed for an unknown reason. Please check your inputs.")

def bayespairing(input, input_file_type = None, input_file = None):
    '''
    Processes input and runs BayesPairing.

    :param input: a dictionary (ImmutableMultiDict) containing a mapping of arguments to their value
    :param input_file_type: the type of file if one is provided (fasta or stockholm)
    :param input_file: a path to a stockholm or fasta file if no sequence string is provided
    :returns: a jsonified dictionary of modules and their hits
    '''
    try:
        arguments = {}

        sequence = input.get("sequence")
        if (not sequence and not input_file):
            abort(400, description="Did not receive a sequence as an argument or a file.")
        elif (not sequence):
            sequence = input_file.name
  
        secondary_structure = input.get("secondary_structure", "")
        secondary_structure_infile = input.get("secondary_structure_infile", 0)

        # a secondary structure can be provided via a fasta file. This means a string cannot be provided, so verify that only one was received
        # moreover, if a fasta file is provided, a secondary structure should not be provided via string
        if (secondary_structure_infile and secondary_structure != ""):
            if (input_file):
                input_file.close()
            abort(400, description="Indicated that the secondary structure is provided in the fasta file, but also provided a separate structure.")
        elif ((input_file_type == "fa" or input_file_type == "fasta") and secondary_structure != ""):
            input_file.close()
            abort(400, description="Received a secondary structure as a string, but a fasta file as a sequence. Please provide the single sequence as a string.")           
        elif (secondary_structure_infile):
            secondary_structure = "infile"

        # the default dataset is always used for now
        dataset = input.get("dataset", "3dMotifAtlas_ALL")

        # load all arguments if they were provided or set the default value
        arguments["t"] = input.get("score_threshold", -2.3)
        arguments["samplesize"] = input.get("sample_size", 20000)
        arguments["theta"] = input.get("theta", 1)
        arguments["lambda"] = input.get("lambda", 0.35)
        arguments["w"] = input.get("window_length", 200)
        arguments["s"] = input.get("step_size", 100)
        arguments["mod"] = input.get("modules", "")
        arguments["constraints"] = input.get("constraints", "")
        arguments["aln"] = input.get("alignment_provided", 0)
        arguments["verbose"]=BAYESPAIRING_VERBOSE

        # we load the modules from the dataset to get the number of modules available.
        file_path = os.path.join(CURRENT_DIRECTORY, "../core/models/" + dataset + "_one_of_each_graph.cPickle")
        if os.path.exists(file_path):
            while True:
                try:
                    nets = pickle.load(open(file_path, "rb"))
                    break
                except Exception as e:
                    if (input_file):
                        input_file.close()
                    abort(400, description="Could not process dataset due to error: " + str(e))
            graphs = nets
        else:
            if (input_file):
                input_file.close()
            abort(400, description="The provided dataset does not appear to exist.")

        if (arguments["mod"] != ""):
            modules_to_check = [int(input_number) for input_number in arguments["mod"]]          
        else:
            modules_to_check = range(len(graphs))

        # executes BayesPairing on the sequence
        output = bp.run_fasta(sequence, modules_to_check, dataset, secondary_structure, arguments, input_file_type)

        if (input_file):
            input_file.close()
        if not output:
            abort(400, description="BayesPairing failed to produce output. Please ensure sure all input variables are valid. If a fasta or stockholm file were provided, they were preserved.")           

        return jsonify({"result": output})

    except ValueError as e:
        if (input_file):
            input_file.close()
        abort(400, description="Received an invalid argument: " + str(e))
    except Exception as e:
        if (input_file):
            input_file.close()
        abort(400, description="BayesPairing failed to produce a result: " + str(e))    

def allowed_file(filename):
    '''
    Verifies that a file has the correct format.

    :param filename: the name of the file to check
    :returns: true if the file format is valid
    '''
    return "." in filename and \
           filename.rsplit(".", 1)[1].lower() in ALLOWED_EXTENSIONS
