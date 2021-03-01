import pickle
from flask import jsonify, abort, Blueprint, request, flash, request, redirect, url_for, current_app, send_from_directory
from werkzeug.utils import secure_filename
from core.src import parse_sequences
from utils import counter
from networkx.readwrite import json_graph
import networkx as nx
import json
import sys
import os
import threading
import atexit
routes = Blueprint("routes", __name__)
atexit.register(exit)
counter = counter.Counter()

CURRENT_DIRECTORY = os.path.dirname(__file__)
ALLOWED_EXTENSIONS = {"fa", "sl"}
DEBUG_MODE = False
BAYESPAIRING_VERBOSE = False

@routes.errorhandler(400)
def resource_not_found(e):
    return jsonify(error=str(e)), 400

@routes.route("/file", methods=["GET"])
def bayespairing_file():
    '''
    Represents the primary endpoint for calling BayesPairing with file input.

    This is primarily done for explicit distinction from the perspective of the client. 
    There is no functional difference between file and string input. Note that you must upload a file first via the upload_input endpoint.

    Inputs are provided in the form of query strings and are not explicitly defined to allow for flexibility.

    :returns: TODO
    '''
    return bayespairing(request.args)

@routes.route("/string", methods=["GET"])
def bayespairing_string():
    '''
    Represents the primary endpoint for calling BayesPairing with string input.

    This is primarily done for explicit distinction from the perspective of the client. 
    There is no functional difference between file and string input. Note that you must upload a file first via the upload_input endpoint.

    Inputs are provided in the form of query strings and are not explicitly defined to allow for flexibility.

    :returns: TODO
    '''
    return bayespairing(request.args)

@routes.route("/file_upload", methods=["POST"])
def upload_input():
    '''
    Uploads a sequence or alignment to be used as user input.

    Permitted formats include fasta and stockholm. 
    This must be done separately as a post to remain within spec (although there is no explicit rule barring doing this within a GET request).

    :returns: the name of the uploaded file so it can be passed as an argument to BayesPairing
    '''
    # check if the post request has the file part
    if "file" not in request.files:
        flash("No file part")
        abort(400, description="Did not receive a file.")
    file = request.files["file"]
    if file.filename == "":
        abort(400, description="Received an invalid file.")
        return redirect(request.url)
    if file and allowed_file(file.filename):
        count = counter.increment_and_get()
        filename = secure_filename(file.filename)
        extension = filename.rsplit(".", 1)[1]
        filename = str(count) + "." + extension
        file.save(os.path.join(current_app.config["UPLOAD_FOLDER"], filename))
        return jsonify({"filename": filename})
        
    return '''
    <!doctype html>
    <title>Upload new File</title>
    <h1>Upload new File</h1>
    <form method=post enctype=multipart/form-data>
      <input type=file name=file>
      <input type=submit value=Upload>
    </form>
    '''

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
            module_graph_mapping[module] = nx.to_dict_of_dicts(graph)
        
        return jsonify({"graphs": module_graph_mapping})
    except Exception as e:
        abort(400, "Received an invalid module: " + str(e))

def bayespairing(input):
    '''
    Processes input and runs BayesPairing.

    :param input: a dictionary containing a mapping of arguments to their value
    :returns: TODO
    '''
    try:
        arguments = {}

        sequence = input.get("sequence", default="", type = str)
        if (sequence == ""):
            abort(400, description="Did not receive a sequence as an argument.")
  
        secondary_structure = input.get("secondary_structure", default="", type = str)
        secondary_structure_infile = input.get("secondary_structure_infile", default = 0, type=int)

        # a secondary structure can be provided via a fasta file. This means a string cannot be provided, so verify that only one was received
        # moreover, if a fasta file is provided, a secondary structure should not be provided via string
        if (secondary_structure_infile and secondary_structure != ""):
            abort(400, description="Indicated that the secondary structure is provided in the fasta file, but also provided a separate structure.")
        elif (".fa" in sequence and secondary_structure != ""):
             abort(400, description="Received a secondary structure as a string, but a fasta file as a sequence. Please provide the single sequence as a string.")           
        elif (secondary_structure_infile):
            secondary_structure = "infile"

        # the default dataset is always used for now
        dataset = input.get("dataset", default="3dMotifAtlas_ALL", type = str)

        arguments["t"] = input.get("score_threshold", default=-2.3, type = float)

        delta = input.get("delta", default=-1, type = int)
        if (delta != -1):
            arguments["delta"] = delta

        arguments["samplesize"] = input.get("sample_size", default=20000, type = int)
        arguments["theta"] = input.get("theta", default=1, type = int)
        arguments["lambda"] = input.get("lambda", default=0.35, type = float)
        arguments["w"] = input.get("window_length", default=200, type = int)
        arguments["s"] = input.get("step_size", default=100, type = int)
        arguments["mod"] = input.get("modules", default="", type = str)
        arguments["constraints"] = input.get("constraints", default="", type = str)
        arguments["aln"] = input.get("alignment_provided", default=0, type = int)
        arguments["verbose"]=BAYESPAIRING_VERBOSE

        # we load the modules from the dataset to get the number of modules available.
        #graphs = pickle.load(open("../models/" + dataset + "_one_of_each_graph.cPickle", "rb"))
        file_path = os.path.join(CURRENT_DIRECTORY, "../core/models/" + dataset + "_one_of_each_graph.cPickle")
        if os.path.exists(file_path):
            while True:
                try:
                    nets = pickle.load(open(file_path, "rb"))
                    break
                except:
                    abort(400, description="Could not process dataset.")
            graphs = nets
        else:
            abort(400, description="The provided dataset does not appear to exist.")

        if (arguments["mod"] != ""):
            modules_to_check = [int(input_number) for input_number in arguments["mod"]]          
        else:
            modules_to_check = range(len(graphs))

        # executes BayesPairing on the sequence
        if (DEBUG_MODE):
            print(parse_sequences.run_fasta(sequence, modules_to_check, dataset, secondary_structure, arguments)[0])
        else:
            parse_sequences.run_fasta(sequence, modules_to_check, dataset, secondary_structure, arguments)

    except ValueError as e:
        abort(400, description="Received an invalid argument: " + str(e))
    except Exception as e:
        abort(400, description="BayesPairing failed to produce a result: " + str(e))

    # TODO, return output directly instead of having a file
    output_file = os.path.join(CURRENT_DIRECTORY, "../core/output/output.pickle")
    # check whether BayesPairing successfully produced output
    if not os.path.exists(output_file):
        abort(400, description="BayesPairing failed to produce output. Please ensure sure all input variables are valid.")

    output = pickle.load(open(output_file, "rb"))
    os.remove(output_file)
    
    input_path = os.path.join(CURRENT_DIRECTORY, "../core/input")
    if (".fa" in sequence or "st" in sequence):
        os.remove(os.path.join(input_path, sequence))


    return jsonify({"result": output})
    
def allowed_file(filename):
    '''
    Verifies that a file has the correct format.
    '''
    return "." in filename and \
           filename.rsplit(".", 1)[1].lower() in ALLOWED_EXTENSIONS

def exit():
    '''
    Handles shutdown in the event that the server terminates.
    '''
    counter.close_file()