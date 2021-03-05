import pickle
from flask import jsonify, abort, Blueprint, request, flash, request, current_app
from werkzeug.utils import secure_filename
from core.src.pipeline import pipeline_bp as bp
from core.src.pipeline import pipeline_chefschoice as chefschoice
import sys
import os
import threading
import atexit
import tempfile
import base64
import json
routes = Blueprint("routes", __name__)
# TODO: exit handling
# atexit.register(exit)

CURRENT_DIRECTORY = os.path.dirname(__file__)
ALLOWED_EXTENSIONS = {"fa", "fasta", "sl"}
DEBUG_MODE = False
BAYESPAIRING_VERBOSE = False


@routes.errorhandler(400)
def resource_not_found(e):
    return jsonify(error=str(e)), 400


@routes.route("/graphs", methods=["POST"])
def get_graph_per_module():
    '''
    Given a list of modules, or a single module, returns graphs for each from the default dataset.

    This is a POST because of the possibility of uploading a dataset.
    (TODO issue #5 Implement support for uploading a user specified dataset)

    :returns: a mapping of module IDs to their representative graphs.
    '''
    try:
        modules = request.form.get("modules")
        if (not modules):
            abort(400, "Did not receive any modules to fetch graphs for.")
        modules = eval(modules)

        if (request.files.get("dataset")):
            flash(
                "User provided datasets are not yet supported. The default will be used instead.")

        module_graph_mapping = retrieve_graphs(modules)

        return jsonify({"graphs": module_graph_mapping})
    except Exception as e:
        abort(400, "Received an invalid module: " + str(e))


@routes.route("/string", methods=["POST"])
def bayespairing_string():
    '''
    Represents the primary endpoint for calling BayesPairing with string input.

    Inputs are provided in the form of query strings and are not explicitly defined to allow for flexibility.
    (TODO issue #5 Implement support for uploading a user specified dataset)

    :returns: a jsonified dictionary of modules and their hits.
    '''
    if not (request.form):
        abort(400, description="Did not receive any arguments.")
    return bayespairing(request.form)


@routes.route("/file", methods=["POST"])
def bayespairing_file():
    '''
    Represents the primary endpoint for calling BayesPairing with file input.

    Inputs are provided in the request body. The file is provided as an individual argument with key input, 
    with arguments provided as a JSON dictionary with key "arguments". Arguments are the same as in the string input case.

    :returns: a jsonified dictionary of modules and their hits.
    '''
    if "input" not in request.files:
        abort(400, description="Did not receive an input file.")

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
            print("BayesPairing File: Stored temporary file: %s\n" % (temp.name))

        # for consistency, we convert from a string to an ImmutableMultiDict (and a dictionary in-between)
        # so we can reuse the BaysesPairing method
        return bayespairing(request.form, file.filename.rsplit(".")[1], temp)

    abort(400, description="BayesPairing failed for an unknown reason. Please check your inputs.")


def bayespairing(input, input_file_type=None, input_file=None):
    '''
    Processes input and runs BayesPairing.

    :param input: a dictionary (ImmutableMultiDict) containing a mapping of arguments to their value.
    :param input_file_type: the type of file if one is provided (fasta or stockholm).
    :param input_file: a path to a stockholm or fasta file if no sequence string is provided.
    :returns: a jsonified dictionary of modules and their hits.
    '''
    try:
        arguments = {}

        if (request.files.get("dataset")):
            flash(
                "User provided datasets are not yet supported. The default will be used instead.")

        # this primarily exits for use in pipelining data to VeRNAl & RNAMigos to avoid additional calls to the server
        get_graphs = input.get("get_graphs", default=0, type=int)

        sequence = input.get("sequence", type=str)
        if (not sequence and not input_file):
            abort(400, description="Did not receive a sequence as an argument or a file.")
        elif (sequence and input_file):
            abort(
                400, description="Received both a sequence and a file input. Are you using the right endpoint?")
        elif (not sequence):
            sequence = input_file.name

        secondary_structure = input.get(
            "secondary_structure", default="", type=str)
        secondary_structure_infile = input.get(
            "secondary_structure_infile", default=0, type=int)

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
        dataset = input.get("dataset", default="reliable_dataset", type=str)

        # load all arguments if they were provided or set the default value
        arguments["t"] = input.get("score_threshold", default=-2.3, type=float)
        arguments["samplesize"] = input.get(
            "sample_size", default=20000, type=int)
        arguments["theta"] = input.get("theta", default=1, type=int)
        arguments["lambda"] = input.get("lambda", default=0.35, type=float)
        arguments["w"] = input.get("window_length", default=200, type=int)
        arguments["s"] = input.get("step_size", default=100, type=int)
        if (input.get("modules")):
            arguments["mod"] = eval(input.get("modules"))
        else:
            arguments["mod"] = None
        arguments["constraints"] = input.get(
            "constraints", default="", type=str)
        arguments["aln"] = input.get("alignment_provided", default=0, type=int)
        arguments["verbose"] = BAYESPAIRING_VERBOSE

        if (arguments["aln"] and get_graphs):
            flash("Received an alignment and a request to produce representative graphs. Graphs cannot be produced, but BayesPairing will be run.")

        if (arguments["mod"]):
            modules_to_check = [int(input_number)
                                for input_number in arguments["mod"]]
        else:
            # we load the modules from the dataset to get the number of modules available.
            try:
                with open(os.path.join(CURRENT_DIRECTORY, "../core/models/" + dataset + ".json")) as f:
                    modules = json.load(f)
            except Exception as e:
                abort(
                    400, description="Could not process dataset due to error: " + str(e))
            modules_to_check = range(len(modules))

        # run BP2
        sequences, all_results = bp.run_fasta(
            sequence, modules_to_check, dataset, secondary_structure, arguments, input_file_type)

        # if an alignment is not provided, we can generate an SVG
        # in the future, we may generate an SVG with an alignment but this is not possible right now
        if not arguments["aln"]:
            all_svg_hits = {}
            all_chef_ss = []

            # generate SVG
            for seqCounter, inputSeqKey in enumerate(list(all_results.keys())):
                (modules_in_svg, chef_ss), temp = chefschoice.bp_chefs_choice(
                    all_results[inputSeqKey], sequences[seqCounter], arguments["t"])

                # now we need to fill svg hits
                svg_hits = {}
                for hit in modules_in_svg:
                    modID, modInfo = hit
                    if modID not in svg_hits:
                        svg_hits[modID] = [modInfo]
                    else:
                        svg_hits[modID].append(modInfo)
                all_svg_hits[inputSeqKey] = svg_hits
                all_chef_ss.append(chef_ss)

                # we only return an SVG for the first sequence for now (this is purely because this is the expected behaviour in the pipeline)
                # once ready to receive all SVGs, this is a trivial fix
                if (seqCounter == 0):
                    with open(temp.name, "rb") as image_file:
                        b64_encoded_svg = base64.b64encode(
                            image_file.read()).decode("utf-8")
                # destroy the temporarily stored svg
                temp.close()

            # TODO: SVG as XML
            output_dict = {"bp_input": sequences, "bp_params": arguments, "chefs_choice_struct": all_chef_ss,
                           "bp_all_hits": all_results, "bp_svg_hits": all_svg_hits, "bp_svg_b64":  b64_encoded_svg}
        else:  # if the input is an alignment, then no SVG
            output_dict = {"bp_input": sequences, "bp_params": arguments, "chefs_choice_struct": [
            ], "bp_all_hits": all_results,  "bp_svg_hits": {}, "bp_svg_b64":  ""}

        if (input_file):
            input_file.close()
        if not output_dict:
            abort(400, description="BayesPairing failed to produce output. Please ensure sure all input variables are valid. If a fasta or stockholm file were provided, they were preserved.")

        motif_graphs = {}
        # get representative graphs for motifs from chefschoice, this is as a result not possible for alignments yet
        if (get_graphs and not arguments["aln"]):
            for sequence in all_svg_hits.keys():
                module_ids = list(all_svg_hits[sequence].keys())
                print(module_ids)
                graphs = retrieve_graphs(module_ids)
                motif_graphs[sequence] = graphs

        output_dict["motif_graphs"] = motif_graphs
        return jsonify(output_dict)

    except ValueError as e:
        if (input_file):
            input_file.close()
        abort(400, description="Received an invalid argument: " + str(e))
    except Exception as e:
        if (input_file):
            input_file.close()
        abort(400, description="BayesPairing failed to produce a result: " + str(e))


def retrieve_graphs(modules):
    '''
    Retrieves representative graphs given a list of modules.

    :param modules: a list of modules (can also be a list containing a single module)
    :returns: a dictionary containing a mapping of modules to their respective graphs
    '''
    try:
        dataset_path = os.path.join(
            CURRENT_DIRECTORY, "../core/models/reliable_dataset.json")
        dataset = json.load(open(dataset_path, "r"))
        module_graph_mapping = {}

        for module in modules:
            graph = dataset[str(module)]["master_graph"]
            graph_dict = {}
            graph_dict["nodes"] = list(graph["nodes"])
            graph_dict["edges"] = list(graph["edges"])
            module_graph_mapping[module] = graph_dict

        dataset.close()
        return module_graph_mapping
    except ValueError as e:
        abort(
            400, description="Failed to process module IDs due to a type error: " + str(e))
    except Exception as e:
        abort(400, description="Failed to generate representative graphs for modules: " + str(e))


def allowed_file(filename):
    '''
    Verifies that a file has the correct format.

    :param filename: the name of the file to check.
    :returns: true if the file format is valid.
    '''
    return "." in filename and \
           filename.rsplit(".")[1].lower() in ALLOWED_EXTENSIONS
