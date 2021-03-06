import pickle
from flask import jsonify, abort, Blueprint, request, flash, request, current_app
from werkzeug.utils import secure_filename
from core.src.pipeline import chefs_assistant
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

@routes.errorhandler(400)
def resource_not_found(e):
    return jsonify(error=str(e)), 400

@routes.route("/graphs", methods=["POST"])
def get_graph_per_module():
    '''
    Given a list of modules, or a single module, returns graphs for each from the default dataset.

    This is a POST because of the possibility of uploading a dataset.

    :returns: a mapping of module IDs to their representative graphs.
    '''
    try:
        dataset = request.form.get("dataset", default = chefs_assistant.DEFAULT_DATASET, type = str)
        modules = request.form.get("modules")
        include_pdb = request.form.get("pdb", default = 1, type = int)
        if (not modules):
            abort(400, "Did not receive any modules to fetch graphs for.")
        modules = eval(modules)

        module_graph_mapping = chefs_assistant.retrieve_graphs(dataset, modules, include_pdb)

        return jsonify({"graphs": module_graph_mapping})
    except ValueError as e:
        raise ValueError(
            "Failed to process module IDs due to a type error: " + str(e))
    except Exception as e:
        abort(400, "Failed to generate representative graphs for modules: " + str(e))


@routes.route("/string", methods=["POST"])
def bayespairing_string():
    '''
    Represents the primary endpoint for calling BayesPairing with string input.

    Inputs are provided in the form of query strings and are not explicitly defined to allow for flexibility.

    :returns: a jsonified dictionary of modules and their hits.
    '''
    if not (request.form):
        abort(400, description="Did not receive any arguments.")
    try:
        result = chefs_assistant.bayespairing(request.form)
        return jsonify(result)
    except Exception as e:
        abort(400, description=str(e))


@routes.route("/file", methods=["POST"])
def bayespairing_file():
    '''
    Represents the primary endpoint for calling BayesPairing with file input.

    Inputs are provided in the request body. The file is provided as an individual argument with key input, 
    with arguments provided as a JSON dictionary with key "arguments". Arguments are the same as in the string input case.

    :returns: a jsonified dictionary of modules and their hits.
    '''
    if "bp_input" not in request.files:
        abort(400, description="Did not receive an input file.")

    file = request.files["bp_input"]

    if file.filename == "":
        abort(400, description="Received an invalid file.")
    if file and chefs_assistant.allowed_file(file.filename):
        # we store the file as a temporary file
        # this has inherent benefits, it handles uniqueness of naming and also destroys the file once we are done
        file.seek(0)
        temp = tempfile.NamedTemporaryFile()
        temp.write(file.read())
        temp.seek(0)

        try:
            result = chefs_assistant.bayespairing(
                request.form, file.filename.rsplit(".")[1], temp)
            return jsonify(result)
        except Exception as e:
            abort(400, description=str(e))

    abort(400, description="BayesPairing failed for an unknown reason. Please check your inputs.")
