import os
import json
import base64
from . import pipeline_bp as bp
from . import pipeline_chefschoice as chefs_choice

CURRENT_DIRECTORY = os.path.dirname(__file__)
ALLOWED_EXTENSIONS = {"fa", "fasta", "sl"}
ALLOWED_DATASETS = {"RELIABLE", "ALL"}
DEFAULT_DATASET = "ALL"
VERBOSE = False


def bayespairing(input, input_file_type=None, input_file=None):
    '''
    Processes input and runs BayesPairing.

    :param input: a dictionary (ImmutableMultiDict) containing a mapping of arguments to their value.
    :param input_file_type: the type of file if one is provided (fasta or stockholm).
    :param input_file: a path to a stockholm or fasta file if no sequence string is provided.
    :returns: a dictionary containing BayesPairing output
    '''
    try:
        # this primarily exits for use in pipelining data to VeRNAl & RNAMigos to avoid additional calls to the server
        get_graphs = input.get("get_graphs", default=0, type=int)
        sequence = input.get("sequence", type=str)

        if (not sequence and not input_file):
            raise Exception(
                "Did not receive a sequence as an argument or a file.")
        elif (sequence and input_file):
            raise Exception(
                "Received both a sequence and a file input. Are you using the right endpoint?")
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
            raise Exception(
                "Indicated that the secondary structure is provided in the fasta file, but also provided a separate structure.")
        elif ((input_file_type == "fa" or input_file_type == "fasta") and secondary_structure != ""):
            input_file.close()
            raise Exception(
                "Received a secondary structure as a string, but a fasta file as a sequence. Please provide the single sequence as a string.")
        elif (secondary_structure_infile):
            secondary_structure = "infile"

        # a user can select from a locally stored dataset
        dataset = input.get("dataset", default=DEFAULT_DATASET, type=str)
        if not (dataset in ALLOWED_DATASETS):
            raise Exception(
                "BayesPairing received an invalid dataset as an argument.")

        score_threshold = input.get(
            "score_threshold", default=-2.3, type=float)
        is_alignment = input_file_type == "sl"

        if (input.get("modules")):
            modules_to_check = [int(input_number)
                                for input_number in eval(input.get("modules"))]
        else:
            # we load the modules from the dataset to get the number of modules available.
            try:
                with open(os.path.join(CURRENT_DIRECTORY, f"../../models/{dataset}.json")) as f:
                    modules = json.load(f)
            except Exception as e:
                raise Exception(
                    "Could not process dataset due to error: " + str(e))
            modules_to_check = range(len(modules))

        # run BayesPairing
        sequences, all_results = bp.run_fasta(
            sequence, input, modules_to_check, dataset, is_alignment, score_threshold, secondary_structure, input_file_type, VERBOSE)

        # if an alignment is not provided, we can generate an SVG
        # in the future, we may generate an SVG with an alignment but this is not possible right now
        if not is_alignment:
            all_svg_hits = {}
            all_chef_ss = []

            # generate SVG
            for seqCounter, inputSeqKey in enumerate(list(all_results.keys())):
                (modules_in_svg, chef_ss), temp = chefs_choice.bp_chefs_choice(
                    all_results[inputSeqKey], sequences[seqCounter], score_threshold)

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
            output_dict = {"input": sequences, "params": input, "chefs_choice_struct": all_chef_ss,
                           "all_hits": all_results, "svg_hits": all_svg_hits, "svg_b64":  b64_encoded_svg}
        else:  # if the input is an alignment, then no SVG
            output_dict = {"input": sequences, "params": input, "chefs_choice_struct": [
            ], "all_hits": all_results,  "svg_hits": {}, "svg_b64":  ""}

        if (input_file):
            input_file.close()
        if not output_dict:
            raise Exception(
                "BayesPairing failed to produce output. Please ensure sure all input variables are valid. If a fasta or stockholm file were provided, they were preserved.")

        motif_graphs = {}
        # get representative graphs for motifs from chefschoice, this is as a result not possible for alignments yet
        if (get_graphs and not is_alignment):
            for sequence in all_svg_hits.keys():
                module_ids = list(all_svg_hits[sequence].keys())
                print(module_ids)
                graphs = retrieve_graphs(dataset, module_ids)
                motif_graphs[sequence] = graphs

        output_dict["motif_graphs"] = motif_graphs
        return output_dict
    except ValueError as e:
        if (input_file):
            input_file.close()
        raise Exception("Received an invalid argument: " + str(e))
    except Exception as e:
        if (input_file):
            input_file.close()
        raise Exception("BayesPairing failed to produce a result: " + str(e))


def retrieve_graphs(dataset, modules, include_pdb):
    '''
    Retrieves representative graphs given a list of modules.

    :param dataset: the name of the dataset to use
    :param modules: a list of modules (can also be a list containing a single module)
    :param include_pdb: an int indicating whether pdbs should be included in the graph output
    :returns: a dictionary containing a mapping of modules to their respective graphs
    '''
    if not (dataset in ALLOWED_DATASETS):
        raise Exception(
            "BayesPairing received an invalid dataset as an argument.")
    dataset_path = os.path.join(
        CURRENT_DIRECTORY, f"../../models/{dataset}.json")
    dataset = json.load(open(dataset_path, "r"))
    module_graph_mapping = {}

    for module in modules:
        graph = dataset[str(module)]["master_graph"]
        graph_dict = {}
        graph_dict["nodes"] = graph["nodes"]
        graph_dict["edges"] = graph["edges"]
        if (include_pdb):
            graph_dict["PDBs"] = dataset[str(module)]["PDBs"]
        module_graph_mapping[module] = graph_dict

    return module_graph_mapping


def allowed_file(filename):
    '''
    Verifies that a file has the correct format.

    :param filename: the name of the file to check.
    :returns: true if the file format is valid.
    '''
    return "." in filename and \
           filename.rsplit(".")[1].lower() in ALLOWED_EXTENSIONS
