# `BayesPairing2`: Identifying RNA 3D modules in sequence and alignments 
[![Build Status](https://travis-ci.com/TeamW-P/RNABayesPairing2.svg?branch=main)](https://travis-ci.com/TeamW-P/RNABayesPairing2)

## This is a fork of BayesPairing 2 to be used within the context of a larger Bioionformatics pipeline. For info on BayesPairing, please click [here](https://jwgitlab.cs.mcgill.ca/sarrazin/rnabayespairing2/).

## Major Changes
Changes to the official branch of BayesPairing (continually updated)

* BayesPairing previously operated by use of absolute imports/references to file paths. These were switched to relative imports in-order to allow it to fit within the project architecture.

## What's New

* A flask layer that allows BayesPairing to be called from an endpoint. This includes endpoints for file input & string input as well as another to retrieve representative graphs from a module database.
* BayesPairing has also been dockerized, and a dockerfile is provided.
* Three additional classes are created within core/src/pipeline.
  * pipeline_bp: A modification of run_fasta from parse_sequences that does not save output to files
  * pipeline_chefschoice: A modification of key methods from chefs_choice that handle SVG generation
  * chefs_assistant: a new class that handles virtually all heavy lifting and processing of the BP service
* Testing! Located within the tests directory! Because BP results are non-deterministic, these mainly validate that results are successful or a failure and whether keys exist, not the exact contents of a specific result.
* TravisCI will automatically run tests after each commit

The command line program remains functional, but must now be run from the root directory instead of src. Here is a sample run: 

`python3 -m core.src.parse_sequences -seq "UUUUUUAAGGAAGAUCUGGCCUUCCCACAAGGGAAGGCCAAAGAAUUUCCUU" -samplesize 1000 -t 4 -d ALL`

Note that you must first set-up and activate the Anaconda environment (requirements provided in environment.yml).

This package includes tools for:

* Identifying RNA 3D modules in a sequence
* Creating your own BayesPairing models
* RNA 3D module plotting 
* Searching for modules in many sequences separately (fasta format)
* Searching for modules in alignments (stockholm format)
* Identifying modules in sequences with known secondary structure.

### BayesPairing2 was accepted to RECOMB 2020

The article is available on biorxiv:
https://www.biorxiv.org/content/10.1101/834762v1.abstract


## Requirements (installed via conda)

* Python 3.6+ 
* Docker
* Flask
* Gunicorn
* Networkx 2.1+
* BioPython
* matplotlib
* weblogo
* wrapt
* anytree
* pgmpy (included)  
* RNAlib (viennaRNA)
* LaTeX (optional)

### Contact

For inquiries on BayesPairing specifically, please contact:

Roman Sarrazin-Gendron
roman.sarrazingendron@mail.mcgill.ca

For inquiries regarding the overall pipeline as well as top level code including the Flask/Docker layer, please open an issue.
