# `BayesPairing2`: Identifying RNA 3D modules in sequence and alignments 

## This is a fork of BayesPairing 2 to be used within the context of a larger Bioionformatics pipeline. For info on BayesPairing, please click [here](https://jwgitlab.cs.mcgill.ca/sarrazin/rnabayespairing2/).

## Major Changes
Changes to the official branch of BayesPairing (continually updated)

* BayesPairing previously operated by use of absolute imports/references to file paths. These were switched to relative imports in-order to allow it to fit within the project architecture.

## What's New

* A flask layer that allows BayesPairing to be called from an endpoint.
* BayesPairing has also been dockerized, and a dockerfile is provided.

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
