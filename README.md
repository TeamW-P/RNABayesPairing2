# `BayesPairing2`: Identifying RNA 3D modules in sequence and alignments 

This package includes tools for:

* Identifying RNA 3D modules in a sequence
* Creating your own BayesPairing models
* RNA 3D module plotting 
* Searching for modules in many sequences separately (fasta format)
* Searching for modules in alignments (stockholm format)
* Identifying modules in sequences with known secondary structure.

NOTE: The user friendly API is currently being built and tested. Sample usage described below is still being improved.

### BayesPairing2 was accepted to RECOMB 2020

The article is available on biorxiv:
https://www.biorxiv.org/content/10.1101/834762v1.abstract


## Requirements

* Python 3.6+ 

* Python modules (installed by our installer)
    * Networkx 2.1+
    * BioPython
    * matplotlib
    * weblogo
    * wrapt
    * anytree
    * pgmpy (included)
    
* RNAlib (install with conda : https://anaconda.org/bioconda/viennarna)
* LaTeX (optional)

## Installing

```
git clone git@jwgitlab.cs.mcgill.ca:sarrazin/rnabayespairing2.git
cd rnabayespairing2
pip install .
```

``pip install .`` in this directory will install the required python libraries, but you will need to install RNAlib separately, and have it in your PATH to use BayesPairing2. You can install RNAlib via conda: https://anaconda.org/bioconda/viennarna)

## Provided datasets

BayesPairing2 comes with three pre-assembled datasets you can immediately start searching sequences with. To search with a specific dataset, use the ``-d`` option with the name of the dataset.

* ``3DmotifAtlas_RELIABLE``: A subset of 60 modules from the RNA 3D Motif Atlas with the highest number of occurrences and highest sequence variation. We are confident in the prediction of those modules given the high quality data we have to train them. This is the default dataset we use.
* ``3DmotifAtlas_ALL``: A dataset containing all the modules we were able to convert from the 3D Motif Atlas to BayesPairing2 models (426). Some of those only had one occurrence and/or may have been trained on limited/incomplete data.
* ``rna3dmotif`` : A dataset containing the 75 most recurrent modules as identified via an exhaustive search of loops in the full PDB database with rna3dmotif.

#### Interpreting dataset-specific output

BayesPairing2 returns results by index, where the indexes correspond to modules in the relevant database. The rna3dmotif modules are described by graphs and sequence logos found in ``bayespairing/DBData/rna3dmotif``. 

The 3D Motif Atlas modules match to entries in that database. The correspondences between indexes of the two BayesPairing2 Atlas databases and the online 3D motif atlas database (with link to each relevant model) are found in the file ``bayespairing/DBData/3DmotifAtlas/3DmotifAtlas_info.csv``.

In this csv file, you can observe that sometimes, more than one module maps to the same entry of the Atlas; this is because the 3D Motif Atlas modules are clustered in 3D, and sometimes it is not possible to represent all occurrences accurately with the same graph, so they must be searched separately.


## Identifying 3D modules in a sequence

The core function of the BayesPairing package is ``parse_sequences.py``. It takes as input a sequence and outputs the position of putative motifs, as well as a probabilistic score evaluating the relative likelihood of finding them in the 3D structure associated with this sequence. It should be executed in the ``bayespairing/src`` folder

``parse_sequences.py [-h] [--verbose] [-seq SEQ] [-ss SS] [-d D] [-mod MOD [MOD ...]] [--aln] [-t T] [-w W] [-s S]
                          [-lambda LAMBDA] [-o O] [--interm] [--pretrained] [-delta DELTA] [-samplesize SAMPLESIZE]
                          [-theta THETA] [-constraints CONSTRAINTS] [--init]``

  * ``-h, --help``           show this help message and exit
  * ``--verbose``             increase output verbosity
  * ``-seq SEQ``               sequences to parse, string (one sequence), or FASTA file (many sequences), or STOCKHOLM (alignment)
  * ``-ss SS``              include secondary structure as a string, or write -ss infile if you include structures to the FASTA   
  * ``-d D ``               Name of the dataset. As a default, the dataset presented in the paper will be used.
  * ``-mod MOD [MOD ...]``    If you only want to parse for specific modules, list them as their dataset index. ex -mod 1 4 12
  * ``--aln``                 if an alignment file (stockholm) is given to -seq
  * ``-t T``                  Score threshold for a module to be called. [-10 to 10]. Default:-2.3 
  * ``-lambda LAMBDA``        weight of the secondary structure score(lambda).[0 to 1] Default:0.35
  * ``-samplesize SAMPLESIZE`` Size of the structure sample. Default value: 20000, although 1000 is ~10 times faster and usually works fine.
  * ``-w W``                   Window length when parsing long sequences [50 to 300]. Default:200
  * ``-s S``                  Step size between windows [10 to 200]. Default:100
  * ``-constraints CONSTRAINTS`` RNAfold hard constraints on input sequence
  * ``--init``          To reset all trained models on the modules included in -mod
  * ``-o O``                  Name of the output file. Default: output  

#### Example: searching for motifs on a TPP riboswitch sequence
 
The scripts described in this section should be run from the ``bayespairing/src`` directory. 
  
The first time you use BayesPairing with a full dataset, it will train all its models before searching a sequence. Those models will not need to be trained again. If you want to reset those models, you can use the ``init`` option. With the ``-d`` option, we are using the rna3dmotif dataset, which includes 75 pre-trained modules.

``python3 parse_sequences.py -seq "UUUUUUAAGGAAGAUCUGGCCUUCCCACAAGGGAAGGCCAAAGAAUUUCCUU" -samplesize 1000 -d rna3dmotif``

The output is very large, so we can raise the threshold to have a better idea of the dominating modules.

``python3 parse_sequences.py -seq "UUUUUUAAGGAAGAUCUGGCCUUCCCACAAGGGAAGGCCAAAGAAUUUCCUU" -samplesize 1000 -t 4 -d rna3dmotif``

```
=========================================================================================
|       MODULE   SCORE                MODULE POSITIONS                   SEQUENCE       |
|          9     4.196                     12-14,43-45                    AAU&GAU       |
|         13     4.721                     15-16,39-42                    CU&AAAG       |
|         22     4.256                     15-16,40-42                     AAG&CU       |
|         28     8.107                           25-30                     CACAAG       |
=========================================================================================


TOTAL TIME: 5.775
```

If we have a secondary structure, the search is much faster

``python3 parse_sequences.py -seq "UUUUUUAAGGAAGAUCUGGCCUUCCCACAAGGGAAGGCCAAAGAAUUUCCUU" -t 4 -d rna3dmotif -ss "......(((((((.((((((((((((....)))))))))..))).)))))))"``

```
=========================================================================================
|       MODULE   SCORE                MODULE POSITIONS                   SEQUENCE       |
|          9     4.254                     12-14,43-45                    AAU&GAU       |
|         28     8.108                           25-30                     CACAAG       |
=========================================================================================


TOTAL TIME: 2.581
```

### Displaying a dataset

To assess what module the module ID matches, we can generate graphs and sequence logos for all modules and store them in the Graphs directory.

``python3 display_modules.py -n "rna3dmotif"``

![](bayespairing/DBData/rna3dmotif/default_logo28.png)
![](bayespairing/DBData/rna3dmotif/default_graph28.png)

Module 28 is a hairpin with a trans sugar-hoogstein non-canonical base pair, as well as a stacking.

As for the RNA 3D Motif Atlas datasets, you can either use the website links in the file ``bayespairing/DBData/3DmotifAtlas/3DmotifAtlas_info.csv`` or generate the graphs and logos with ``display_modules.py``.

### Building your own dataset

For building new datasets, full dataset cross-validation, results presented in the paper and other advance uses, see the python notebook ``BayesPairing2_usage.ipynb`` in the ``bayespairing/test`` directory.


### Contact

Roman Sarrazin-Gendron
roman.sarrazingendron@mail.mcgill.ca
