FASTA_EXTENSIONS = ["afa", "fa", "fasta"]
STOCKHOLM_EXTENSIONS = ["stk", "sto"]
ALLOWED_EXTENSIONS = FASTA_EXTENSIONS + STOCKHOLM_EXTENSIONS
ALLOWED_DATASETS = ["RELIABLE", "ALL"]
# bp expects to return a list of params used. remove the ones not relevant (i.e. params used internally only)
PARAMS_TO_REMOVE = ["get_graphs", "sequence", "pdb", "vernal", "rnamigos"]
DEFAULT_DATASET = "ALL"