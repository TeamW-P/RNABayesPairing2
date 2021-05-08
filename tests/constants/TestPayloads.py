# Stores all payloads for use in testing

# STRING INPUT

STRING_INPUT_NO_SS = {'dataset': 'ALL',
                      'sequence': 'UUUUUUAAGGAAGAUCUGGCCUUCCCACAAGGGAAGGCCAAAGAAUUUCCUU',
                      'score_threshold': '0',
                      'sample_size': '1000',
                      'theta': '1',
                      'lambda': '0.35',
                      'window_length': '200',
                      'step_size': '100',
                      'modules': '',
                      'constraints': '',
                      'get_graphs': '1',
                      'pdb': '1'}
STRING_INPUT_WITH_SS = {'dataset': 'ALL',
                        'sequence': 'UUUUUUAAGGAAGAUCUGGCCUUCCCACAAGGGAAGGCCAAAGAAUUUCCUU',
                        'secondary_structure': '......(((((((.((((((((((((....)))))))))..))).)))))))',
                        'score_threshold': '0',
                        'sample_size': '1000',
                        'theta': '1',
                        'lambda': '0.35',
                        'window_length': '200',
                        'step_size': '100',
                        'modules': '',
                        'constraints': '',
                        'get_graphs': '1',
                        'pdb': '1'}
STRING_INPUT_CUSTOM_MODULES = {'dataset': 'ALL',
                               'sequence': 'UUUUUUAAGGAAGAUCUGGCCUUCCCACAAGGGAAGGCCAAAGAAUUUCCUU',
                               'score_threshold': '0',
                               'sample_size': '1000',
                               'theta': '1',
                               'lambda': '0.35',
                               'window_length': '200',
                               'step_size': '100',
                               'modules': '["36", "48"]',
                               'constraints': '',
                               'get_graphs': '1',
                               'pdb': '1'}
STRING_INPUT_SEQUENCE_ONLY = {
    'sequence': 'UUUUUUAAGGAAGAUCUGGCCUUCCCACAAGGGAAGGCCAAAGAAUUUCCUU'}
STRING_INPUT_INCORRECT_DATASET = {'dataset': 'BAD',
                                  'sequence': 'UUUUUUAAGGAAGAUCUGGCCUUCCCACAAGGGAAGGCCAAAGAAUUUCCUU'}
STRING_INPUT_EMPTY_MODULES = {'dataset': 'ALL',
                              'sequence': 'UUUUUUAAGGAAGAUCUGGCCUUCCCACAAGGGAAGGCCAAAGAAUUUCCUU',
                              'score_threshold': '0',
                              'sample_size': '1000',
                              'theta': '1',
                              'lambda': '0.35',
                              'window_length': '200',
                              'step_size': '100',
                              'modules': '[]',
                              'constraints': '',
                              'get_graphs': '1',
                              'pdb': '1'}

# FILE INPUT

FILE_INPUT_SEQUENCE_NO_SS = {'dataset': 'ALL',
                             'secondary_structure_infile': '0',
                             'score_threshold': '0',
                             'sample_size': '1000',
                             'theta': '1',
                             'lambda': '0.35',
                             'window_length': '200',
                             'step_size': '100',
                             'modules': '',
                             'constraints': '',
                             'get_graphs': '1',
                             'pdb': '1'}
FILE_INPUT_SEQUENCE_SS_INFILE = {'dataset': 'ALL',
                                 'secondary_structure_infile': '1',
                                 'score_threshold': '0',
                                 'sample_size': '1000',
                                 'theta': '1',
                                 'lambda': '0.35',
                                 'window_length': '200',
                                 'step_size': '100',
                                 'modules': '',
                                 'constraints': '',
                                 'get_graphs': '1',
                                 'pdb': '1'}
FILE_INPUT_SEQUENCE_CUSTOM_MODULES = {'dataset': 'ALL',
                                      'secondary_structure_infile': '0',
                                      'score_threshold': '0',
                                      'sample_size': '1000',
                                      'theta': '1',
                                      'lambda': '0.35',
                                      'window_length': '200',
                                      'step_size': '100',
                                      'modules': '["36", "48"]',
                                      'constraints': '',
                                      'get_graphs': '1',
                                      'pdb': '1'}
FILE_INPUT_STOCKHOLM_ALIGNMENT = {'dataset': 'ALL',
                                  'secondary_structure_infile': '0',
                                  'score_threshold': '0',
                                  'sample_size': '1000',
                                  'theta': '1',
                                  'lambda': '0.35',
                                  'window_length': '200',
                                  'step_size': '100',
                                  'modules': '',
                                  'constraints': '',
                                  'aln': '1',
                                  'get_graphs': '1',
                                  'pdb': '1'}
FILE_INPUT_ALIGNMENT_SS = {'dataset': 'ALL',
                           'secondary_structure_infile': '1',
                           'score_threshold': '0',
                           'sample_size': '1000',
                           'theta': '1',
                           'lambda': '0.35',
                           'window_length': '200',
                           'step_size': '100',
                           'modules': '',
                           'constraints': '',
                           'get_graphs': '1',
                           'pdb': '1'}
FILE_INPUT_STRING_AND_FASTA = {'sequence': 'UUUUUUAAGGAAGAUCUGGCCUUCCCACAAGGGAAGGCCAAAGAAUUUCCUU',
                               'dataset': 'ALL',
                               'secondary_structure_infile': '0',
                               'score_threshold': '0',
                               'sample_size': '1000',
                               'theta': '1',
                               'lambda': '0.35',
                               'window_length': '200',
                               'step_size': '100',
                               'modules': '',
                               'constraints': '',
                               'get_graphs': '1',
                               'pdb': '1'}
FILE_INPUT_SS_STRING_PROVIDED = {'secondary_structure': '......(((((((.((((((((((((....)))))))))..))).)))))))',
                                 'dataset': 'ALL',
                                 'secondary_structure_infile': '0',
                                 'score_threshold': '0',
                                 'sample_size': '1000',
                                 'theta': '1',
                                 'lambda': '0.35',
                                 'window_length': '200',
                                 'step_size': '100',
                                 'modules': '',
                                 'constraints': '',
                                 'get_graphs': '1',
                                 'pdb': '1'}
FILE_INPUT_INCORRECT_DATASET = {'dataset': 'SPECIAL',
                                'secondary_structure_infile': '0',
                                'score_threshold': '0',
                                'sample_size': '1000',
                                'theta': '1',
                                'lambda': '0.35',
                                'window_length': '200',
                                'step_size': '100',
                                'modules': '',
                                'constraints': '',
                                'get_graphs': '1',
                                'pdb': '1'}

# GRAPH RETRIEVAL

GRAPH_RETRIEVAL_VALID_MODULES = "[\"36\", \"48\"]"
GRAPH_RETRIEVAL_INVALID_MODULE = "[\"36\", \"1000\"]"
GRAPH_RETRIEVAL_EMPTY_MODULE = []
GRAPH_RETRIEVAL_INCLUDE_PDB = 1
GRAPH_RETRIEVAL_NO_PDB = 0

# MODULE INFO RETRIEVAL

MODULE_INFO_RETRIEVAL_VALID_MODULES = "[\"36\", \"48\"]"
MODULE_INFO_RETRIEVAL_INVALID_MODULES = "[\"36\", \"1000\"]"
MODULE_INFO_RETRIEVAL_EMPTY_MODULE = "[]"

# MISC

ALL_DATASET_NAME = "ALL"
ILLEGAL_DATASET_NAME = "ILLEGAL"
