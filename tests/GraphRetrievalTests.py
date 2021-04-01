import unittest
import json
import os
from .constants.TestPayloads import *
from .constants.TestEnvironment import *
# relative paths refuse to work because Python is terrible!
from app import app
import io

CURRENT_DIRECTORY = os.path.dirname(__file__)

# NOTE: Results here can be considered determnistic but note that dataset changes mean that the responses also need to be changed.
# NOTE: For boolean based comparisons, we still use assertEquals so there is no need to update code in the case of a new response (i.e. you only need to change the expected response JSONs)


class GraphRetrievalTests(unittest.TestCase):
    '''
    Implements unit testing for the graph retrieval endpoint.
    '''

    def setUp(self):
        self.app = app.test_client()

    def test_succesful_request_pdb(self):
        '''
        Tests for a successful run of graph retrieval with PDBs included in output.
        '''
        payload = GRAPH_RETRIEVAL_PDB
        headers = {}

        response = self.app.post(
            GRAPHS_URL, content_type='multipart/form-data', headers=headers, data=payload)

        with open(os.path.join(CURRENT_DIRECTORY, "responses/GRAPH_RETRIEVAL_SUCCESS_PDB.json")) as f:
            expected_response = json.load(f)
        f.close()

        self.assertEqual(200, response.status_code)
        self.assertEqual(expected_response, response.json)

    def test_succesful_request_no_pdb(self):
        '''
        Tests for a successful run of graph retrieval where no PDBs are included in output.
        '''
        payload = GRAPH_RETRIEVAL_NO_PDB
        headers = {}

        response = self.app.post(
            GRAPHS_URL, content_type='multipart/form-data', headers=headers, data=payload)

        with open(os.path.join(CURRENT_DIRECTORY, "responses/GRAPH_RETRIEVAL_SUCCESS_NO_PDB.json")) as f:
            expected_response = json.load(f)
        f.close()

        self.assertEqual(200, response.status_code)
        self.assertEqual(expected_response, response.json)

    def test_failure_empty_modules(self):
        '''
        Tests for a failure run of graph retrieval where an empty list of modules is provided.
        '''
        payload = GRAPH_RETRIEVAL_EMPTY_MODULES
        headers = {}

        response = self.app.post(
            GRAPHS_URL, content_type='multipart/form-data', headers=headers, data=payload)

        self.assertEqual(400, response.status_code)

    def test_failure_bad_dataset(self):
        '''
        Tests for a failure run of graph retrieval where an invalid dataset is provided.
        '''
        payload = GRAPH_RETRIEVAL_BAD_DATASET
        headers = {}

        response = self.app.post(
            GRAPHS_URL, content_type='multipart/form-data', headers=headers, data=payload)

        self.assertEqual(400, response.status_code)

    def test_failure_invalid_module(self):
        '''
        Tests for a failure run of graph retrieval where an invalid module is provided.
        '''
        payload = GRAPH_RETRIEVAL_INVALID_MODULE
        headers = {}

        response = self.app.post(
            GRAPHS_URL, content_type='multipart/form-data', headers=headers, data=payload)

        self.assertEqual(400, response.status_code)

    def test_failure_no_arguments(self):
        '''
        Tests for a failure run of graph retrieval where no arguments are provided.
        '''
        payload = {}
        headers = {}

        response = self.app.post(
            GRAPHS_URL, content_type='multipart/form-data', headers=headers, data=payload)

        self.assertEqual(400, response.status_code)
