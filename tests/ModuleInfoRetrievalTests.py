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


class ModuleInfoRetrieval(unittest.TestCase):
    '''
    Implements unit testing for module info retrieval.
    '''

    def setUp(self):
        self.app = app.test_client()

    def test_succesful_module_info_retrieval(self):
        '''
        Tests for a successful run of module info retrieval with PDBs included in output.
        '''
        payload = {}
        headers = {}
        url = f"{MODULE_INFO_RETRIEVAL_URL}?{MODULE_INFO_RETRIEVAL_DATASET_KEY}={ALL_DATASET_NAME}&{MODULE_INFO_RETRIEVAL_MODULES_KEY}={MODULE_INFO_RETRIEVAL_VALID_MODULES}"

        response = self.app.get(
            url, headers=headers, data=payload)

        with open(os.path.join(CURRENT_DIRECTORY, "responses/MODULE_INFO_RETRIEVAL_SUCCESS.json")) as f:
            expected_response = json.load(f)
        f.close()

        self.assertEqual(200, response.status_code)
        self.assertEqual(expected_response, response.json)

    def test_failure_empty_modules(self):
        '''
        Tests for a failure run of module info retrieval where an empty list of modules is provided.
        '''
        payload = {}
        headers = {}
        url = f"{MODULE_INFO_RETRIEVAL_URL}?{MODULE_INFO_RETRIEVAL_DATASET_KEY}={ALL_DATASET_NAME}&{MODULE_INFO_RETRIEVAL_MODULES_KEY}={MODULE_INFO_RETRIEVAL_EMPTY_MODULE}"

        response = self.app.get(
            url, headers=headers, data=payload)

        self.assertEqual(400, response.status_code)

    def test_failure_bad_dataset(self):
        '''
        Tests for a failure run of module info retrieval where an invalid dataset is provided.
        '''
        payload = {}
        headers = {}
        url = f"{MODULE_INFO_RETRIEVAL_URL}?{MODULE_INFO_RETRIEVAL_DATASET_KEY}={ILLEGAL_DATASET_NAME}&{MODULE_INFO_RETRIEVAL_MODULES_KEY}={MODULE_INFO_RETRIEVAL_VALID_MODULES}"

        response = self.app.get(
            url, headers=headers, data=payload)

        self.assertEqual(400, response.status_code)

    def test_failure_invalid_module(self):
        '''
        Tests for a failure run of module info retrieval where an invalid module is provided.
        '''
        payload = {}
        headers = {}
        url = f"{MODULE_INFO_RETRIEVAL_URL}?{MODULE_INFO_RETRIEVAL_DATASET_KEY}={ALL_DATASET_NAME}&{MODULE_INFO_RETRIEVAL_MODULES_KEY}={MODULE_INFO_RETRIEVAL_INVALID_MODULES}"

        response = self.app.get(
            url, headers=headers, data=payload)

        self.assertEqual(400, response.status_code)

    def test_failure_no_arguments(self):
        '''
        Tests for a failure run of module info retrieval where no arguments are provided.
        '''
        payload = {}
        headers = {}

        response = self.app.get(
            MODULE_INFO_RETRIEVAL_URL, headers=headers, data=payload)

        self.assertEqual(400, response.status_code)
