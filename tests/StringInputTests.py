import unittest
import json
import os
from .constants.TestPayloads import *
from .constants.TestEnvironment import *
# relative paths refuse to work because Python is terrible!
from app import app

CURRENT_DIRECTORY = os.path.dirname(__file__)

# IMPORTANT NOTE: BP Results are not deterministic. To that effect, these tests validate the existence of keys and whether responses are successful or a failure
# NOTE: For boolean based comparisons, we still use assertEquals so there is no need to update code in the case of a new response (i.e. you only need to change the expected response JSONs)


class StringInputTests(unittest.TestCase):
    '''
    Implements unit testing for the string input endpoint.
    '''

    def setUp(self):
        self.app = app.test_client()

    def test_successful_request_no_ss(self):
        '''
        Tests for a successful run of BayesPairing where a sequence, and no secondary structure are provided.
        '''
        payload = STRING_INPUT_NO_SS
        headers = {}

        response = self.app.post(
            STRING_URL, content_type='multipart/form-data', headers=headers, data=payload)
        with open(os.path.join(CURRENT_DIRECTORY, "responses/RESPONSE_STRING_NO_SS.json")) as f:
            expected_response = json.load(f)
        f.close()

        self.assertEqual(200, response.status_code)
        self.assertEqual("motif_graphs" in expected_response,
                         "motif_graphs" in response.json)
        self.compare_core(expected_response, response.json)

    def test_successful_request_ss(self):
        '''
        Tests for a successful run of BayesPairing where a sequence, and a secondary structure are provided.
        '''
        payload = STRING_INPUT_WITH_SS
        headers = {}

        response = self.app.post(
            STRING_URL, content_type='multipart/form-data', headers=headers, data=payload)
        with open(os.path.join(CURRENT_DIRECTORY, "responses/RESPONSE_STRING_WITH_SS.json")) as f:
            expected_response = json.load(f)
        f.close()

        self.assertEqual(200, response.status_code)
        # BP results are NOT determnisitc, the best we can do is validate success (except for params & input)
        self.assertEqual("motif_graphs" in expected_response,
                         "motif_graphs" in response.json)
        self.compare_core(expected_response, response.json)

    def test_successful_request_custom_modules(self):
        '''
        Tests for a successful run of BayesPairing where a sequence, and a specified list of modules are provided.
        '''
        payload = STRING_INPUT_CUSTOM_MODULES
        headers = {}

        response = self.app.post(
            STRING_URL, content_type='multipart/form-data', headers=headers, data=payload)
        with open(os.path.join(CURRENT_DIRECTORY, "responses/RESPONSE_STRING_CUSTOM_MODULES.json")) as f:
            expected_response = json.load(f)
        f.close()

        self.assertEqual(200, response.status_code)
        # BP results are NOT determnisitc, the best we can do is validate success (except for params & input)
        self.assertEqual("motif_graphs" in expected_response,
                         "motif_graphs" in response.json)
        self.compare_core(expected_response, response.json)

    def test_successful_request_sequence_only(self):
        '''
        Tests for a successful run of BayesPairing where only a sequence, and no other optional params
        '''
        payload = STRING_INPUT_SEQUENCE_ONLY
        headers = {}

        response = self.app.post(
            STRING_URL, content_type='multipart/form-data', headers=headers, data=payload)
        with open(os.path.join(CURRENT_DIRECTORY, "responses/RESPONSE_STRING_SEQUENCE_ONLY.json")) as f:
            expected_response = json.load(f)
        f.close()

        self.assertEqual(200, response.status_code)
        # BP results are NOT determnisitc, the best we can do is validate success (except for params & input)
        self.assertEqual("motif_graphs" not in expected_response,
                         "motif_graphs" not in response.json)
        self.compare_core(expected_response, response.json)

    def test_failure_empty_modules(self):
        '''
        Tests for a failure run of BayesPairing where an empty list of modules is provided.
        '''
        payload = STRING_INPUT_EMPTY_MODULES
        headers = {}

        response = self.app.post(
            STRING_URL, content_type='multipart/form-data', headers=headers, data=payload)

        self.assertEqual(400, response.status_code)

    def test_failure_incorrect_dataset(self):
        '''
        Tests for a failure run of BayesPairing where an incorrect dataset is provided.
        '''
        payload = STRING_INPUT_INCORRECT_DATASET
        headers = {}

        response = self.app.post(
            STRING_URL, content_type='multipart/form-data', headers=headers, data=payload)

        self.assertEqual(400, response.status_code)

    def test_failure_no_arguments(self):
        '''
        Tests for a failure run of BayesPairing where no arguments are provided.
        '''
        payload = {}
        headers = {}

        response = self.app.post(
            STRING_URL, content_type='multipart/form-data', headers=headers, data=payload)

        self.assertEqual(400, response.status_code)

    def compare_core(self, expected, received):
        self.assertEqual("all_hits" in expected, "all_hits" in received)
        self.assertEqual("chefs_choice_struct" in expected,
                         "chefs_choice_struct" in received)
        self.assertEqual("svg" in expected, "svg" in received)
        self.assertEqual("svg_hits" in expected,
                         "svg_hits" in received)
        self.assertEqual(expected["input"], received["input"])
        self.assertEqual(expected["params"], received["params"])
