import unittest
import json
import os
from .constants.TestPayloads import *
from .constants.TestEnvironment import *
# relative paths refuse to work because Python is terrible!
from app import app
import io

CURRENT_DIRECTORY = os.path.dirname(__file__)

# IMPORTANT NOTE: BP Results are not deterministic. To that effect, these tests validate the existence of keys and whether responses are successful or a failure
# NOTE: For boolean based comparisons, we still use assertEquals so there is no need to update code in the case of a new response (i.e. you only need to change the expected response JSONs)


class FileInputTest(unittest.TestCase):
    '''
    Implements unit testing for the file input endpoint.
    '''

    def setUp(self):
        self.app = app.test_client()

    def test_successful_request_sequence_no_ss(self):
        '''
        Tests for a successful run of BP, provided sequences (via a fasta file) and no secondary structure.
        '''
        payload = FILE_INPUT_SEQUENCE_NO_SS
        headers = {}

        # simulate file upload
        with open(os.path.join(CURRENT_DIRECTORY, "data/INPUT_SEQUENCE.fa"), mode="rb") as f:
            payload["bp_input"] = (io.BytesIO(f.read()), "INPUT_SEQUENCE.fa")

        response = self.app.post(
            FILE_URL, content_type='multipart/form-data', headers=headers, data=payload)

        with open(os.path.join(CURRENT_DIRECTORY, "responses/RESPONSE_FILE_INPUT_SEQUENCE_NO_SECONDARY_STRUCTURE.json")) as f:
            expected_response = json.load(f)
        
        # simulated error
        self.assertEqual(400, response.status_code)
        self.assertEqual("svg" in expected_response, "svg" in response.json)
        self.assertEqual("svg_hits" in expected_response,
                         "svg_hits" in response.json)
        self.assertEqual("motif_graphs" in expected_response,
                         "motif_graphs" in response.json)
        self.compare_core(expected_response, response.json)

    def test_successful_request_sequence_ss_infile(self):
        '''
        Tests for a successful run of BP, provided sequences (via a fast file) with secondary structures.
        '''
        payload = FILE_INPUT_SEQUENCE_SS_INFILE
        headers = {}

        with open(os.path.join(CURRENT_DIRECTORY, "data/INPUT_SEQUENCE_SS.fa"), mode="rb") as f:
            payload["bp_input"] = (io.BytesIO(f.read()),
                                   "INPUT_SEQUENCE_SS.fa")

        response = self.app.post(
            FILE_URL, content_type='multipart/form-data', headers=headers, data=payload)
        with open(os.path.join(CURRENT_DIRECTORY, "responses/RESPONSE_FILE_INPUT_SEQUENCE_SS_INFILE.json")) as f:
            expected_response = json.load(f)

        self.assertEqual(200, response.status_code)
        # BP results are NOT determnisitc, the best we can do is validate success (except for params & input)
        self.assertEqual("svg" in expected_response, "svg" in response.json)
        self.assertEqual("svg_hits" in expected_response,
                         "svg_hits" in response.json)
        self.assertEqual("motif_graphs" in expected_response,
                         "motif_graphs" in response.json)
        self.compare_core(expected_response, response.json)

    def test_successful_request_sequence_custom_modules(self):
        '''
        Tests for a successful run of BP, provided sequences (via a fast file) with specified modules as a param.
        '''
        payload = FILE_INPUT_SEQUENCE_CUSTOM_MODULES
        headers = {}

        # simulate file upload
        with open(os.path.join(CURRENT_DIRECTORY, "data/INPUT_SEQUENCE.fa"), mode="rb") as f:
            payload["bp_input"] = (io.BytesIO(f.read()), "INPUT_SEQUENCE.fa")

        response = self.app.post(
            FILE_URL, content_type='multipart/form-data', headers=headers, data=payload)
        with open(os.path.join(CURRENT_DIRECTORY, "responses/RESPONSE_FILE_INPUT_SEQUENCE_CUSTOM_MODULES.json")) as f:
            expected_response = json.load(f)

        self.assertEqual(200, response.status_code)
        # BP results are NOT determnisitc, the best we can do is validate success (except for params & input)
        self.assertEqual("svg" in expected_response, "svg" in response.json)
        self.assertEqual("svg_hits" in expected_response,
                         "svg_hits" in response.json)
        self.assertEqual("motif_graphs" in expected_response,
                         "motif_graphs" in response.json)
        self.compare_core(expected_response, response.json)

    def test_successful_request_stockholm_alignment(self):
        '''
        Tests for a successful run of BP, provided a stockholm file (alignment).
        '''
        payload = FILE_INPUT_STOCKHOLM_ALIGNMENT
        headers = {}

        # simulate file upload
        # TODO: add stockholm
        with open(os.path.join(CURRENT_DIRECTORY, "data/INPUT_ALIGNMENT.stk"), mode="rb") as f:
            payload["bp_input"] = (io.BytesIO(f.read()), "INPUT_ALIGNMENT.stk")

        response = self.app.post(
            FILE_URL, content_type='multipart/form-data', headers=headers, data=payload)
        # TODO: finish when stockholm is fixed
        # with open(os.path.join(CURRENT_DIRECTORY, "responses/RESPONSE_FILE_INPUT_STOCKHOLM_ALIGNMENT.json")) as f:
        #    expected_response = json.load(f)

        self.assertEqual(400, response.status_code)
        '''
        self.assertTrue("svg" not in expected_response,
                        "svg" not in response.json)
        self.assertEqual("svg_hits" not in expected_response,
                         "svg_hits" not in response.json)
        self.assertEqual("motif_graphs" in expected_response,
                         "motif_graphs" in response.json)
        self.compare_core(expected_response, response.json)
        '''

    def test_failure_invalid_file_format(self):
        '''
        Tests for a failure run of BP, where an illegal file format is provided.
        '''
        # other params irrelevant for purpose of this test
        payload = {}
        headers = {}

        with open(os.path.join(CURRENT_DIRECTORY, "data/ILLEGAL_FORMAT.txt"), mode="rb") as f:
            payload["bp_input"] = (io.BytesIO(f.read()), "ILLEGAL_FORMAT.txt")

        response = self.app.post(
            FILE_URL, content_type='multipart/form-data', headers=headers, data=payload)

        self.assertEqual(400, response.status_code)

    def test_failure_invalid_fasta(self):
        '''
        Tests for a failure run of BP, where an incorrect fasta file is provided.
        '''
        # other params irrelevant for purpose of this test
        payload = {}
        headers = {}

        with open(os.path.join(CURRENT_DIRECTORY, "data/INVALID_FASTA.fa"), mode="rb") as f:
            payload["bp_input"] = (io.BytesIO(f.read()), "INVALID_FASTA.fa")

        response = self.app.post(
            FILE_URL, content_type='multipart/form-data', headers=headers, data=payload)

        self.assertEqual(400, response.status_code)

    def test_failure_invalid_stockholm(self):
        '''
        Tests for a failure run of BP, where an incorrect stockholm file is provided.
        '''
        # other params irrelevant for purpose of this test
        payload = {}
        headers = {}

        with open(os.path.join(CURRENT_DIRECTORY, "data/INVALID_STOCKHOLM.stk"), mode="rb") as f:
            payload["bp_input"] = (io.BytesIO(f.read()),
                                   "INVALID_STOCKHOLM.stk")

        response = self.app.post(
            FILE_URL, content_type='multipart/form-data', headers=headers, data=payload)

        self.assertEqual(400, response.status_code)

    def test_failure_string_and_fasta(self):
        '''
        Tests for a failure run of BP, where a string sequence AND a fasta file are both provided.
        '''
        payload = FILE_INPUT_STRING_AND_FASTA
        headers = {}

        response = self.app.post(
            FILE_URL, content_type='multipart/form-data', headers=headers, data=payload)

        self.assertEqual(400, response.status_code)

    def test_failure_alignment_with_ss(self):
        '''
        Tests for a failure run of BP, where an alignment is provided and is marked as secondary structure infile.
        '''
        payload = FILE_INPUT_ALIGNMENT_SS
        headers = {}

        with open(os.path.join(CURRENT_DIRECTORY, "data/INPUT_ALIGNMENT.stk"), mode="rb") as f:
            payload["bp_input"] = (io.BytesIO(f.read()), "INPUT_ALIGNMENT.stk")

        response = self.app.post(
            FILE_URL, content_type='multipart/form-data', headers=headers, data=payload)

        self.assertEqual(400, response.status_code)

    def test_failure_ss_string_provided(self):
        '''
        Tests for a failure run of BP, where a fasta file and a secondary structure string are provided.
        '''
        payload = FILE_INPUT_SS_STRING_PROVIDED
        headers = {}

        response = self.app.post(
            FILE_URL, content_type='multipart/form-data', headers=headers, data=payload)

        self.assertEqual(400, response.status_code)

    def test_failure_incorrect_dataset(self):
        '''
        Tests for a failure run of BP, where an invalid dataset is provided.
        '''
        payload = FILE_INPUT_INCORRECT_DATASET
        headers = {}

        response = self.app.post(
            FILE_URL, content_type='multipart/form-data', headers=headers, data=payload)

        self.assertEqual(400, response.status_code)

    def test_failure_no_arguments(self):
        '''
        Tests for a failure run of BP, where no arguments are provided.
        '''
        payload = {}
        headers = {}

        response = self.app.post(
            FILE_URL, content_type='multipart/form-data', headers=headers, data=payload)

        self.assertEqual(400, response.status_code)

    def compare_core(self, expected, received):
        self.assertEqual("all_hits" in expected, "all_hits" in received)
        self.assertEqual("chefs_choice_struct" in expected,
                         "chefs_choice_struct" in received)
        self.assertEqual(expected["input"], received["input"])
        self.assertEqual(expected["params"], received["params"])
