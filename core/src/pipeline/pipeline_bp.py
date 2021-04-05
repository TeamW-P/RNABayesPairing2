import os
import pickle
from multiprocessing import Process, Manager
from .. import testSS
from .. import BayesPairing
from ..parse_sequences import unpick, run_BP, process_siblings, parse_sequences, check_already_there, seq_ranges
from .constants import *
import random
from Bio import SeqIO
import argparse
import time as timer
import collections
from functools import reduce
from Bio import AlignIO
from operator import itemgetter

# this class is a modification of parse_sequences.py. It switches to direct outputs instead of saving files locally,
# and also accepts file input as an argument.
# helper methods are still called from parse_sequences


def run_fasta(input, arguments, modules_to_parse, dataset, aln, t, ss="", input_file_type=None, verbose=False):
    '''
    Runs BayesPairing given a set of arguments and a sequence or alignment.

    :param input: a string sequence or the path to a fasta or stockholm file
    :param arguments: an ImmutableMuliDict containing optional arguments to be used by BayesPairing
    :param modules_to_parse: a list of modules to parse
    :param dataset: the name of the dataset to use
    :param aln: whether an alignment (stockholm) is provided
    :param t: the score threshold
    :param ss: a string secondary structure
    :param input_file_type: the type of input file given (fasta or stockholm)
    :param verbose: whether to print verbose logging

    :returns: a tuple containing the input sequence or sequences and BayesPairing output
    '''
    # load all arguments if they were provided or set the default value
    samplesize = arguments.get(
        "sample_size", default=20000, type=int)
    theta = arguments.get("theta", default=1, type=int)
    Lambda = arguments.get("lambda", default=0.35, type=float)
    w = arguments.get("window_length", default=200, type=int)
    s = arguments.get("step_size", default=100, type=int)
    constraints = arguments.get(
        "constraints", default="", type=str)

    # note: these arguments are not expected to be received from the user and are here for future expendability
    Theta = arguments.get("theta", default=1, type=int)
    Delta = arguments.get("delta")
    fuzzy = arguments.get("fuzzy")

    pretrained = False
    interm = False
    first_run = False

    if input_file_type in STOCKHOLM_EXTENSIONS:
        if verbose:
            print("Alignment file detected, scanning alignment", input)
        prediction_scores = {}
        sequences = []
        with open(input, "r") as f:
            for num, record in enumerate(SeqIO.parse(f, "stockholm")):
                if verbose:
                    print("scanning record", record.id)
                id = record.id
                seq = str(record.seq)
                seq = seq.upper()
                if "T" in seq:
                    seq = str(seq).replace("T", "U")
                gseq = seq
                ugseq = gseq.replace("-", "")
                sequences.append((seq, ugseq))
            # print("SEQUENCES",len(sequences))

            if len(sequences[0]) < 250:
                maxs = run_BP(sequences, ss, modules_to_parse, dataset, "NONE", aln=aln, t=t, samplesize=samplesize,
                              pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=first_run)
                first_run = False
                prediction_scores[id] = maxs
                if interm:
                    print(maxs)

            else:
                all_maxes = {}
                index = 0
                while index + w < len(seq):
                    if verbose:
                        print("Running BayesPairing on sequence window:",
                              index, index + w)
                    bf = max(0, index)

                    if len(ss) > (index+w):
                        ss1 = ss[index:index+w]
                    else:
                        ss1 = ""

                    cutseq = [(seq[0][index:index + w], str(seq[0]
                                                            [index:index + w]).replace("-", "")) for seq in sequences]
                    maxs = run_BP(cutseq, ss1, modules_to_parse, dataset, "NONE", aln=aln, t=t, samplesize=samplesize,
                                  pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=first_run)
                    first_run = False
                    for mod in maxs:
                        for cand in maxs[mod]:
                            cand[1] = [[k + bf for k in l] for l in cand[1]]
                    # all_maxes.append(maxs)
                    all_maxes = parse_sequences.aggregate(maxs, all_maxes)
                    index = index + s
                if verbose:
                    print("Running BayesPairing on sequence window:",
                          index, len(seq))
                if len(ss) > (index+w):
                    ss2 = ss[index:]
                else:
                    ss2 = ""
                cutseq = [(seq[0][index:], str(seq[0][index:]).replace("-", ""))
                          for seq in sequences]
                maxs = run_BP(cutseq, ss2, modules_to_parse, dataset, "NONE", aln=aln, t=t, samplesize=samplesize,
                              pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=False)
                if interm:
                    print(maxs)
                for mod in maxs:
                    for cand in maxs[mod]:
                        cand[1] = [[k + bf for k in l] for l in cand[1]]
                # all_maxes.append(maxs)
                all_maxes = parse_sequences.aggregate(maxs, all_maxes)
                prediction_scores[id] = all_maxes
    elif (input_file_type in FASTA_EXTENSIONS) and aln:
        if verbose:
            print("Alignment file detected, scanning alignment", input)
        prediction_scores = {}
        sequences = []
        with open(input, "r") as f:
            for num, record in enumerate(SeqIO.parse(f, "fasta")):
                if verbose:
                    print("scanning record", record.id)
                id = record.id
                seq = str(record.seq)
                seq = seq.upper()
                if "T" in seq:
                    seq = str(seq).replace("T", "U")
                gseq = seq
                ugseq = gseq.replace("-", "")
                sequences.append((seq, ugseq))
            # print("SEQUENCES",len(sequences))

            if len(sequences[0]) < 250:
                maxs = run_BP(sequences, ss, modules_to_parse, dataset, "NONE", aln=aln, t=t, samplesize=samplesize,
                              pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=first_run)
                first_run = False
                prediction_scores[id] = maxs
                if interm:
                    print(maxs)

            else:
                all_maxes = {}
                index = 0
                while index + w < len(seq):
                    if verbose:
                        print("Running BayesPairing on sequence window:",
                              index, index + w)
                    bf = max(0, index)

                    if len(ss) > (index+w):
                        ss1 = ss[index:index+w]
                    else:
                        ss1 = ""

                    cutseq = [(seq[0][index:index + w], str(seq[0]
                                                            [index:index + w]).replace("-", "")) for seq in sequences]
                    maxs = run_BP(cutseq, ss1, modules_to_parse, dataset, "NONE", aln=aln, t=t, samplesize=samplesize,
                                  pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=first_run)
                    first_run = False
                    for mod in maxs:
                        for cand in maxs[mod]:
                            cand[1] = [[k + bf for k in l] for l in cand[1]]
                    # all_maxes.append(maxs)
                    all_maxes = parse_sequences.aggregate(maxs, all_maxes)
                    index = index + s
                if verbose:
                    print("Running BayesPairing on sequence window:",
                          index, len(seq))
                if len(ss) > (index+w):
                    ss2 = ss[index:]
                else:
                    ss2 = ""

                cutseq = [(seq[0][index:], str(seq[0][index:]).replace("-", ""))
                          for seq in sequences]
                maxs = run_BP(cutseq, ss2, modules_to_parse, dataset, "NONE", aln=aln, t=t, samplesize=samplesize,
                              pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=False)
                if interm:
                    print(maxs)
                for mod in maxs:
                    for cand in maxs[mod]:
                        cand[1] = [[k + bf for k in l] for l in cand[1]]
                # all_maxes.append(maxs)
                all_maxes = parse_sequences.aggregate(maxs, all_maxes)
                prediction_scores[id] = all_maxes

    elif (input_file_type in FASTA_EXTENSIONS) and ss == "":
        if verbose:
            print("FASTA file detected, scanning", input)
        prediction_scores = {}
        sequences = []
        with open(input, "r") as f:
            for num, record in enumerate(SeqIO.parse(f, "fasta")):
                id = record.id
                seq = str(record.seq)
                seq = seq.upper()
                sequences.append(seq)
                if "T" in seq:
                    seq = str(seq).replace("T", "U")
                if len(seq) < 250:
                    maxs = run_BP(seq, ss, modules_to_parse, dataset, "NONE", aln=aln, t=t, samplesize=samplesize, pretrained=pretrained,
                                  Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=first_run)
                    first_run = False
                    if interm:
                        print(maxs)
                    prediction_scores[id] = maxs

                else:
                    all_maxes = {}
                    index = 0
                    while index + w < len(seq):
                        if verbose:
                            print(
                                "Running BayesPairing on sequence window:", index, index + w)
                        bf = max(0, index)
                        maxs = run_BP(seq[index:index + w], ss, modules_to_parse, dataset, "NONE", aln=aln, t=t, samplesize=samplesize,
                                      pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=first_run)
                        first_run = False
                        if interm:
                            print(maxs)
                        for mod in maxs:
                            for cand in maxs[mod]:
                                cand[1] = [[k + bf for k in l]
                                           for l in cand[1]]
                        # all_maxes.append(maxs)
                        all_maxes = parse_sequences.aggregate(maxs, all_maxes)
                        index = index + s
                    if verbose:
                        print("Running BayesPairing on sequence window:",
                              index, len(seq))
                    maxs = run_BP(seq[index:], ss, modules_to_parse, dataset, "NONE", aln=aln, t=t, samplesize=samplesize,
                                  pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=False)
                    if interm:
                        print(maxs)
                    for mod in maxs:
                        for cand in maxs[mod]:
                            cand[1] = [[k + bf for k in l] for l in cand[1]]
                    # all_maxes.append(maxs)
                    all_maxes = parse_sequences.aggregate(maxs, all_maxes)
                    prediction_scores[id] = all_maxes

    elif (input_file_type in FASTA_EXTENSIONS) and ss == "infile":
        if verbose:
            print("FASTA file file with secondary structure detected, scanning", input)
        prediction_scores = {}
        sequences = []
        with open(input, "r") as f:
            lines = f.readlines()
            for line_n in range(0, len(lines), 3):
                id = lines[line_n].strip(">").strip("\n")
                seq = lines[line_n+1].strip("\n").strip()
                seq = seq.upper()
                ss = lines[line_n+2].strip("\n").strip()
                sequences.append(seq)
                if "T" in seq:
                    seq = str(seq).replace("T", "U")
                if len(seq) < 300:
                    maxs = run_BP(seq, ss, modules_to_parse, dataset, "NONE", aln=aln, t=t, samplesize=samplesize, pretrained=pretrained,
                                  Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=first_run)
                    first_run = False
                    prediction_scores[id] = maxs

                else:
                    all_maxes = {}
                    index = 0
                    while index + w < len(seq):
                        if verbose:
                            print(
                                "Running BayesPairing on sequence window:", index, index + w)
                        bf = max(0, index)
                        maxs = run_BP(seq[index:index + w], ss, modules_to_parse, dataset, "NONE", aln=aln, t=t, samplesize=samplesize,
                                      pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=first_run)
                        first_run = False
                        if interm:
                            print(maxs)
                        for mod in maxs:
                            for cand in maxs[mod]:
                                cand[1] = [[k + bf for k in l]
                                           for l in cand[1]]
                        # all_maxes.append(maxs)
                        all_maxes = parse_sequences.aggregate(maxs, all_maxes)
                        index = index + s
                    if verbose:
                        print("Running BayesPairing on sequence window:",
                              index, index + w)
                    maxs = run_BP(seq[index:], ss, modules_to_parse, dataset, "NONE", aln=aln, t=t, samplesize=samplesize,
                                  pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=False)
                    if interm:
                        print(maxs)
                    for mod in maxs:
                        for cand in maxs[mod]:
                            cand[1] = [[k + bf for k in l] for l in cand[1]]
                    # all_maxes.append(maxs)
                    all_maxes = parse_sequences.aggregate(maxs, all_maxes)
                    prediction_scores[id] = all_maxes
    else:
        if verbose:
            print("Scanning input sequence", input)
        if "T" in input:
            input = input.upper()
            input = str(input).replace("T", "U")
        sequences = [input]
        if len(input) <= w:
            maxs = run_BP(input, ss, modules_to_parse, dataset, "NONE", aln=aln, t=t, samplesize=samplesize, pretrained=pretrained,
                          Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=first_run, constraints=constraints)
            # first_run=False
            if interm:
                print(maxs)
            prediction_scores = {"input_seq": maxs}
        else:
            all_maxes = {}
            index = 0
            while index + w < len(input):
                if verbose:
                    print("Running Bayespairing on sequence window:",
                          index, index + w)
                #print(input[index:index + w])
                bf = max(0, index)
                maxs = run_BP(input[index:index + w], ss, modules_to_parse, dataset, "NONE", aln=aln, t=t, samplesize=samplesize,
                              pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=first_run)
                first_run = False
                if interm:
                    print(maxs)
                for mod in maxs:
                    for cand in maxs[mod]:
                        cand[1] = [[k + bf for k in l] for l in cand[1]]
                # all_maxes.append(maxs)
                all_maxes = parse_sequences.aggregate(maxs, all_maxes)
                index = index + s
            if verbose:
                print("Running Bayespairing on sequence window:", index, len(input))
            maxs = run_BP(input[index:], ss, modules_to_parse, dataset, "NONE", aln=aln, t=t, samplesize=samplesize,
                          pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=False)
            if interm:
                print(maxs)
            for mod in maxs:
                for cand in maxs[mod]:
                    cand[1] = [[k + bf for k in l] for l in cand[1]]
            # all_maxes.append(maxs)
            all_maxes = parse_sequences.aggregate(maxs, all_maxes)

    return sequences, prediction_scores
