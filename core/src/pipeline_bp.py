import os
import pickle
from multiprocessing import Process, Manager
from . import testSS
from . import BayesPairing
from .parse_sequences import unpick, run_BP, process_siblings, parse_sequences, check_already_there, seq_ranges
import random
from Bio import SeqIO
import argparse
import time as timer
import collections
from functools import reduce
from Bio import AlignIO
import os.path
from operator import itemgetter
import sys

CURRENT_DIRECTORY = os.path.dirname(__file__)
INPUT_DIRECTORY = os.path.join(CURRENT_DIRECTORY, "../input")

def run_fasta(input, modules_to_parse, dataset, ss="", arguments={}):
    if "aln" in arguments:
        aln = arguments["aln"]
    else:
        aln = False

    if "theta" in arguments:
        Theta = float(arguments["theta"])
    else:
        Theta = 1

    if "Lambda" in arguments:
        Lambda = float(arguments["Lambda"])
    else:
        Lambda = 0.35

    if "delta" in arguments:
        Delta = float(arguments["delta"])
    else:
        Delta = None

    if "fuzzy" in arguments:
        fuzzy = arguments["fuzzy"]
    else:
        fuzzy = None


    if "verbose" in arguments:
        verbose = arguments["verbose"]
    else:
        verbose = False

    if "pretrained" in arguments:
        pretrained = arguments["pretrained"]
    else:
        pretrained = False

    if "samplesize" in arguments:
        samplesize = int(arguments["samplesize"])
    else:
        samplesize = 20000

    if "t" in arguments:
        t = float(arguments["t"])
    else:
        t = -2.3
    if "w" in arguments:
        w = int(arguments["w"])
    else:
        w = 200
    if "s" in arguments:
        s = int(arguments["s"])
    else:
        s = 100
    if "o" in arguments:
        o = arguments["o"]
    else:
        o = "output"


    if "interm" in arguments:
         interm = arguments["interm"]
    else:
         interm = False

            
    if "constraints" in arguments:
         constraints = arguments["constraints"]
    else:
         constraints = ""


    if 'init' in arguments:
        first_run=arguments["init"]
    else:
        first_run= False
        
        
    if ".st" in input.lower():
        if verbose:
            print("Alignment file detected, scanning alignment", input)
        prediction_scores = {}
        sequences = []
        with open(os.path.join(INPUT_DIRECTORY, input), "r") as f:
            for num, record in enumerate(SeqIO.parse(f, "stockholm")):
                if verbose:
                    print("scanning record", record.id)
                id = record.id
                seq = str(record.seq)
                seq = seq.upper()
                if "T" in seq:
                    seq = str(seq).replace("T", "U")
                gseq = seq
                ugseq= gseq.replace("-","")
                sequences.append((seq,ugseq))
            #print("SEQUENCES",len(sequences))

            if len(sequences[0]) < 3000:
                maxs = run_BP(sequences, ss, modules_to_parse, dataset, "NONE", aln= aln, t= t, samplesize=samplesize, pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=first_run)
                first_run=False
                prediction_scores[id] = [maxs]
                if interm:
                    print(maxs)

            else:
                all_maxes = []
                index = 0
                while index + w < len(seq):
                    if verbose:
                        print("Running BayesPairing on sequence window:", index, index + w)
                    bf = max(0, index)
                    
                    if len(ss)>(index+w):
                        ss1 = ss[index:index+w]
                    else:
                        ss1 = ""
                        
                    cutseq = [(seq[0][index:index + w],str(seq[0][index:index + w]).replace("-","")) for seq in sequences]
                    maxs = run_BP(cutseq, ss1, modules_to_parse, dataset, "NONE", aln= aln, t= t, samplesize=samplesize, pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=first_run)
                    first_run=False
                    for mod in maxs:
                        for cand in maxs[mod]:
                            cand[1] = [[k + bf for k in l] for l in cand[1]]
                    all_maxes.append(maxs)
                    index = index + s
                if verbose:
                    print("Running BayesPairing on sequence window:", index, len(seq))
                if len(ss)>(index+w):
                    ss2 = ss[index:]
                else:
                    ss2 = ""
                cutseq = [(seq[0][index:],str(seq[0][index:]).replace("-","")) for seq in sequences]
                maxs = run_BP(cutseq, ss2, modules_to_parse, dataset, "NONE", aln= aln, t= t, samplesize=samplesize, pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=False)
                if interm:
                    print(maxs)
                for mod in maxs:
                    for cand in maxs[mod]:
                        cand[1] = [[k + bf for k in l] for l in cand[1]]
                all_maxes.append(maxs)
                prediction_scores[id] = all_maxes
        # print("PREDICTION_SCORES",prediction_scores)
        return prediction_scores


    
    elif ".fa" in input.lower() and aln:
        if verbose:
            print("Alignment file detected, scanning alignment", input)
        prediction_scores = {}
        sequences = []
        with open(os.path.join(INPUT_DIRECTORY, input), "r") as f:
            for num, record in enumerate(SeqIO.parse(f, "fasta")):
                if verbose:
                    print("scanning record", record.id)
                id = record.id
                seq = str(record.seq)
                seq = seq.upper()
                if "T" in seq:
                    seq = str(seq).replace("T", "U")
                gseq = seq
                ugseq= gseq.replace("-","")
                sequences.append((seq,ugseq))
            #print("SEQUENCES",len(sequences))

            if len(sequences[0]) < 3000:
                maxs = run_BP(sequences, ss, modules_to_parse, dataset, "NONE", aln= aln, t= t, samplesize=samplesize, pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=first_run)
                first_run=False
                prediction_scores[id] = [maxs]
                if interm:
                    print(maxs)

            else:
                all_maxes = []
                index = 0
                while index + w < len(seq):
                    if verbose:
                        print("Running BayesPairing on sequence window:", index, index + w)
                    bf = max(0, index)
                    
                    if len(ss)>(index+w):
                        ss1 = ss[index:index+w]
                    else:
                        ss1 = ""
                        
                    cutseq = [(seq[0][index:index + w],str(seq[0][index:index + w]).replace("-","")) for seq in sequences]
                    maxs = run_BP(cutseq, ss1, modules_to_parse, dataset, "NONE", aln= aln, t= t, samplesize=samplesize, pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=first_run)
                    first_run=False
                    for mod in maxs:
                        for cand in maxs[mod]:
                            cand[1] = [[k + bf for k in l] for l in cand[1]]
                    all_maxes.append(maxs)
                    index = index + s
                if verbose:
                    print("Running BayesPairing on sequence window:", index, len(seq))
                if len(ss)>(index+w):
                    ss2 = ss[index:]
                else:
                    ss2 = ""
                
                cutseq = [(seq[0][index:],str(seq[0][index:]).replace("-","")) for seq in sequences]
                maxs = run_BP(cutseq, ss2, modules_to_parse, dataset, "NONE", aln= aln, t= t, samplesize=samplesize, pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=False)
                if interm:
                    print(maxs)
                for mod in maxs:
                    for cand in maxs[mod]:
                        cand[1] = [[k + bf for k in l] for l in cand[1]]
                all_maxes.append(maxs)
                prediction_scores[id] = all_maxes
        return prediction_scores


    elif ".fa" in input.lower() and ss=="":
        if verbose:
            print("FASTA file detected, scanning", input)
        prediction_scores = {}
        sequences = []
        with open(os.path.join(INPUT_DIRECTORY, input), "rU") as f:
            for num, record in enumerate(SeqIO.parse(f, "fasta")):
                id = record.id
                seq = str(record.seq)
                seq = seq.upper()
                sequences.append(seq)
                if "T" in seq:
                    seq = str(seq).replace("T", "U")
                if len(seq)<300:
                    maxs = run_BP(seq, ss, modules_to_parse, dataset, "NONE", aln= aln, t= t, samplesize=samplesize, pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=first_run)
                    first_run=False
                    if interm:
                        print(maxs)
                    prediction_scores[id] = [maxs]

                else:
                    all_maxes = []
                    index = 0
                    while index + w < len(seq):
                        if verbose:
                            print("Running BayesPairing on sequence window:", index, index + w)
                        bf = max(0, index)
                        maxs = run_BP(seq[index:index + w], ss, modules_to_parse, dataset, "NONE", aln= aln, t= t, samplesize=samplesize, pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=first_run)
                        first_run=False
                        if interm:
                            print(maxs)
                        for mod in maxs:
                            for cand in maxs[mod]:
                                cand[1] = [[k + bf for k in l] for l in cand[1]]
                        all_maxes.append(maxs)
                        index = index + s
                    if verbose:
                        print("Running BayesPairing on sequence window:", index, len(seq))
                    maxs = run_BP(seq[index:], ss, modules_to_parse, dataset, "NONE", aln= aln, t= t, samplesize=samplesize, pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=False)
                    if interm:
                        print(maxs)
                    for mod in maxs:
                        for cand in maxs[mod]:
                            cand[1] = [[k + bf for k in l] for l in cand[1]]
                    all_maxes.append(maxs)
                    prediction_scores[id] = all_maxes

        return prediction_scores

    elif "fa" in input.lower() and ss=="infile":
        if verbose:
            print("FASTA file file with secondary structure detected, scanning", input)
        prediction_scores = {}
        sequences = []
        with open(os.path.join(INPUT_DIRECTORY, input), "rU") as f:
            lines = f.readlines()
            for line_n in range(0, len(lines), 3):
                id = lines[line_n].strip(">").strip("\n")
                seq = lines[line_n+1].strip("\n").strip()
                seq = seq.upper()
                ss = lines[line_n+2].strip("\n").strip()
                sequences.append(seq)
                if "T" in seq:
                    seq = str(seq).replace("T", "U")
                if len(seq)<300:
                    maxs = run_BP(seq, ss, modules_to_parse, dataset, "NONE", aln= aln, t= t, samplesize=samplesize, pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=first_run)
                    first_run=False
                    prediction_scores[id] = [maxs]

                else:
                    all_maxes = []
                    index = 0
                    while index + w < len(seq):
                        if verbose:
                            print("Running BayesPairing on sequence window:", index, index + w)                        
                        bf = max(0, index)
                        maxs = run_BP(seq[index:index + w], ss, modules_to_parse, dataset, "NONE", aln= aln, t= t, samplesize=samplesize, pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=first_run)
                        first_run=False
                        if interm:
                            print(maxs)
                        for mod in maxs:
                            for cand in maxs[mod]:
                                cand[1] = [[k + bf for k in l] for l in cand[1]]
                        all_maxes.append(maxs)
                        index = index + s
                    if verbose:
                        print("Running BayesPairing on sequence window:", index, index + w)                   
                    maxs = run_BP(seq[index:], ss, modules_to_parse, dataset, "NONE", aln= aln, t= t, samplesize=samplesize, pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=False)
                    if interm:
                        print(maxs)
                    for mod in maxs:
                        for cand in maxs[mod]:
                            cand[1] = [[k + bf for k in l] for l in cand[1]]
                    all_maxes.append(maxs)
                    prediction_scores[id] = all_maxes
        return prediction_scores

    else:
        if verbose:
            print("Scanning input sequence", input)
        if "T" in input:
            input = input.upper()
            input = str(input).replace("T", "U")
        sequences = [input]
        if len(input) <= w:
            maxs = run_BP(input, ss, modules_to_parse, dataset, "NONE", aln= aln, t= t, samplesize=samplesize, pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=first_run, constraints=constraints)
            #first_run=False
            if interm:
                print(maxs)
            return maxs
        else:
            all_maxes = []
            index = 0
            while index + w < len(input):
                if verbose:
                    print("Running Bayespairing on sequence window:", index, index + w)
                #print(input[index:index + w])
                bf = max(0, index)
                maxs = run_BP(input[index:index + w], ss, modules_to_parse, dataset, "NONE", aln= aln, t= t, samplesize=samplesize, pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=first_run)
                first_run=False
                if interm:
                    print(maxs)
                for mod in maxs:
                    for cand in maxs[mod]:
                        cand[1] = [[k + bf for k in l] for l in cand[1]]
                all_maxes.append(maxs)
                index = index + s
            if verbose:
                print("Running Bayespairing on sequence window:", index, len(input))
            maxs = run_BP(input[index:], ss, modules_to_parse, dataset, "NONE", aln= aln, t= t, samplesize=samplesize, pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=False)
            if interm:
                print(maxs)
            for mod in maxs:
                for cand in maxs[mod]:
                    cand[1] = [[k + bf for k in l] for l in cand[1]]
            all_maxes.append(maxs)
            return all_maxes

    return []