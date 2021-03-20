import os
import pickle
from multiprocessing import Process, Manager
from . import testSS
from . import BayesPairing
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
import json
from .chefschoice import bp_chefs_choice

CURRENT_DIRECTORY = os.path.dirname(__file__)
INPUT_DIRECTORY = os.path.join(CURRENT_DIRECTORY, "../input")

def unpick(dataset,direc,typ):
    file_path = os.path.join(CURRENT_DIRECTORY, "../"+direc+"/" + dataset + "_"+typ)
    if os.path.exists(file_path):
        print("Dataset found!", file=sys.stderr)
        while True:
            try:
                nets = pickle.load(open(file_path, "rb"))
                break
            except:
                pass
        return nets
    else:
        print("Dataset not found!", file=sys.stderr)
        return None

def run_BP(seq, ss, modules_to_parse, dataset, left_out, aln=False, t=-5, samplesize=20000, pretrained=False, Lambda=0.35, Theta=1,Delta=None, leave_out_sequence=False,left_out_sequence="",fuzzy=False,verbose=False, first_run=False, constraints = '',indexes=[]):
    file_path = os.path.join(CURRENT_DIRECTORY, "../models/" + dataset + "_models.pickle")
    if not (os.path.isfile(file_path)) or left_out!="NONE":
        if verbose:
            print("running cross-validation")
        nets = collections.OrderedDict()
    else:
        #nets = pickle.load(open("../models/" + dataset + "_models.pickle", "rb"))
        nets = unpick(dataset,"models","models.pickle")
    if not pretrained:
        if verbose and first_run:
            print("reseting models, training with Theta",str(Theta))
        if first_run:
            nets = collections.OrderedDict()
        #print("SETTING UP BAYES NET")
        BNs = BayesPairing.setup_BN(modules_to_parse, nets, dataset, left_out,leave_out_sequence=leave_out_sequence,left_out_sequence=left_out_sequence, Lambda=Lambda,Theta=Theta, verbose=verbose,indexes=indexes)
        pickle.dump(BNs,open(os.path.join(CURRENT_DIRECTORY, "../models/" + dataset + "_models.pickle"), "wb"))
    else:
        #nets = pickle.load(open("../models/" + dataset + "_models.pickle", "rb"))
        nets = unpick(dataset,"models","models.pickle")
        if verbose:
            print('PRE-LOADED MODELS:',str(list(nets.keys())))
        BNs = nets
    #print("PARSING SEQUENCE")
    if aln:
        return_dict = BayesPairing.parse_alignment(seq, modules_to_parse, ss, dataset, BNs, t, samplesize, Lambda, Theta, Delta, fuzzy=fuzzy, verbose=verbose)
    else:
        return_dict = BayesPairing.parse_sequence(seq, modules_to_parse, ss, dataset, BNs,t, samplesize, Lambda, Theta, Delta, fuzzy=fuzzy, verbose=verbose, constraints = constraints)
    #siblings = pickle.load(open("../models/" + dataset + "_siblings.pickle", "rb"))    
    #siblings = unpick(dataset,"models","siblings.pickle")

    # load the siblings from the json dataset, siblings should be a dict in the format {0: [1,2,3], 1:[0,2,3]}...
    with open(os.path.join(CURRENT_DIRECTORY, "../models/" + dataset + ".json")) as f:
        module_siblings = json.load(f)

    # map the json dict to proper format
    siblings = {int(k): v["siblings"] for k, v in module_siblings.items()}

    # instead of mapping the dict to the sibling dict, pass the dataset through and process after
    noSiblings = process_siblings(return_dict,siblings)
    
    return noSiblings

def process_siblings(results,siblings):
    real_results = {}
    
    
    for mod in results:
        if len(results[mod])==0:
            real_results[mod] = results[mod]
            continue
        
        if mod in siblings:
            sibs = siblings[mod]
        else:
            real_results[mod] = results[mod]
            continue
        #print("RESULTS",results[mod])
        bestCall = results[mod][0][2]
        sibRes = []
        for x in sibs:
            if x in results:
                if len(results[x])>0:
                    sibRes.append(results[x][0][2])
        if not any(x>bestCall+10 for x in sibRes):
            real_results[mod] = results[mod]
    return real_results

def parse_sequences(input,modules_to_check=[],dataset="",ss="",m=4,n=3,sm=0.3,mc=3,p=25000,sw=1,t=15.7,w=200,s=100,sscons=2):
    if dataset=="":
        dataset="rna3dmotif"
    graphs = pickle.load(open(os.path.join(CURRENT_DIRECTORY, "../models/" + dataset + "_one_of_each_graph.cPickle"), "rb"))
    if len(modules_to_check) == 0:
        modules_to_check = range(len(graphs))

    argum = {"m":m, "n":n, "sm":sm,"mc":mc,"p":p,"sw":sw,"t":t,"w":w,"s":s,"sscons":2}
    this_output = run_fasta(input,modules_to_check,dataset,ss,argum)
    return this_output



#to debug
def get_stats(prediction_scores,modules_to_parse,threshold=-5):
    #print("getting stats",prediction_scores)
    n_sequences = 0
    n_hits = {}
    for m in modules_to_parse:
        n_hits[m]=0
    for sequence in prediction_scores:
        hit_dict = {}
        n_sequences=n_sequences+1
        #for window in prediction_scores[sequence]:
        window = prediction_scores[sequence] 
        scored_positions = {}
        #print(window)
        for module in window:
            scored_positions[module] = []
            #print("SCANNING",module,window[module])
            

            mod_number = module
            if len(window[module])>0:
                max_score = sorted(window[module],key=itemgetter(2),reverse=True)[0][2]
                if max_score>threshold and window[module][0][1] not in scored_positions[module]:
                    #print("admissible score")
                    scored_positions[module].append(window[module][0][1])
                    hit_dict[module]=True
        for mod in hit_dict:
            if hit_dict[mod]:
                n_hits[mod]+=1
    OUTPUT_STRING = ""
    output = []
    #print("HITS",n_hits)
    for m in n_hits:
        avgHits = round(n_hits[m]/n_sequences,2)
        if avgHits >0 :
            output.append(["|", m, n_hits[m],avgHits, "|"])
    OUTPUT_STRING=OUTPUT_STRING+("=========================================\n")
    headers = ["|", "MODULE", "N HITS", "PERCENTAGE", "|"]
    output.insert(0, headers)
    for row in output:
        OUTPUT_STRING=OUTPUT_STRING+("{: >1}\t{: >4}\t{: >6}\t{: >10}\t{: >1}".format(*row))+"\n"
    OUTPUT_STRING=OUTPUT_STRING+"=========================================\n"
    return OUTPUT_STRING

def aggregate(maxs,all_maxes):
    #print("ALL",all_maxes)
    #print("maxs",maxs)
    for mod in maxs:
        if mod not in all_maxes:
            all_maxes[mod] = maxs[mod]
        else:
            for j in maxs[mod]:
                all_maxes[mod].append(j)
    return all_maxes

def run_fasta(input, modules_to_parse, dataset, ss="", arguments={}):
    fOUTPUT = ""

    if not os.path.exists(os.path.join(CURRENT_DIRECTORY, '../output')):
        os.makedirs(os.path.join(CURRENT_DIRECTORY, '../output'))

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

            if len(sequences[0]) < 250:
                maxs = run_BP(sequences, ss, modules_to_parse, dataset, "NONE", aln= aln, t= t, samplesize=samplesize, pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=first_run)
                first_run=False
                prediction_scores[id] = maxs
                if interm:
                    print(maxs)

            else:
                all_maxes = {}
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
                    #all_maxes.append(maxs)
                    all_maxes = aggregate(maxs,all_maxes)
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
                #all_maxes.append(maxs)
                all_maxes = aggregate(maxs,all_maxes)
                prediction_scores[id] = all_maxes
        # print("PREDICTION_SCORES",prediction_scores)
        for id in prediction_scores:
            fOUTPUT = fOUTPUT + "\nRESULTS FOR ID " + id + "\n"
            fOUTPUT = fOUTPUT + present_output(prediction_scores[id], t)
        stats = get_stats(prediction_scores, modules_to_parse, t)
        fOUTPUT = fOUTPUT + "\n" + stats
        pickle.dump(prediction_scores, open(os.path.join(CURRENT_DIRECTORY, "../output/" + o + ".pickle"), "wb"))


    
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

            if len(sequences[0]) < 250:
                maxs = run_BP(sequences, ss, modules_to_parse, dataset, "NONE", aln= aln, t= t, samplesize=samplesize, pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=first_run)
                first_run=False
                prediction_scores[id] = maxs
                if interm:
                    print(maxs)

            else:
                all_maxes = {}
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
                    #all_maxes.append(maxs)
                    all_maxes = aggregate(maxs,all_maxes)
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
                #all_maxes.append(maxs)
                all_maxes = aggregate(maxs,all_maxes)
                prediction_scores[id] = all_maxes
        # print("PREDICTION_SCORES",prediction_scores)
        for id in prediction_scores:
            fOUTPUT = fOUTPUT + "\nRESULTS FOR ID " + id + "\n"
            fOUTPUT = fOUTPUT + present_output(prediction_scores[id], t)
        stats = get_stats(prediction_scores, modules_to_parse, t)
        fOUTPUT = fOUTPUT + "\n" + stats
        pickle.dump(prediction_scores, open(os.path.join(CURRENT_DIRECTORY, "../output/" + o + ".pickle"), "wb"))


    elif ".fa" in input.lower() and ss=="":
        if verbose:
            print("FASTA file detected, scanning", input)
        prediction_scores = {}
        sequences = []
        with open(os.path.join(INPUT_DIRECTORY, input), "r") as f:
            for num, record in enumerate(SeqIO.parse(f, "fasta")):
                id = record.id
                seq = str(record.seq)
                seq = seq.upper()
                sequences.append(seq)
                if "T" in seq:
                    seq = str(seq).replace("T", "U")
                if len(seq)<250:
                    maxs = run_BP(seq, ss, modules_to_parse, dataset, "NONE", aln= aln, t= t, samplesize=samplesize, pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=first_run)
                    first_run=False
                    if interm:
                        print(maxs)
                    prediction_scores[id] = maxs

                else:
                    all_maxes = {}
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
                        #all_maxes.append(maxs)
                        all_maxes = aggregate(maxs,all_maxes)
                        index = index + s
                    if verbose:
                        print("Running BayesPairing on sequence window:", index, len(seq))
                    maxs = run_BP(seq[index:], ss, modules_to_parse, dataset, "NONE", aln= aln, t= t, samplesize=samplesize, pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=False)
                    if interm:
                        print(maxs)
                    for mod in maxs:
                        for cand in maxs[mod]:
                            cand[1] = [[k + bf for k in l] for l in cand[1]]
                    #all_maxes.append(maxs)
                    all_maxes = aggregate(maxs,all_maxes)
                    prediction_scores[id] = all_maxes

        for id in prediction_scores:
            fOUTPUT=fOUTPUT+"\nRESULTS FOR ID "+id+"\n"
            fOUTPUT=fOUTPUT+present_output(prediction_scores[id], t)
        stats=get_stats(prediction_scores,modules_to_parse,t)
        fOUTPUT= fOUTPUT+"\n"+stats
        pickle.dump(prediction_scores,open(os.path.join(CURRENT_DIRECTORY, "../output/"+o+".pickle"),"wb"))

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
                    prediction_scores[id] = maxs

                else:
                    all_maxes = {}
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
                        #all_maxes.append(maxs)
                        all_maxes = aggregate(maxs,all_maxes)
                        index = index + s
                    if verbose:
                        print("Running BayesPairing on sequence window:", index, index + w)                   
                    maxs = run_BP(seq[index:], ss, modules_to_parse, dataset, "NONE", aln= aln, t= t, samplesize=samplesize, pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=False)
                    if interm:
                        print(maxs)
                    for mod in maxs:
                        for cand in maxs[mod]:
                            cand[1] = [[k + bf for k in l] for l in cand[1]]
                    #all_maxes.append(maxs)
                    all_maxes = aggregate(maxs,all_maxes)
                    prediction_scores[id] = all_maxes
        for id in prediction_scores:
            fOUTPUT=fOUTPUT+"\nRESULTS FOR ID "+id+"\n"
            fOUTPUT=fOUTPUT+present_output(prediction_scores[id], t)+"\n"
        stats=get_stats(prediction_scores,modules_to_parse,t)
        fOUTPUT= fOUTPUT+"\n"+stats
        pickle.dump(prediction_scores,open(os.path.join(CURRENT_DIRECTORY, "../output/"+o+".pickle"),"wb"))

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
            fOUTPUT=fOUTPUT+present_output(maxs, t)+"\n"
            prediction_scores = {"input_seq":maxs}
            pickle.dump(prediction_scores,open(os.path.join(CURRENT_DIRECTORY, "../output/"+o+".pickle"),"wb"))
        else:
            all_maxes = {}
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
                #all_maxes.append(maxs)
                all_maxes = aggregate(maxs,all_maxes)
                index = index + s
            if verbose:
                print("Running Bayespairing on sequence window:", index, len(input))
            maxs = run_BP(input[index:], ss, modules_to_parse, dataset, "NONE", aln= aln, t= t, samplesize=samplesize, pretrained=pretrained, Lambda=Lambda, Theta=Theta, Delta=Delta, fuzzy=fuzzy, verbose=verbose, first_run=False)
            if interm:
                print(maxs)
            for mod in maxs:
                for cand in maxs[mod]:
                    cand[1] = [[k + bf for k in l] for l in cand[1]]
            #all_maxes.append(maxs)
            all_maxes = aggregate(maxs,all_maxes)

            fOUTPUT=fOUTPUT+present_output(all_maxes, t)+"\n"
            prediction_scores = {"input_seq":all_maxes}
            pickle.dump(prediction_scores,open(os.path.join(CURRENT_DIRECTORY, "../output/"+o+".pickle"),"wb"))
    return fOUTPUT,sequences,prediction_scores



def seq_ranges(all_pos):
    all_pos = reduce(lambda x,y: x+y,all_pos)
    output_string = ""
    for ind, pos in enumerate(all_pos):
        if ind == 0:
            output_string = str(pos)
            continue
        if ind == 1:
            if all_pos[ind - 1] != pos - 1:
                output_string = output_string + "," + str(pos)
                continue
        if all_pos[ind - 1] != pos - 1:
            if all_pos[ind - 1] == (all_pos[ind - 2] + 1):
                output_string = output_string + "-" + str(all_pos[ind - 1]) + "," + str(pos)
            else:
                output_string = output_string + "," + str(pos)
            continue
        if ind == len(all_pos) - 1:
            if all_pos[ind] == (all_pos[ind - 1] + 1):
                output_string = output_string + "-" + str(pos)
            else:
                output_string = output_string + "," + str(pos)
    return output_string


def check_already_there(results, new_range, new_mod):
    new_seq_range = new_range[1]
    for current_results in results:
        old_range = current_results[3]
        old_module = current_results[1]
        if old_range.strip() == new_range[1].strip() and int(old_module) == int(new_mod):
            return True

    return False


def present_output(all_maxes, threshold, offset=0):
    OUTPUT_STRING = ""
    output = []
    #for m in all_maxes:
        #print("all_maxes",m)
    m = all_maxes
    for current_module in sorted(m.keys()):
        if len(m[current_module])<1:
            continue
        for sub_max in m[current_module]:
            this_max = [round(sub_max[2], 3), seq_ranges(sub_max[1]),sub_max[0]]
            if this_max[0] > threshold and not check_already_there(output, this_max, current_module):
                output.append(["|", current_module, *this_max, "|"])
    output = sorted(output)
    #print(output)
    OUTPUT_STRING=OUTPUT_STRING+("=========================================================================================\n")
    headers = ["|", "MODULE", "SCORE", "MODULE POSITIONS", "SEQUENCE", "|"]
    output.insert(0, headers)
    for row in output:
        OUTPUT_STRING=OUTPUT_STRING+("{: >1}\t{: >4}\t{: >6}\t{: >30}\t{: >25}\t{: >1}".format(*row))+"\n"
    OUTPUT_STRING=OUTPUT_STRING+"=========================================================================================\n"
    return OUTPUT_STRING


if __name__ == "__main__":
    # Manages arguments
    START_TIME = timer.time()

    arguments = {}
    parser = argparse.ArgumentParser()
    parser.add_argument("--verbose", help="increase output verbosity", action="store_true")
    #parser.add_argument("-sm",help="sample more minor sequences by making the regex more flexible (include more nucleotides[0.1 to 0.9,0.9 being the most strict]")
    #parser.add_argument("-mc", help="Allow more substitutions for the regex [1 to 3]")
    #parser.add_argument("-m", help="number of folded candidates, must be larger than n")
    parser.add_argument("-seq", help="sequences to parse, FASTA file")
    parser.add_argument("-ss", help="secondary structure, as a string for single sequence, or write -ss infile if ss in FASTA")
    parser.add_argument("-d", help="Dataset, as the pickle file name without the extension. Default will be the dataset presented in the paper")
    parser.add_argument("-mod", nargs='+', help="If you only want to parse specific modules, list them. ex -mod 1 4 12")
    parser.add_argument("--aln", help="to input an alignment file (stockholm)", action="store_true")
    parser.add_argument("-t", help="Score threshold for a module to be called. [-5 to 5]. Default:0", default=0)
    parser.add_argument("-w", help="Window Length [50 to 300]. Default:200 ")
    parser.add_argument("-s", help="Step size between windows [10 to 200]. Default:100 ")
    parser.add_argument("-lambda", help="weight of the secondary structure weight(lambda). Default:0.2 ", dest="Lambda")
    parser.add_argument("-o", help="Name of the output file. Default: output ",default="output")
    parser.add_argument("--interm", help="output the best intermediate results.",action="store_true")
    #parser.add_argument("-sscons", help="Constraint level of the module-secondary structure match. Integer from 0 to 5, 0 being most constrained")
    parser.add_argument("--pretrained", help="Use this option if you have already trained all relevant models", action="store_true",dest="pretrained")
    parser.add_argument("-delta", help="If you do not want to sample the whole secondary structure landscape, the maximum delta energy.")
    parser.add_argument("-samplesize", help="Size of the structure sample. Default value:5000", dest="samplesize")
    parser.add_argument("-theta", help="Constant to multiply the sequence score by. Default 1")
    parser.add_argument("-constraints", help="RNAfold constraints on input sequence")
    parser.add_argument("--init", help="If you want to reset models, or this is your first time running BayesPairing", action="store_true")
    
    args = parser.parse_args()
    if args.verbose:
        arguments["verbose"] = True
    if args.pretrained:
        arguments["pretrained"] = True
    if args.delta != None:
        arguments["delta"] = args.delta
    if args.samplesize:
        arguments["samplesize"] = args.samplesize
    if args.t != None:
        arguments["t"] = args.t
    if args.theta != None:
        arguments["theta"] = args.theta
    if args.Lambda != None:
        arguments["lambda"] = args.Lambda
    if args.w != None:
        arguments["w"] = args.w
    if args.s != None:
        arguments["s"] = args.s
    if args.aln:
        arguments["aln"] = True
    else:
        arguments["aln"] = False
    if args.o != None:
        arguments["o"] = args.o
    if args.mod != None:
        arguments["mod"] = args.mod
    if args.interm != None:
        arguments["interm"] = args.interm
    if args.constraints != None:
        arguments["constraints"] = args.constraints    
     
    if args.init != None:
        arguments["init"] = args.init
        
    if args.seq != None:
        seq = args.seq
    else:
        print("No sequence input, please use python parse_sequences.py -seqs SEQUENCES_FASTA_FILE")
        exit()

    if args.ss:
        ss = args.ss
    else:
        ss = ""

    if args.d != None:
        dataset = args.d
        if dataset == "3dMotifAtlas_RELIABLE":
            dataset = RELIABLE
        elif dataset == "3dMotifAtlas_ALL":
            dataset = ALL
    else:
        dataset = "ALL"
    # the default dataset is rna3dmotif

    # we load the modules from the dataset to get the number of modules available.
    #graphs = pickle.load(open("../models/" + dataset + "_one_of_each_graph.cPickle", "rb"))
    #graphs = unpick(dataset,"models","one_of_each_graph.cPickle")

    # updated load from json
    with open(os.path.join(CURRENT_DIRECTORY, "../models/" + dataset +".json")) as f:
        modules = json.load(f)

    if "mod" in arguments:
        modules_to_check = [int(input_number) for input_number in arguments["mod"]]
    else:
        # By default, run all modules.
        #modules_to_check = list(range(len(graphs)))
        #excluded_modules = [23,36,46,62]
        #for i in excluded_modules:
        #    if i in modules_to_check:
        #        modules_to_check.remove(i)
        modules_to_check = range(len(modules))


    #timer.sleep(5)
    # executes BayesPairing on the sequence
    if not arguments["aln"]:
        all_svg_hits = {}
        all_chef_ss = []
        #run BP2 and print results
        toPrint, seqInfo, all_results = run_fasta(seq, modules_to_check, dataset, ss, arguments)
        print(toPrint)

        #generate SVG
        outName = arguments["o"]
        for seqCounter,inputSeqKey in enumerate(list(all_results.keys())):
            if seqCounter>0:
                finalName = outName+str(seqCounter)
            else:
                finalName = outName

            #print("THE HITS")
            #print(all_results[inputSeqKey])
            modules_in_svg, chef_ss = bp_chefs_choice(all_results[inputSeqKey],seqInfo[seqCounter],float(arguments["t"]),finalName)

            #now we need to fill svg hits
            svg_hits = {}
            for hit in modules_in_svg:
                modID, modInfo = hit
                if modID not in svg_hits:
                    svg_hits[modID]=[modInfo]
                else:
                    svg_hits[modID].append(modInfo)
            all_svg_hits[inputSeqKey] = svg_hits
            all_chef_ss.append(chef_ss)
        #enter data in json
        output_dict = {"input": seqInfo, "params": arguments, "chefs_choice_struct": all_chef_ss, "all_hits":all_results, "svg_hits" : all_svg_hits }
    

    else: #if the input is an alignment, then no SVG
        toPrint, seqInfo, all_results = run_fasta(seq, modules_to_check, dataset, ss, arguments)
        print(toPrint)
        output_dict = {"input": seq, "params": arguments, "all_hits":all_results }


    outFileName = os.path.join(CURRENT_DIRECTORY, "../output/"+arguments["o"] +".json")
    with open(outFileName,"w+") as f:
        json.dump(output_dict,f)




    END_TIME = timer.time()
    print("TOTAL TIME:",round(END_TIME-START_TIME,3))

#if __name__ != "__main__":
#    parse_sequences("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
