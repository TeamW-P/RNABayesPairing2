import os
import pickle
import BayesPairing
import random
from Bio import SeqIO
from random import shuffle
import parse_sequences
import random
DATASET_NAME = "default"
graphs = pickle.load(open("../models/default_one_of_each_graph.cPickle", "rb"))
pdbs = pickle.load(open("../models/default_PDB_names.cPickle", "rb"))
mod_sequences,test_sequences, all_indexes = pickle.load(open("../models/default_sequences.pickle", "rb"))

fasta_path = "pdb_seqres.txt"
listd = os.listdir("../Graphs/all_graphs_pickled")

PDBlist = set([x[0:4] for x in listd])


import sys, random


def computeCountAndLists(s):
    # WARNING: Use of function count(s,'UU') returns 1 on word UUU
    # since it apparently counts only nonoverlapping words UU
    # For this reason, we work with the indices.

    # Initialize lists and mono- and dinucleotide dictionaries
    List = {}  # List is a dictionary of lists
    List['A'] = [];
    List['C'] = [];
    List['G'] = [];
    List['U'] = [];
    List['.'] = []
    nuclList = ["A", "C", "G", "U", "."]
    s = s.upper()
    nuclCnt = {}  # empty dictionary
    dinuclCnt = {}  # empty dictionary
    for x in nuclList:
        nuclCnt[x] = 0
        dinuclCnt[x] = {}
        for y in nuclList:
            dinuclCnt[x][y] = 0

    # Compute count and lists
    nuclCnt[s[0]] = 1
    nuclTotal = 1
    dinuclTotal = 0
    for i in range(len(s) - 1):
        x = s[i];
        y = s[i + 1]
        List[x].append(y)
        nuclCnt[y] += 1;
        nuclTotal += 1
        dinuclCnt[x][y] += 1;
        dinuclTotal += 1
    assert (nuclTotal == len(s))
    assert (dinuclTotal == len(s) - 1)
    return nuclCnt, dinuclCnt, List


def chooseEdge(x, dinuclCnt):
    numInList = 0
    for y in ['A', 'C', 'G', 'U']:
        numInList += dinuclCnt[x][y]
    z = random.random()
    denom = dinuclCnt[x]['A'] + dinuclCnt[x]['C'] + dinuclCnt[x]['G'] + dinuclCnt[x]['U']
    numerator = dinuclCnt[x]['A']
    if z < float(numerator) / float(denom):
        dinuclCnt[x]['A'] -= 1
        return 'A'
    numerator += dinuclCnt[x]['C']
    if z < float(numerator) / float(denom):
        dinuclCnt[x]['C'] -= 1
        return 'C'
    numerator += dinuclCnt[x]['G']
    if z < float(numerator) / float(denom):
        dinuclCnt[x]['G'] -= 1
        return 'G'
    dinuclCnt[x]['U'] -= 1
    return 'U'


def connectedToLast(edgeList, nuclList, lastCh):
    D = {}
    for x in nuclList: D[x] = 0
    for edge in edgeList:
        a = edge[0];
        b = edge[1]
        if b == lastCh: D[a] = 1
    for i in range(2):
        for edge in edgeList:
            a = edge[0];
            b = edge[1]
            if D[b] == 1: D[a] = 1
    ok = 0
    for x in nuclList:
        if x != lastCh and D[x] == 0: return 0
    return 1


def eulerian(s):
    nuclCnt, dinuclCnt, List = computeCountAndLists(s)
    # compute nucleotides appearing in s
    nuclList = []
    for x in ["A", "C", "G", "U"]:
        if x in s: nuclList.append(x)
    # compute numInList[x] = number of dinucleotides beginning with x
    numInList = {}
    for x in nuclList:
        numInList[x] = 0
        for y in nuclList:
            numInList[x] += dinuclCnt[x][y]
    # create dinucleotide shuffle L
    firstCh = s[0]  # start with first letter of s
    lastCh = s[-1]
    edgeList = []
    for x in nuclList:
        if x != lastCh: edgeList.append([x, chooseEdge(x, dinuclCnt)])
    ok = connectedToLast(edgeList, nuclList, lastCh)
    return ok, edgeList, nuclList, lastCh


def shuffleEdgeList(L):
    n = len(L);
    barrier = n
    for i in range(n - 1):
        z = int(random.random() * barrier)
        tmp = L[z]
        L[z] = L[barrier - 1]
        L[barrier - 1] = tmp
        barrier -= 1
    return L


def dinuclShuffle(s):
    if len(s)==0:
        return("")
    ok = 0
    while not ok:
        ok, edgeList, nuclList, lastCh = eulerian(s)
    nuclCnt, dinuclCnt, List = computeCountAndLists(s)

    # remove last edges from each vertex list, shuffle, then add back
    # the removed edges at end of vertex lists.
    for [x, y] in edgeList: List[x].remove(y)
    for x in nuclList: shuffleEdgeList(List[x])
    for [x, y] in edgeList: List[x].append(y)

    # construct the eulerian path
    L = [s[0]];
    prevCh = s[0]
    for i in range(len(s) - 2):
        ch = List[prevCh][0]
        L.append(ch)
        del List[prevCh][0]
        prevCh = ch
    L.append(s[-1])
    t = ''.join(L)
    return t


# print(PDBlist)
def run_BP(seq, ss, modules_to_parse, dataset, left_out):
    maxs = BayesPairing.parse_sequence(seq, modules_to_parse,ss, DATASET_NAME, left_out,p=10000)

    return maxs

def shuffle_seq(seq):
    seq = list(seq)
    shuffle(seq)
    return "".join(seq)

def get_constraints_from_BN(positions,graph):
    if len(positions) > 0:
        constraints = []
        bps = []
        ncbps = []
        bp_types = []
        real_bps = []
        real_ncbps = []
        for i in graph.edges():
            #   print(i,graph.get_edge_data(*i))
            if (graph.get_edge_data(*i)['label'].upper() == "CWW") and i[0] < i[1]:
                bps.append(i)
            elif (graph.get_edge_data(*i)['label'].upper() not in ["B53","S33","S55"]) and i[0] < i[1]:
                # print(graph.get_edge_data(*i))
                ncbps.append((i, graph.get_edge_data(*i)['label'].upper()))
        #print('BASE PAIRS')
        #print(bps)
        #print(ncbps)
        nodes = []
        for i in graph.nodes():
            nodes.append(int(i))
        sortednodes = sorted(nodes)
        #print(sortednodes)
        for j in range(len(sortednodes)):
            n = sortednodes[j]
            for bp in bps:
                (a, b) = bp
                if n == a:
                    pairing_node = positions[j]
                    partner_ind = sortednodes.index(b)
                    partner_node = positions[partner_ind]
                    real_bps.append((pairing_node, partner_node))
                elif n == b:
                    pairing_node = positions[j]
                    partner_ind = sortednodes.index(a)
                    partner_node = positions[partner_ind]
                    real_bps.append((partner_node, pairing_node))
            for bp in ncbps:
                (a, b) = bp[0]
                # print(a,b)
                if n == a:
                    pairing_node = positions[j]
                    partner_ind = sortednodes.index(b)
                    partner_node = positions[partner_ind]
                    real_ncbps.append(((pairing_node, partner_node), bp[1]))
                elif n == b:
                    pairing_node = positions[j]
                    partner_ind = sortednodes.index(a)
                    partner_node = positions[partner_ind]
                    real_ncbps.append(((partner_node, pairing_node), bp[1]))
        return (set(real_bps), set(real_ncbps))
    else:
        return ([], [])
def parse_FR3D(positions, bps, ncbps,aiming_for):
    print("AIMING FOR",aiming_for)
    print("FOUND",positions)
    #PDB_name = PDB.upper() + ".nxpickled"
    #chain = get_chain_from_PDB(PDB,positions[0])
    if len(positions)==0:
        print("POSITIONS ARE NULL, ERROR")
        return 0
    max_score =len(aiming_for)
    score = 0
    for i in positions:
        if i in aiming_for:
            score = score +1

    score = score/max_score
    print("SCORE :", score)
    return score


def compare_to_FR3D(score, positions, module_graph, chain, aiming_for):
    #print("GETTING CONSTRAINTS :", positions,module_graph)
    #print(module_graph.edges(data=True))
    #exit()
    bps, ncbps = get_constraints_from_BN(positions, module_graph)
    score = parse_FR3D(positions, bps, ncbps,aiming_for)

    return score


def get_seq_ss(PDBid,ex):
    #print(PDBid)
    PDB, chain = PDBid.split("_")[0:2]
    #print(PDB)
    # print("../all_graphs_pickled/" + PDB + ".nxpickled")
    try:
        #g = pickle.load(open("../Graphs/all_graphs_pickled/" + PDB + ".nxpickled", "rb"))
        
        g = pickle.load(open("../../../../BayesPairing1/rnabayespairing/bayespairing/models/all_graphs_pickled/" + PDB + ".nxpickled", "rb"))
    except FileNotFoundError:
        print("PDB FILE NOT FOUND")
        return ("", 0,0)
    seq = ""
    nodes = []
    for node in g.nodes(data=True):
        #print(node)
        # print(node[0][0],chain)
        if node[0][0] == chain:
            nodecode = node[0][1]
            if node[1]["nt"]!= "TU":
                nodes.append((int(nodecode), node[1]["nt"]))
            else:
                nodes.append((int(nodecode), "U"))
    sortednodes = sorted(list(nodes))
    #print("FIRST NODE:",sortednodes[0])
    nuc_by_node = {}
    missing_nuc = False
    # print("NODES")
    for i in sortednodes:
        nuc_by_node[i[0]] = i[1]
    #print(sortednodes)
    try:
        for i in range(1, int(sortednodes[-1][0]) + 1):
            if i not in nuc_by_node.keys() :
                if ("A" in seq or "G" in seq or "C" in seq or "U" in seq):
                    seq = seq + "" #should be N or gap, trying not ot crash shit.
                    #seq = seq + "N"
            else:
                seq = seq + nuc_by_node[i]
        if chain in g.graph["ss"]:
            ss = g.graph['ss'][chain]
        else:
            ss = ""
        # print(seq)
        # print("MISSING_NUC",PDBid,missing_nuc)
        if "T" in seq:
            seq = seq.replace("T","U")
        #exit()
        #print(seq)
        #exit(0)
    except:
        return ("","","")
    return (seq, ss, chain)


def run_validation(module_to_test,full_seq,LEFT_OUT,left_out_seq, leave_out_seq, aiming_for, indexes):
    results = {}

    i = module_to_test
    results[i] = []
    
    print("PARSING :",full_seq, "LEFT OUT",left_out_seq)

    maxs = parse_sequences.run_BP(full_seq, "", [i], DATASET_NAME, left_out=LEFT_OUT,leave_out_sequence=leave_out_seq,left_out_sequence=left_out_seq,indexes=indexes,samplesize=20000)

    print("ALL_RESULTS :",maxs)
    k = modules_to_test.index(i)

    for el in maxs[module_to_test]:
        print("COMPARING TO FR3D")
        modseq,positions,bp_score,sse_positions,seq_score = el
        positions2 = []
        for ii in positions:
            for ji in ii:
                positions2.append(ji)
        positions = positions2
        score1 = compare_to_FR3D(el[2], positions,
                graphs[i][0],chain="",aiming_for=aiming_for) 

        candidate_results = [score1, modseq, bp_score,seq_score]
        results[i].append(candidate_results)


    return results

def convert_seq(seq, aiming_for):
    ungapped_seq = ""
    ungapped_aiming_for = []
    
    gapped_pos = []
    for ind,nuc in enumerate(seq):
        if nuc in ["-","_","."]:
            gapped_pos.append(ind)
        else:
            ungapped_seq+=nuc
    
    for pos in aiming_for:
        n_gaps_before = sum([1 if x<pos else 0 for x in gapped_pos])
        
        ungapped_aiming_for.append(pos-n_gaps_before)
    return ungapped_seq,ungapped_aiming_for

def get_random_indexes(n):
    to_pick = list(range(n))
    samp = random.sample(to_pick,100)
    return samp

def remove_at(i, s):
    return s[:i] + s[i+1:]


def correct_aiming_for(seq,aiming_for):
    print("CORRECTING AIMING_FOR",aiming_for)
    n_removed = 0
    pos_to_remove = []
    new_target = []
    for ind,pos in enumerate(aiming_for):
        if ind==0:
            new_target.append(pos)
            continue
        dist = pos-aiming_for[ind-1]
        if 1<dist<5:
            for i in range(dist-1):
                pos_to_remove.append(aiming_for[ind-1]+i+1)
                n_removed+=1
        new_target.append(pos-n_removed)
    to_cut = sorted(pos_to_remove,reverse=True)    
    for cut in to_cut:
        seq = remove_at(cut,seq)
    print("CORRECTED AIMING_FOR",new_target)
    return seq,new_target



def cross_validation(modules_to_test,it):
    
    crossval_results = {}
    CV_SEQUENCES = []
    for i in modules_to_test:
        done_seqs = []
        TESTED_SEQUENCES = []
        crossval_results[i] = []
        MODULE_SEQUENCES = mod_sequences[i]
        THESE_SEQUENCES = test_sequences[i]
        print("NUMBER OF MODULE SEQS",len(MODULE_SEQUENCES),"NUMBER OF TEST SEQS",len(THESE_SEQUENCES))
        print("CURRENTLY WORKING ON MODULE",i,"THERE ARE",len(MODULE_SEQUENCES),"OCCURRENCES")
        
        print("MODULE SEQS")
        for zz in MODULE_SEQUENCES:
            print(zz)
            

        indexes = all_indexes[i]  
            
        print("INDEXES",len(indexes),indexes)
        for ind in indexes:
            entry = THESE_SEQUENCES[ind]
            if "-" in MODULE_SEQUENCES[ind]:
                print("GAP IN MODULE SEQUENCE")
                continue
            fseq,faiming_for = entry
            seq,aiming_for = convert_seq(fseq,faiming_for)
            print("DOING",seq,aiming_for,MODULE_SEQUENCES[ind],flush=True)

            pdb_len = len(seq)

            if pdb_len in range(10, 200) and "T" not in seq and "-" not in seq and "N" not in seq:

                new_seq = ""
                seq = list(seq)
                for s in seq:
                    s = s.upper()
                    if s not in ["A","C","G","U"] :
                        s = "A"
                    new_seq+=s
                seq = new_seq
                
                #if seq in done_seqs:
                #    print("THIS SEQ WAS TESTED BEFORE")
                #    continue
                #done_seqs.append(seq)
                seq,aiming_for=correct_aiming_for(seq,aiming_for)
                scores = run_validation(i,seq.upper(),ind, leave_out_seq=False,left_out_seq=MODULE_SEQUENCES[ind],aiming_for=aiming_for,indexes=indexes)

                print("SCORES",scores)    
                TESTED_SEQUENCES.append(seq)

                for k in scores:
                    scores[k].sort(key=lambda k: (k[2], k[0]), reverse=True)
                    print(scores[k])
                    crossval_results[i].append(scores[k])
                #====================================================
            elif pdb_len > 200:
                first = aiming_for[0]
                last = aiming_for[-1]
                new_seq, new_ss = [], []
                offset = max(0,first-100)
                if last - first < 200:
                    new_seq = seq[max(0,first-100):last+100]
                    
                    print("NEW SEQ LENGTH")
                    print(len(new_seq),[x - offset for x in aiming_for])

                    if len(new_seq)>300:
                        continue

                                    
                    seq = list(new_seq)
                    clean_seq = ""
                    for s in seq:
                        s = s.upper()
                        if s not in ["A","C","G","U"] :
                            s = "A"
                        clean_seq+=s
                    new_seq = clean_seq

                    aiming_for = [x-offset for x in aiming_for]
                    #print("BEFORE ADJUSTMENT;",aiming_for,len(new_seq))
                    new_seq,aiming_for=correct_aiming_for(new_seq,aiming_for)
                    
                    scores = run_validation(i,new_seq,ind,leave_out_seq=False,left_out_seq=MODULE_SEQUENCES[ind], aiming_for=aiming_for,indexes=indexes)

                    #print("SCORES",scores)
                    for k in scores:
                        scores[k].sort(key=lambda k: (k[2], k[0]), reverse=True)
                        crossval_results[i].append(scores[k])
                    TESTED_SEQUENCES.append(new_seq)

        cv_fn = "../output/2fold_"+it + str(i)
        #pickle.dump(crossval_results[i], open(cv_fn, "wb"))
        CV_SEQUENCES.append(TESTED_SEQUENCES)
       
    for i in crossval_results:
        for j in crossval_results[i]:
            print(j)
   #     #print(crossval_results)
    pickle.dump(crossval_results, open(str("2fold"+it+").pickle"), "wb"))
    print('FINAL RESULTS')
    for i in crossval_results:
        print("RESULTS FOR MODULE :", i)
        for j in crossval_results[i]:
            print(sorted(j, reverse=True))

    #pickle.dump(CV_SEQUENCES, open("ALL_SEQUENCE_TESTED_BP_OCT2.pickle", "wb"))
    
def test_fasta(input, modules_to_parse):
    prediction_scores = {}
    with open(input, "rU") as f:
        r_count = 0
        for record in SeqIO.parse(f, "fasta"):
            r_count = r_count + 1
            if r_count % 25 == 0:
                pickle.dump(prediction_scores, open('partial_predictions2.cPickle', 'wb'))
            id = record.id
            prediction_scores[id] = {}
            seq = str(record.seq).replace("T", "U")
            if len(seq) > 300:
                continue
            maxs = run_BP(seq, "", modules_to_parse, "NONE")
            print(maxs)
            for ind, module in enumerate(maxs):
                if len(maxs[ind]) > 0:
                    prediction_scores[id][modules_to_parse[ind]] = (maxs[ind][0], maxs[ind][1])
                else:
                    prediction_scores[id][modules_to_parse[ind]] = (0, [])
    pickle.dump(prediction_scores, open("prediction_score_3.cPickle", "wb"))


if __name__ == "__main__":
     
    modules_to_test = list(range(len(graphs)))
    cross_validation(modules_to_test,it)


