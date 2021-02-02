import pickle
import random

import BayesPairing
import make_BN_from_carnaval
import parse_sequences
import collections
import argparse
import testSS


def computeCountAndLists(s):
    """
    Constructs nucleotide dictionary, dinucleotide dictionary and list of dicnucleotides

    :param s: Nucleotide sequence String consisting of A,C,G,U or .
    :return: Tuple of nucleotide counts, dinucleotide counts and dict of dinucleotides in List form
    """
    # WARNING: Use of function count(s,'UU') returns 1 on word UUU
    # since it apparently counts only nonoverlapping words UU
    # For this reason, we work with the indices.

    # Initialize lists and mono- and dinucleotide dictionaries
    # List is a dictionary of lists
    List = {'A': [], 'C': [], 'G': [], 'U': [], '.': []}
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
    """
    Chooses an adjacent nucleotide from the given dinucleotide count Dictionary and subtracts 1 from the count

    :param x: Nucleotide A,C,G,U
    :param dinuclCnt: Dict of dinucleotide counts
    :return:
    """

    z = random.random()
    #print(z)
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
    for x in nuclList:
        D[x] = 0
    for edge in edgeList:
        a = edge[0];
        b = edge[1]
        if b == lastCh:
            D[a] = 1
    for i in range(2):
        for edge in edgeList:
            a = edge[0];
            b = edge[1]
            if D[b] == 1: D[a] = 1
    for x in nuclList:
        if x != lastCh and D[x] == 0: return 0
    return 1


def eulerian(s):
    """

    :param s:
    :return:
    """
    nuclCnt, dinuclCnt, List = computeCountAndLists(s)
    # compute nucleotides appearing in s
    nuclList = []
    for x in ["A", "C", "G", "U"]:
        if x in s:
            nuclList.append(x)
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
        if x != lastCh:
            edgeList.append([x, chooseEdge(x, dinuclCnt)])
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

def shuffle_seq(seq):
    try:
        shuffled = dinuclShuffle(seq)
        return shuffled
    except:
        seq = list(seq)
        random.shuffle(seq)
    return "".join(seq)

# Functions above to shuffle - TODO: Look over and understand code

# Helper functions

def convert_seq(seq, aiming_for):
    """
    Function to convert sequence to ungapped sequence and corresponding positions for module
    :param seq: Target nucleotide sequence
    :param aiming_for: Module positions in sequence
    :return:
    """

    ungapped_seq = ""
    ungapped_aiming_for = []

    gapped_pos = []
    for ind, nuc in enumerate(seq):
        if nuc in ["-", "_", "."]:
            gapped_pos.append(ind)
        else:
            ungapped_seq += nuc

    for pos in aiming_for:
        n_gaps_before = sum([1 if x < pos else 0 for x in gapped_pos])

        ungapped_aiming_for.append(pos - n_gaps_before)
    return ungapped_seq, ungapped_aiming_for

# Functions to compute scores

def get_constraints_from_BN(positions, graph):
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
            elif (graph.get_edge_data(*i)['label'].upper() not in ["B53", "S33", "S55"]) and i[0] < i[1]:
                # print(graph.get_edge_data(*i))
                ncbps.append((i, graph.get_edge_data(*i)['label'].upper()))
        # print('BASE PAIRS')
        # print(bps)
        # print(ncbps)
        nodes = []
        for i in graph.nodes():
            nodes.append(int(i))
        sortednodes = sorted(nodes)
        # print(sortednodes)
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


def parse_FR3D(positions, bps, ncbps, aiming_for):
    print("AIMING FOR", aiming_for)
    print("FOUND", positions)
    # PDB_name = PDB.upper() + ".nxpickled"
    # chain = get_chain_from_PDB(PDB,positions[0])
    if len(positions) == 0:
        print("POSITIONS ARE NULL, ERROR")
        return 0
    max_score = len(aiming_for)
    score = 0
    for i in positions:
        if i in aiming_for:
            score = score + 1

    score = score / max_score
    print("SCORE :", score)
    return score


def compare_to_FR3D(score, positions, module_graph, chain, aiming_for):
    # print("GETTING CONSTRAINTS :", positions,module_graph)
    # print(module_graph.edges(data=True))
    # exit()
    bps, ncbps = get_constraints_from_BN(positions, module_graph)
    score = parse_FR3D(positions, bps, ncbps, aiming_for)

    return score

# Functions to run BP2

def train_BN(module_number, motif_sequences):
    """
    Call make_BN to create a BN from the module number and given motif sequences
    :param module_number: Module number to train
    :param motif_sequences: Module sequences to train on
    :return: Dict of BNs for each module: Only 1 module in this case
    """

    nets = collections.OrderedDict()
    print("making Bayes Net for module", module_number)
    nets[module_number] = make_BN_from_carnaval.make_BN(module=module_number, dataset=DATASET_NAME, graphs=graphs, motif_sequences=motif_sequences)
    return nets


def run_BP(module_number, sequence, BNs, ss):
    """
    Parse the input sequence and return the scores, based on parse_sequences.py
    :param module_number: Module number to test
    :param sequence: Test sequence
    :param BNs: Dict of Bayesian Networks
    :param ss: Secondary sequence of input sequence
    :return: List of module positions and scores
    """
    return_dict = BayesPairing.parse_sequence(seq=sequence, modules=[module_number], ss=ss, dataset=DATASET_NAME, BNs=BNs, t=-5, samplesize=20000, Lambda=0.35, Theta=1,
                                              Delta=None, fuzzy=False, verbose=False, constraints="")

    return parse_sequences.process_siblings(return_dict, siblings)

# Validation of BP2

# Validate on single sequence in module
def validate_sequence(module_number, sequence, BNs, aiming_for, ss):
    """
    Run BP2 with the given sequence and trained BNs and evaluate the scores of each result.

    :param module_number: Module number to test
    :param sequence: Sequence to test
    :param BNs: Bayesian network for trained modules
    :param aiming_for: Target positions of module in the sequence
    :param ss: Secondary sequence of target sequence
    :return: List of scores of each potential module found in the sequence
    """

    results = {module_number: []}

    maxs = run_BP(module_number, sequence, BNs, ss)

    print("ALL_RESULTS :", maxs)

    # best_result = list(maxs.keys())[0]
    if module_number in siblings:
        print("SEARCHING FOR GROUP OF MODULE", module_number, "SIBLINGS", siblings[module_number], "FOUND", maxs.keys())
    else:
        print("SEARCHING FOR GROUP OF MODULE", module_number, "NO SIBLINGS", module_number, "FOUND", maxs.keys())
    for best_result in list(maxs.keys()):
        for el in maxs[best_result][:100]:
            print("COMPARING TO FR3D", "MODULE", best_result, el)
            modseq, positions, bp_score, sse_positions, seq_score = el
            positions2 = []
            for ii in positions:
                for ji in ii:
                    positions2.append(ji)
            positions = positions2
            score1 = compare_to_FR3D(el[2], positions,
                                     graphs[module_number][0], chain="", aiming_for=aiming_for)

            candidate_results = [score1, modseq, bp_score, seq_score]
            results[module_number].append(candidate_results)

    return results

def run_validation(module_number, BNs, module_sequences, target_sequences,test_indexes, shuffle, ss):
    """
    Run the validation to obtain scores for each sequence and module to test.

    :param module_number: Module number to test
    :param BNs: Bayesian Network for trained modules
    :param module_sequences: Module sequences to test
    :param target_sequences: Target sequences to test
    :param test_indexes: Indexes of test sequences, used for ss creation
    :param shuffle: Boolean to shuffle sequence
    :param ss: Boolean to use secondary sequence as extra input
    :return: List of results for each sequence tested
    """
    results = []
    # Validate each module and target sequence in test split
    for i in range(0, len(module_sequences)):
        entry = target_sequences[i]
        if "-" in module_sequences[i]:
            print("GAP IN MODULE SEQUENCE")
            continue
        fseq, faiming_for = entry
        seq, aiming_for = convert_seq(fseq, faiming_for)
        print("DOING", seq, aiming_for, module_sequences[i], flush=True)

        pdb_len = len(seq)

        # length of sequence between 10 and 200
        if pdb_len in range(10, 200) and "T" not in seq and "-" not in seq and "N" not in seq:

            new_seq = ""
            seq = list(seq)
            for s in seq:
                s = s.upper()
                if s not in ["A", "C", "G", "U"]:
                    s = "A"
                new_seq += s

        elif pdb_len > 200:

            # NEED to sort aiming_for or else we will have incorrect pos
            aiming_for = sorted(aiming_for)
            first = aiming_for[0]
            last = aiming_for[-1]


            offset = max(0, first - 80)
            if last - first < 300:
                new_seq = seq[max(0, first - 80):last + 80]

                print("NEW SEQ LENGTH")
                print(len(new_seq), [x - offset for x in aiming_for])

                if len(new_seq) > 400:
                    continue

                seq = list(new_seq)
                clean_seq = ""
                for s in seq:
                    s = s.upper()
                    if s not in ["A", "C", "G", "U"]:
                        s = "A"
                    clean_seq += s
                new_seq = clean_seq

                aiming_for = [x - offset for x in aiming_for]

            else:
                continue

        if shuffle:
            new_seq = shuffle_seq(new_seq)

        custom_ss = ""
        if ss:
            modpos = aiming_for
            constraints = testSS.get_constraints_from_BN(graphs[module_number][test_indexes[i]], modpos)
            print(constraints, modpos)
            constrained_ss = ["."] * len(new_seq)
            for pp, cons in enumerate(constraints):

                # if out of bounds, skip this sequence for now
                if len(constrained_ss) < modpos[pp]:
                    continue
                constrained_ss[modpos[pp]] = constraints[pp]
            constraints = "".join(constrained_ss)
            print("CALLING RNAFOLD WITH CONSTRAINTS", constraints)
            custom_ss = testSS.call_rnafold(new_seq, cons=constraints, just_ss=True)
            print(custom_ss)

        scores = validate_sequence(module_number=module_number, sequence=new_seq, BNs=BNs, aiming_for=aiming_for,ss=custom_ss)

        for k in scores:
            scores[k].sort(key=lambda k: (k[2], k[0]), reverse=True)
            print(scores[k])
            results.append(scores[k])

    return results

# Could change to k_fold
def two_fold_CV(modules_to_test, output_file, shuffle, ss):
    """
    Function to run 2 fold CV. For each module, split the module and full length sequences into two splits, train on one and test on the other. Then repeat by switching the splits.
    Saves the output into a Dict of Lists of scores for each module found for each sequence in a pickle file.

    :param modules_to_test: List of modules to validate
    :param output_file: Output file string
    :param shuffle: Boolean to shuffle sequences
    :param ss: Boolean to include secondary sequence
    :return:
    """


    TO_SKIP = []

    for module_number in modules_to_test:
        if module_number in TO_SKIP:  # we have done a sibling
            print("Module", module_number, "skipped.")
            continue

        # Idea - open pickle results, load previous and continue
        try:
            cv_results = pickle.load(open(output_file, "rb"))
        except IOError:
            cv_results = {}


        # If key doesn't exists, validate, otherwise skip it
        if module_number not in cv_results:


            module_results = []



            module_sequences = mod_sequences[module_number]
            target_sequences = test_sequences[module_number]

            num_sequences = len(module_sequences)
            if num_sequences < 2:
                print("Module", module_number, "does not have enough data for Cross Validation, skipped.")

            indexes = list(range(0, num_sequences))

            # Shuffle and split indexes in half
            random.shuffle(indexes)

            # TODO: Fix this so it's more logical with higher K-folds
            split = [indexes[:num_sequences//2], indexes[num_sequences//2:]]

            for j in [0,1]:

                test_indexes = split[1 - j]
                train_indexes = split[0 + 1]
                modules_train = [module_sequences[i] for i in train_indexes]
                modules_test = [module_sequences[i] for i in test_indexes]
                target_test = [target_sequences[i] for i in test_indexes]

                # Train the BN
                BNs = train_BN(module_number, motif_sequences=modules_train)

                # Results = list of scores
                module_results = module_results + run_validation(module_number, BNs=BNs, module_sequences=modules_test, target_sequences=target_test,test_indexes=test_indexes, shuffle=shuffle, ss=ss)

            # Set results in dict
            cv_results[module_number] = module_results
            # Save to output file

            pickle.dump(cv_results, open(output_file, "wb"))

        # If module is in siblings, skip the rest
        if module_number in siblings:
            for sib in siblings[module_number]:
                TO_SKIP.append(sib)




if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-shuffle", help="Di-nucleotide shuffle of sequences for false discovery validation. Default: False.", action="store_true")
    parser.add_argument("-ss", help="Validate with secondary sequences. Default: False.", action="store_true")
    parser.add_argument("-d", help="Dataset to validate. Default: 3dMotifAtlas_ALL")

    args = parser.parse_args()

    if args.d is None:
        # Default Dataset to validate
        args.d = "3dMotifAtlas_RELIABLE"
    if args.shuffle is None:
        args.shuffle = False
    if args.ss is None:
        args.ss = False

    DATASET_NAME = args.d

    # Open the Dataset files
    graphs = pickle.load(open("../models/" + DATASET_NAME + "_one_of_each_graph.cPickle", "rb"))
    pdbs = pickle.load(open("../models/" + DATASET_NAME + "_PDB_names.cPickle", "rb"))
    mod_sequences, test_sequences = pickle.load(open("../models/" + DATASET_NAME + "_sequences.pickle", "rb"))
    siblings = pickle.load(open("../models/" + DATASET_NAME + "_siblings.pickle", "rb"))

    # Test all the modules
    modules_to_test = list(range(0, len(graphs)))
    #modules_to_test = [182]

    # Output file name
    file_name = "2fold_cv_" + DATASET_NAME
    if args.shuffle:
        file_name += "_shuffled"
    if args.ss:
        file_name += "_ss"
    file_name += ".pickle"


    # Run 2 fold CV
    two_fold_CV(modules_to_test=modules_to_test, output_file=file_name, shuffle=args.shuffle, ss=args.ss)

