
import sys
import pickle
import re
import random
from itertools import combinations
import networkx as nx

import RNA

import varnaapi as va
import os

CURRENT_DIRECTORY = os.path.dirname(__file__)

#tmp = pickle.load(open('../models/3dMotifAtlas_ALL_one_of_each_graph.cPickle', 'rb'))
#MODULES = [t[0] for t in tmp]

def random_hex():
    """Light color generator"""
    rand = lambda: random.randint(100, 255)
    return '#{:02X}{:02X}{:02X}'.format(rand(), rand(), rand())


def energy_of_score(score):
    return - int(score * 100)


def strands_from_graph(G):
    #print(G.edges(data=True))
    backbones = [(i, j) for (i, j, data) in G.edges(data=True) if data['label'] == 'B53']
    H = G.edge_subgraph(backbones)
    strands = [sorted(list(s)) for s in sorted(nx.connected_components(H.to_undirected()), key=min)]
    return strands

class BPCandidate:
    def __init__(self, name, score, position, seq,moduleInfo):
        #print('-------------')
        #print(name)
        self.name = name
        self.energy = energy_of_score(score)
        self.score = score
        self._parse_position(position)
        self.sequence = seq
        self.moduleInfo = moduleInfo
        # self.import_module()

    def _parse_position(self, pos):
        self.real_pos = pos
        length = len(pos)
        lst = [(t[0],t[-1]) for t in pos]
        bps = [(lst[0][0], lst[-1][-1])]
        for i in range(len(lst) - 1):
            bps.append((lst[i][-1], lst[i + 1][0]))
        self.adjusted_pos = lst
        self.bps = tuple(bps)
        gaps = []
        for p in pos:
            for i in range(p[0], p[-1]+1):
                if i not in p:
                    gaps.append(i)
        self.gaps = gaps

        if length == 1:
            self.decomp = RNA.DECOMP_PAIR_HP
            self.labels = (lst[0][0], lst[0][-1])
        elif length == 2:
            self.decomp = RNA.DECOMP_PAIR_IL
            self.labels = (lst[0][0], lst[1][-1], lst[0][-1], lst[1][0])

    def import_module(self, MODULES):
        self.aux_bps = []
        module = MODULES[self.name]
        strands = strands_from_graph(module)
        pos_map = {}
        #print(self.get_positions())
        #print(strands)
        for i in range(len(strands)):
            start = strands[i][0]
            new_start = self.get_positions()[i][0]
            end = strands[i][-1]
            for j in range(end - start + 1):
                pos_map[start + j] = new_start + j

        non_canonical_bps = [(pos_map[i], pos_map[j], data['label']) for (i, j, data) in module.edges(data=True) if data['label'] not in ['CWW', 'B53']]
        self.non_canonical_bps = non_canonical_bps
        # if len(self.position) == 1:
        #     gap = self.position[0][0] - min(module.nodes())
        #     for i, j in module.edges():
        #         label = module.get_edge_data(i, j)['label']
        #         if label not in ['B53', 'CWW']:
        #             if label[0] == 'T':
        #                 stericity = "trans"
        #             else:
        #                 stericity = "cis"


    def get_positions(self):
        return self.adjusted_pos

    def get_bps(self):
        return self.bps

def energy_of_candidates(candidates):
    dist = {RNA.DECOMP_PAIR_HP: {}, RNA.DECOMP_PAIR_IL: {}}
    for cd in candidates:
        dist[cd.decomp][cd.labels] = cd.energy

    return dist

def module_energy(i, j, k, l, d, data):
    if d == RNA.DECOMP_PAIR_HP:
        return data[d].get((i,j), 0)
    elif d == RNA.DECOMP_PAIR_IL:
        return data[d].get((i,j,k,l), 0)
    else:
        return 0

def soft_constraint(fc, candidates):
    data = energy_of_candidates(candidates)
    fc.sc_add_f(module_energy)
    fc.sc_add_data(data, None)
    return fc

def candidates_generator(lst):
    length = len(lst)
    for i in range(1, length + 1):
        for x in combinations(lst, i):
            if non_overlap(x):
                yield x


def non_overlap(lst):
    """Return True if given junction list is consistent
    """
    bps_lst = [t.get_bps() for t in lst]
    bps_lst.sort(key=lambda t: (t[0][0], -t[0][1]))
    length = len(lst)
    res = []
    for i in range(length):
        for j in range(i + 1, length):
            root = bps_lst[j][0]
            res.append(any([bp[0] < root[0] and bp[1] > root[1] for bp in bps_lst[i][1:]]))
    return all(res)

def gen_hard_constraint(constraint, junctions):
    for junction in junctions:
        for i, j in junction.adjusted_pos:
            for ind in range(i, j+1):
                constraint[ind] = 'x'
        for i, j in junction.get_bps():
            constraint[i] = '('
            constraint[j] = ')'


def post_mfe_adjust(mfe, junctions):
    for junction in junctions:
        mfe += junction.energy/100
    return mfe


def bracket_to_index(inst):
    """
    For a given bracket-dot presented secondary structure S, the function returns an iterger list L.
    L[i] = j if i and j are paired in S.
    """
    res = [-1]*(len(inst)+2)
    tmp = []
    for i,c in enumerate('('+inst+')'):
        if c == '(':
            tmp.append(i)
        elif c == ')':
            j = tmp.pop()
            res[i], res[j] = j, i
    return res


def decomposition(inst):
    """
    Decompose a given bracket-dot presented RNA secondary structure into several basic components 
    in tree-presentation.
    A basic component is presented by a list of its paired bases positions
    """
    index = bracket_to_index(inst)
    def aux(ind):
        """
        A recursive function decomposing a given RNA secondary structure in index list 
        from a given starting position.
        """
        tmp = []
        res = []
        k = ind + 1
        while True:
            if index[k] == -1 :
                k += 1
            elif index[k] > k:
                tmp.append((k, index[k]))
                res += aux(k)
                k = index[k]+1
            else:
                break
        return [[(ind,index[ind])]+tmp]+res

    res = aux(0)
    return [tuple(t) for t in res[1:]]

def draw_structure(ss, junctions, match_dist, sequence=None,NAME="test"):
    loops = decomposition(ss)
    matches = []
    for loop in loops:
        res = match_dist.get(loop)
        if res is not None:
            matches.append(res)
    matches += list(junctions)
    #print(matches)
    v = va.VARNA(structure=ss, seq=sequence)
    v.set_numeric_params(resolution=10)

    # Note that varnaapi is 0-indexed
    for m in matches:
        color = random_hex()
        v.add_annotation(va.LoopAnnotation(str(m.name), anchor=m.get_bps()[0][0]-1, size=6))
        for pos in m.get_positions():
            v.add_highlight_region(pos[0]-1, pos[-1]-1, fill=color, radius=10)

        # for i, j, label in m.non_canonical_bps:
        #     if label[0] == 'S':
        #         v.add_aux_BP(i-1, j-1, edge5='s', edge3='s', color='#F8E854', thickness=0.5)
        #     else:
        #         edge5, edge3 = label[1].lower(), label[2].lower()
        #         if edge5 == 'w':
        #             edge5 = 'wc'
        #         if edge3 == 'w':
        #             edge3 = 'wc'
        #         if label[0] == 'T':
        #             stericity = 'trans'
        #         else:
        #             stericity = 'cis'
        #         v.add_aux_BP(i-1, j-1, stericity=stericity, edge5=edge5, edge3=edge3, color='#800080')
    v.savefig(os.path.join(CURRENT_DIRECTORY, '../output/'+NAME+'.svg'))
    return matches

def candidates_to_matches(candidates):
    match_dist = {}
    for candidate in candidates:
        bps = candidate.get_bps()
        current = match_dist.get(bps)
        if current is None:
            match_dist[bps] = candidate
        else:
            match_dist[bps] = max(current, candidate, key=lambda t: t.score)
    return match_dist


def main(seq, candidateLst, NAME):
    fc = RNA.fold_compound(seq)
    length = len(seq)
    candidates = [BPCandidate(*t) for t in candidateLst if len(t[2]) <= 2]
    fc = soft_constraint(fc, candidates)
    junction_candidates = [BPCandidate(*t) for t in candidateLst if len(t[2]) > 2]

    match_dist = candidates_to_matches(candidates)

    # for lst in match_dist.values():
    #     lst.sort(key=lambda t:t.score, reverse=True)

    results = []
    # No junction included
    #print('Pure soft constraint')
    ss, mfe = fc.mfe()
    results.append([mfe, ss, []])

    # Convert all valid junctions combinaison into hard constraint
    for junctions in candidates_generator(junction_candidates):
        fc.hc_init()
        # Add hard constraint introduced by junctions
        constraint = ['.'] * (length+1)
        gen_hard_constraint(constraint, junctions)
        hd_constraint = ''.join(constraint)[1:]
        # print('hard constraint:')
        # print(hd_constraint)
        # print(seq)
        fc.hc_add_from_db(hd_constraint, RNA.CONSTRAINT_DB_ENFORCE_BP | RNA.CONSTRAINT_DB_DEFAULT)
        ss, mfe = fc.mfe()
        # Add energy of junctions back
        mfe = post_mfe_adjust(mfe, junctions)
        results.append([mfe, ss, junctions])

    # Secondary structure with minimum MFE
    min_result = min(results)
    #print(min_result)
    #print("match_dist",match_dist)
    matches = draw_structure(min_result[1], min_result[2], match_dist, sequence=seq, NAME=NAME)

    s = [str((i+1)%10) for i in range(len(seq))]
    print(''.join(s))
    print(seq)
    print(min_result[1])

    #for i in matches:
        #print(i.name)
    #    print(i.score)
    #    print(i.real_pos)
        #print(i.moduleInfo)
    return [(i.name,i.moduleInfo) for i in matches], min_result[1]

def parse_candidates_from_dist(dist, min_score=2):
    res = []
    for name, lst in dist.items():
        for t in lst:
            if t[2] >= min_score:
                # Note that RNAlib sequence is 1-indexed
                res.append([name, t[2], [[i+1 for i in l] for l in t[3]], t[0], t])
    return res

def parse_sequences_file(fpath):
    seqPtn = re.compile('[ACGUS-]+')
    with open(fpath) as f:
        res = [seqPtn.findall(line)[-1] for line in f.readlines() if not line[0] == '#']
    #print([len(w) for w in res])
    return res


def bp_chefs_choice(dist,seq,min_score=2,NAME="test"):

    candidateLst = parse_candidates_from_dist(dist, min_score)
    modules_in_svg = main(seq, candidateLst,NAME)
    return modules_in_svg

if __name__ == '__main__':
    # s = pickle.load(open('tpp_example.pickle', 'rb'))
    # dist = s['AL935263.2/99249-99354107/99249-99354'][0]
    dist = pickle.load(open(os.path.join(CURRENT_DIRECTORY, 'RF059.pickle'), 'rb'))

    candidateLst = parse_candidates_from_dist(dist, min_score=2)


    #to call this



    print(candidateLst)

    seq = 'AAUAGUUACUGGGGGUGCCCGCUUUCGGGCUGAGAGAGAAGGCAAGCUUCUUAACCCUUUGGACCUGAUCUGGUUCGUACCAGCGUGGGGAAGUAGAGGAAUUGUUUU'
    # seq = 'AAGUUGCACCCGGGGUGCCUGUAUUCUCAACGAUCUCAAGGCCUCUUGUCCUGGAUUGUUGUGAAUUGGGCUGAGCAAGUCCCUAUGGACCUGAACAGGAUAAUGCCUGCGAAGGGAGUGUGCAUUUCUACUUUU'
    # print(len(seq))
    main(seq, candidateLst,"bob")
    #print("INITIAL MODULE RESULTS")
    #print(dist)
    # print('--------------------------')
    # print('AliFold')
    # seqsPath = 'random_test_RF00059.stockholm.txt'
    # seqs = parse_sequences_file(seqsPath)
    # main(seqs, candidateLst)
