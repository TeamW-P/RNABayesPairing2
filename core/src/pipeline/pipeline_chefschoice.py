
import sys
import pickle
import re
import random
from itertools import combinations
import networkx as nx

import RNA
from .. import chefschoice
import varnaapi as va
import os
import tempfile

CURRENT_DIRECTORY = os.path.dirname(__file__)

#tmp = pickle.load(open('../models/3dMotifAtlas_ALL_one_of_each_graph.cPickle', 'rb'))
#MODULES = [t[0] for t in tmp]

# This class is a modification of chefschoice. The primary change is that files are not permanently saved on the disk.
# Instead, a tempfile is generated, and this is returned to the caller, who can then retrieve and destroy the file.

def draw_structure(ss, junctions, match_dist, sequence=None):
    loops = chefschoice.decomposition(ss)
    matches = []
    for loop in loops:
        res = match_dist.get(loop)
        if res is not None:
            matches.append(res)
    matches += list(junctions)
    # print(matches)
    v = va.VARNA(structure=ss, seq=sequence)
    v.set_numeric_params(resolution=10)

    # Note that varnaapi is 0-indexed
    for m in matches:
        color = chefschoice.random_hex()
        v.add_annotation(va.LoopAnnotation(
            str(m.name), anchor=m.get_bps()[0][0]-1, size=6))
        for pos in m.get_positions():
            v.add_highlight_region(pos[0]-1, pos[-1]-1, fill=color, radius=10)

    temp = tempfile.NamedTemporaryFile(suffix=".svg")
    v.savefig(temp.name)
    temp.seek(0)
    return matches, temp


def main(seq, candidateLst):
    print(seq)
    fc = RNA.fold_compound(seq)
    candidates = [chefschoice.BPCandidate(*t)
                  for t in candidateLst if len(t[2]) <= 2]
    fc = chefschoice.soft_constraint(fc, candidates)
    junction_candidates = [chefschoice.BPCandidate(
        *t) for t in candidateLst if len(t[2]) > 2]

    match_dist = chefschoice.candidates_to_matches(candidates)

    # for lst in match_dist.values():
    #     lst.sort(key=lambda t:t.score, reverse=True)

    results = []
    # No junction included
    #print('Pure soft constraint')
    ss, mfe = fc.mfe()
    results.append([mfe, ss, []])

    # Convert all valid junctions combinaison into hard constraint
    for junctions in chefschoice.candidates_generator(junction_candidates):
        fc.hc_init()
        # Add hard constraint introduced by junctions
        constraint = ['.'] * (length+1)
        chefschoice.gen_hard_constraint(constraint, junctions)
        hd_constraint = ''.join(constraint)[1:]
        # print('hard constraint:')
        # print(hd_constraint)
        # print(seq)
        fc.hc_add_from_db(
            hd_constraint, RNA.CONSTRAINT_DB_ENFORCE_BP | RNA.CONSTRAINT_DB_DEFAULT)
        ss, mfe = fc.mfe()
        # Add energy of junctions back
        mfe = chefschoice.post_mfe_adjust(mfe, junctions)
        results.append([mfe, ss, junctions])

    # Secondary structure with minimum MFE
    min_result = min(results)
    # print(min_result)
    # print("match_dist",match_dist)
    matches, temp = draw_structure(
        min_result[1], min_result[2], match_dist, sequence=seq)

    s = [str((i+1) % 10) for i in range(len(seq))]
    print(''.join(s))
    print(seq)
    print(min_result[1])

    for i in matches:
        print(i.name)
    #    print(i.score)
    #    print(i.real_pos)
        print(i.moduleInfo)
    return ([(i.name, i.moduleInfo) for i in matches], min_result[1]), temp


def bp_chefs_choice(dist, seq, min_score=2):
    candidateLst = chefschoice.parse_candidates_from_dist(dist, min_score)
    modules_in_svg, temp = main(seq, candidateLst)
    return (modules_in_svg), temp
