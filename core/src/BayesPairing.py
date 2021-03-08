__author__ = 'Roman'
#import bayes_to_seqpy as bs
import os
import sys
#sys.path.append('/home/mcb/rsarra2/anaconda3/envs/py35/lib/python3.6/site-packages')
import pickle
#import subprocess
from . import testSS
import heapq
from operator import itemgetter
import ast
import math
#import time
from . import make_BN_from_carnaval as makeBN
from .pgmpy.models import RNAModule
import collections
from .scanning.src.classes import SSETree, SSE
from .scanning.scanning import exact_matching, yield_matching
from anytree.iterators import PreOrderIter
import numpy as np
import itertools
import math
from Bio import AlignIO
from .folding import Fold
import RNA
import json
import os

CURRENT_DIRECTORY = os.path.dirname(__file__)

Lambda  = 0.35

md = RNA.md()
md.uniq_ML = 1
def boltzmann(e):
    """Return Boltzmann factor of given energy"""
    kT = (RNA.K0+md.temperature)*RNA.GASCONST*1e-3

    return np.exp(-e/kT)


def setup_BN(modules, BNs, dataset,left_out,leave_out_sequence=False,left_out_sequence="",Lambda=0.35,Theta=1, verbose=False,indexes=[]):
    # Load dataset here for optimmization
    with open(os.path.join(CURRENT_DIRECTORY, "../models/" + dataset + ".json")) as f:
        modules_data = json.load(f)

    for mod in modules:
        if mod not in BNs:
            print("making Bayes Net for module",mod)
            # Pass in the module_data
            BNs[mod] = makeBN.call_makeBN(mod, dataset, modules_data[str(mod)], Lambda, Theta)
    if verbose:
        print("Bayes Net dataset:", BNs.keys())
    return BNs


def get_rotations(strands):
    if len(strands)<2:
        permutations = [strands]
        return permutations
    strands = list(strands)
    current_permutations = []
    is_first = strands[0]
    current_p = strands.copy()
    current_permutations.append(tuple(current_p.copy()))
    for i in range(len(strands)-1):
        moved = current_p[0]
        current_p.append(moved)
        current_p.pop(0)
        current_permutations.append(tuple(current_p.copy()))
    return current_permutations

#by sebastian will
def parseRNAStructure(structure, *, opening = "([{<", closing = ")]}>"):
    stack = { op:list() for op in opening }
    bps = [-1]*len(structure)

    for i,c in enumerate(structure):
        for (op,cl) in zip(opening,closing):
            if c==op:
                stack[op].append(i)
            elif c==cl:
                if len(stack[op]) == 0:
                    raise ParseError("Unbalanced RNA dot-bracket structure reading "+cl+".")
                j = stack[op].pop()
                bps[i] = j
                bps[j] = i

    for op in opening:
        if len(stack[op]) > 0:
            raise ParseError("Unbalanced RNA dot-bracket structure reading "+op+".")

    return bps

def find_significant_columns(aln_sequences,struct):
    align = aln_sequences
    pairs = parseRNAStructure(struct)
    #print("BPs",pairs)
    cantBeGood = []
    good_pos = []
    for pos in range(len(align[0])):
        nuc_count = (len([x[pos] for x in align])-[x[pos] for x in align].count("-") )/len([x[pos] for x in align]) 
        #print("current column",[x[pos] for x in align])
        #print(nuc_count)
        if nuc_count>0.5:
            good_pos.append(pos)
        else:
            cantBeGood.append(pos)
            if pairs[pos]>-1:
                cantBeGood.append(pairs[pos])
        
    only_good_pos = [x for x in good_pos if x not in cantBeGood]
          
    #print("good pos", good_pos)
    return only_good_pos
	
def parse_alignment2(sequences, modules, ss, dataset, BNs, t=-3, samplesize=20000, Lambda=0.35, Theta=1, Delta=None, fuzzy=False, verbose=False):
    #print("ALIGNMENT SEQUENCES",sequences)
    seqs = [x[0] for x in sequences]
    fc = Fold(seqs)
    ss_mfe, mfe, fee = fc.constraint_folding()
    nb = samplesize
    return_dict = {}
    scores = {}
    MODEL_PATH = "../models/"

    MODULES = BNs

    #adding constrained mfe of each module
    #for m in modules:
    #    constraints = testSS.get_constraints_from_BN()
    #print(sequences)
    ss = []
    #while len(ss)<100:
    #    print('rnaalifold failed, trying again')
    #    ss = testSS.call_alifold(aln_file)
    #if verbose:
    #    print("Sampling",str(nb),"secondary structures",flush=True)

    if len(ss)<1:
        #ss = testSS.call_rnasubopt(seq)
        ss = list(fc.non_redundant_sampling(nb,delta=Delta))
    else:
        ss = [(ss,1)]

    if verbose:
        ss = list(ss)

    db = {}
    for i in modules:
        j = BNs[i]
        #print(j.ID,j.n_nodes,j.nodes,j.ss,j.stackings,j.gaps)
        sse_mod = SSE(BNs[i].canonical, root=False)
        #print(sse_mod.strands_len,sse_mod.enclosed_bp,sse_mod.closing_bp)

        strands= sse_mod.strands_len
        rotations = get_rotations(strands)
        for rotated,rotation_form in enumerate(rotations):
            db.setdefault(rotation_form, []).append((BNs[i],rotated))
    if verbose:
        print("Module database",db, flush=True)
    #FOR NOW: dont score same sequence twice. later: dont score same sequence twice, but allow multiple positions! (positions are broken right now)
    modules_predicted = {}
    

    #real_pos = find_significant_columns(seqs)


    for ind,subopt_output in enumerate(ss):
        if ind%500==0 and verbose:
            print("Parsing structure",str(ind),"of",str(len(ss)))
        struct,energy = subopt_output
        sDelta = energy - mfe
        #check threshold

        

        real_ss = "".join([struct[pos] for pos in real_pos])
        real_seq = "".join([sequences[0][0][pos] for pos in real_pos])
        struct = real_ss
		
        
        #print(struct)
        #print("===================================================")
        #print(seq)
        #print(struct)

        #print(struct[20:60])
        try:
            tree = SSETree.from_bracket(struct,seq=real_seq)
        except:
            continue        
        #print(tree.seq)


        #if no match with exact matching above some threshold, start scoring some of the fuzzy matches.

        exact_matching(tree,db, fuzzy=fuzzy)
        for r in yield_matching(tree,fuzzy=fuzzy):

            #if "-" in r.seq:
             #   continue
            #print("Found match!",r.seq,r.pos,r.module.stackings, r.node.seq, r.module.ss, r.fuzzy_form)
            #print(struct)
            if not r.module.ID in modules_predicted:
                r.eval(alignment=True,sequences = [x[0] for x in sequences], aln_sequences=seqs, ungapped_positions = real_pos, fold_compound=fc)
                if r.score > t:
                    modules_predicted.setdefault(r.module.ID, []).append([r.seq, r.pos, r.score-Lambda*sDelta,r.node.positions, r.score])
            else:
                #print(r.pos,  modules_predicted[r.module.ID])

                if len([1 for x in modules_predicted[r.module.ID] if str(x[3])==str(r.node.positions)])<1:
                    if any(r.seq in output for output in modules_predicted[r.module.ID]):
                        r.score = [x[4] for x in modules_predicted[r.module.ID] if x[0] == r.seq][0]
                        if r.score > t:
                            #print(r.module.ID, r.seq, r.pos, r.score)
                            modules_predicted.setdefault(r.module.ID,[]).append([r.seq, r.pos, r.score-Lambda*sDelta,r.node.positions,r.score])
                    else:
                        r.eval(alignment=True, sequences=[x[0] for x in sequences], aln_sequences=seqs, ungapped_positions = real_pos, fold_compound=fc)
                        if r.score>t:
                            #print(r.module.ID, r.seq, r.pos, r.score)
                            modules_predicted.setdefault(r.module.ID, []).append([r.seq, r.pos, r.score-Lambda*sDelta,r.node.positions,r.score])

    sorted(modules_predicted)
    return modules_predicted


def jared(seq, module, cutoff=-111):
    #print("Executing JAR3D-style BayesPairing")
    score = module.odds(seq)
    #print(score)
    if cutoff==-111:
        return score[0]
    else:
        return score[0]>=cutoff,score[0]

def parse_alignment(sequences, modules, ss, dataset, BNs, t=-5, samplesize=20000, Lambda=0.35, Theta=1, Delta=None, fuzzy=False, verbose=False, constraints = ''):
    seqs = [x[0] for x in sequences]
    consensuses = []
    #getting the consensus of each column
    for x in range(len(seqs[0])):
        consensuses.append(consensus(seqs,x))
    
    fc = Fold(seqs)
    ss_mfe, mfe, fee = fc.constraint_folding()
    nb = samplesize
    fuzzy=False

    if len(constraints)<len(seqs[0]):
        if verbose:
            print("calling folding",flush=True)
        ss_mfe, mfe, fee = fc.constraint_folding()
    else:
        ss_mfe, mfe, fee = fc.constraint_folding(c=constraints)
    nb = samplesize

    return_dict = {}
    scores = {}
    MODEL_PATH = "../models/"

    MODULES = BNs
    if verbose:
         print("Sampling",str(nb),"secondary structures", "with maximum Delta",str(Delta),flush=True,)
    if len(ss)<1:
        ss = fc.non_redundant_sampling(nb, delta=Delta)
    else:
        ss = [(ss,None)]
        
      

    
    ss = list(ss)[:samplesize]
    
    if verbose:
        for zz in ss:
            print(zz)
    

    db = {}
    for i in modules:
        j = BNs[i]
        j.folded_positions = []
        sse_mod = SSE(BNs[i].canonical, root=False)
        strands= sse_mod.strands_len
        rotations = get_rotations(strands)
        for rotated,rotation_form in enumerate(rotations):
            db.setdefault(rotation_form, []).append((BNs[i],rotated))
    if verbose:
        print("Module database",db, flush=True)

    modules_predicted = {}
    for mod_p in modules:
        modules_predicted[mod_p]= {}
    BOLTZMANN_SUM = 0
    for ind,subopt_output in enumerate(ss):
        if ind%500==0 and verbose:
            print("Parsing structure",str(ind),"of",str(nb))
        struct,energy = subopt_output
        if energy is not None:
            prob = boltzmann(energy)
            BOLTZMANN_SUM += prob
        else:
            prob=1
            BOLTZMANN_SUM = 1

            
        real_pos = find_significant_columns(seqs,struct)
        
        #skipped_pos contains the number of columns to add to each column ID in the reduced alignment.
        skipped_pos = []
        toAdd = 0
        for i in range(0,max(real_pos)):
            if i not in real_pos:
                toAdd+=1
            else:
                skipped_pos.append(toAdd)
            
        real_ss = "".join([struct[pos] for pos in real_pos])
        real_seq = "".join([sequences[0][0][pos] for pos in real_pos])
        
        
        #print("problematic struct",struct)
        #print("accepted columns",real_pos)
        #print("bps",parseRNAStructure(struct))
        #print("skipped",skipped_pos)   

        
        struct = real_ss

        tree = SSETree.from_bracket(struct,seq=real_seq)


        exact_matching(tree,db, fuzzy=fuzzy)
        for r in yield_matching(tree,fuzzy=fuzzy):
            if not str(r.node.positions) in modules_predicted[r.module.ID]:
                r.eval(alignment=True,sequences = [x[0] for x in sequences], aln_sequences=seqs, ungapped_positions = real_pos, fold_compound=fc)
                if r.score > t:
                    modules_predicted[r.module.ID][str(r.node.positions)] = [r.seq, r.pos,  prob,r.node.positions, r.score]
                    
            else:
                modules_predicted[r.module.ID][str(r.node.positions)][2]= modules_predicted[r.module.ID][str(r.node.positions)][2] + prob

    predicted_list = {}
    for mod in modules_predicted:
        this_mod_list = []
        for matched in modules_predicted[mod]:
            match = modules_predicted[mod][matched]
            fract = match[2]/BOLTZMANN_SUM
            if fract<0:
                print("WRONG PROB SIGN")
                continue
            try:
                energy_term = Lambda * math.log(fract)
            except:
                energy_term = 0

            new_score=match[4]+energy_term
            match[2] = new_score
            
            #this is necessary to re-add the missing columns in the alignment so the output matches the input.
            #print("attempting to insert missing columns in match", match[1])
            newRes = []
            for strand in match[1]:
                strand2 = [x+skipped_pos[x] for x in strand]
                newRes.append(strand2)
            match[1] = newRes
            #print("successfully updated:",match[1])
            #print("updating seq, before:",match[0])
            newSeq = ""
            for strand in match[1]:
                for x in strand:
                    newSeq+=consensuses[x]
            match[0] = newSeq
            #print("successfully updated seq, after:",match[0])
            
            
            
            this_mod_list.append(match)
        predicted_list[mod] = (sorted(this_mod_list,key=lambda k:k[2], reverse=True))
    return predicted_list
    

def consensus(seqs, colIndex):
    lst = [seq[colIndex] for seq in seqs]
    cand = set(lst)
    if "-" in cand and len(cand)>1:
        cand.remove("-")
    return max(cand, key=lst.count)
    
def parse_sequence(seq, modules, ss, dataset, BNs, t=-10, samplesize=20000, Lambda=0.35, Theta=1, Delta=None, fuzzy=False, verbose=False, constraints = ''):
    #t=5
    #print("PARSING SEQUENCES",seq,"LAMBDA",Lambda)
    fuzzy=False
    #print("defining fold compound",flush=True)
    fc = Fold(seq)
    if len(constraints)<len(seq):
        if verbose:
            print("calling folding",flush=True)
        ss_mfe, mfe, fee = fc.constraint_folding()
    else:
        ss_mfe, mfe, fee = fc.constraint_folding(c=constraints)
    nb = samplesize

    return_dict = {}
    scores = {}
    MODEL_PATH = "../models/"

    MODULES = BNs
    #print("ss",ss)
    #print("TESTING SEQ AND SS",seq,ss)
    if verbose:
         print("Sampling",str(nb),"secondary structures", "with maximum Delta",str(Delta),flush=True,)
    if len(ss)<1:
        #ss = testSS.call_rnasubopt(seq)
        ss = fc.non_redundant_sampling(nb, delta=Delta)
    else:
        ss = [(ss,None)]
        #print("PROVIDED SS")
        #print(seq)
        #print(ss[0][0])

    #adding constrained mfe of each module
    #for m in modules:
    

    
    ss = list(ss)[:samplesize]
    #ss = sorted(ss,reverse=True, key = lambda x:x[1])
    
    #if verbose:
    #    for zz in ss:
    #        print(zz)
    

    db = {}
    for i in modules:
        j = BNs[i]
        j.folded_positions = []
        #print(j.ID,j.n_nodes,j.nodes,j.ss,j.stackings,j.gaps,"canonical",j.canonical)
        sse_mod = SSE(BNs[i].canonical, root=False)
        #print(sse_mod.strands_len,sse_mod.enclosed_bp,sse_mod.closing_bp)

        strands= sse_mod.strands_len
        #print("SSE STRAND LEN",t)
        rotations = get_rotations(strands)
        for rotated,rotation_form in enumerate(rotations):
            db.setdefault(rotation_form, []).append((BNs[i],rotated))
    if verbose:
        print("Module database",db, flush=True)

    #FOR NOW: dont score same sequence twice. later: dont score same sequence twice, but allow multiple positions! (positions are broken right now)
    modules_predicted = {}
    for mod_p in modules:
        modules_predicted[mod_p]= {}
    BOLTZMANN_SUM = 0
    for ind,subopt_output in enumerate(ss):
        #print("PARSING STRUCTURE NUMBER",ind)
        if ind%500==0 and verbose:
            print("Parsing structure",str(ind),"of",str(nb))
        struct,energy = subopt_output
        #sDelta = energy - mfe
        #check threshold
        if energy is not None:
            prob = boltzmann(energy)
            BOLTZMANN_SUM += prob
        else:
            prob=1
            BOLTZMANN_SUM = 1

        #print("===================================================")
        #print(seq)
        #print(struct)

        #print(struct[20:60])
        tree = SSETree.from_bracket(struct,seq=seq)
        #print(tree.seq)


        exact_matching(tree,db, fuzzy=fuzzy)
        for r in yield_matching(tree,fuzzy=fuzzy):
            #print("Found match!",r.seq,r.pos,r.module.stackings, r.node.seq, r.module.ss, r.fuzzy_form, "SSE POSITIONS",r.node.positions, r.node.seq)
            #print(seq)
            #print(struct)
            if not str(r.node.positions) in modules_predicted[r.module.ID]:
                r.eval(fold_compound=fc)
                #print("SCORED MODULE",r.module.ID, r.seq, r.pos, r.score)
                #print("PASSING THRESHOLD",r.score,t)
                if r.score > t:
                    #print(r.module.ID, r.seq, r.pos, r.score)
                    #print("Adding module match to results")
                    modules_predicted[r.module.ID][str(r.node.positions)] = [r.seq, r.pos,  prob,r.node.positions, r.score]
                    
            else:
                #print("ADDING PREDICTION TO OTHERS",r.pos,  modules_predicted[r.module.ID])
                #TODO : this doesnt really work with fuzzy matches, so we end up returning the same positions many times.
                #if len([1 for x in modules_predicted[r.module.ID] if str(x[3])==str(r.node.positions)])<1:
                #print("ADDING ENERGY",energy,"TO SUM", modules_predicted[r.module.ID][str(r.node.positions)][2])
                modules_predicted[r.module.ID][str(r.node.positions)][2]= modules_predicted[r.module.ID][str(r.node.positions)][2] + prob
                #print("NEW ENERGY TOTAL FOR MOD",modules_predicted[r.module.ID][str(r.node.positions)][2])
                
                #    elif any(r.seq in output for output in modules_predicted[r.module.ID]):
                #        r.score = [x[4] for x in modules_predicted[r.module.ID] if x[0] == r.seq][0]
                #        #print("PASSING THRESHOLD-SAME SEQ",r.score,t)
                #        if r.score > t:
                #            #print(r.module.ID, r.seq, r.pos, r.score)
                #            modules_predicted.setdefault(r.module.ID,[]).append([r.seq, r.pos, energy,r.node.positions,r.score])
                #   else:
                #       r.eval(fold_compound=fc)
                #        #print("SCORED MODULE", r.module.ID, r.seq, r.pos, r.score)
                #        #print("PASSING THRESHOLD-NEW SEQ",r.score,t)
                #        if r.score>t:
                #            #print(r.module.ID, r.seq, r.pos, r.score)
                #            modules_predicted.setdefault(r.module.ID,[]).append([r.seq, r.pos, energy ,r.node.positions, r.score])
    
    predicted_list = {}
    #print("UPDATING SCORES",modules_predicted)
    for mod in modules_predicted:
        this_mod_list = []
        for matched in modules_predicted[mod]:
            match = modules_predicted[mod][matched]
            fract = match[2]/BOLTZMANN_SUM
            if fract<0:
                print("WRONG PROB SIGN")
                continue
            #print("TO LOG",match[2],ENERGY_SUM,fract)
            try:
                energy_term = Lambda * math.log(fract)
            except:
                energy_term = 0
            #print("ENERGY_TERM",energy_term,"LAMBDA",Lambda,"ENERGY",match[2],"ENERGY SUM",ENERGY_SUM)
            new_score=match[4]+energy_term
            match[2] = new_score
            this_mod_list.append(match)
        predicted_list[mod] = (sorted(this_mod_list,key=lambda k:k[2], reverse=True))
    #sorted(modules_predicted)
    
    
    return predicted_list



if __name__ == "__main__":
    args = sys.argv
    sequence = sys.argv[1]
    dataset = sys.argv[2]
    module = int(sys.argv[3])
    BN = makeBN.call_makeBN(module, dataset,"NONE",False,"", 10, 0.2)
    score = jared(sequence,BN)
    print(score)
    
