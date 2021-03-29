#!/usr/bin/env python3

import itertools
from collections import defaultdict
import logging
from operator import mul

import networkx as nx
import numpy as np
import pandas as pd

from ..base import DirectedGraph
from ..factors.discrete import TabularCPD, JointProbabilityDistribution, DiscreteFactor
from ..independencies import Independencies
from ..extern import six
from ..extern.six.moves import range, reduce
from .MarkovModel import MarkovModel
from . import BayesianModel
import matplotlib
import matplotlib.patches as mpatches
import networkx as nx
from matplotlib import pyplot as plt
from collections import namedtuple
import math

# TO RUN TESTS, USE THIS IMPORT INSTEAD
from ... import testSS


from ...folding import Fold

#Lambda = 0.2
BasePair = namedtuple('BasePair', ('fst','snd'))

#add stacking flag

#self.ss
#implementer un pipeline pour comparer la vitesse



#le plus simple en cas de scanning
#self.ss = structure secondaire sans le stack
#self.stackings = liste ( number (indexed to 0 ) de 0 au nombre de brins avec le nombre de stackings a chaque brin ) 

#have both node numbering ways (start at 0 and true)


class RNAModule(BayesianModel):


    def __init__(self):
        super(RNAModule, self).__init__()

    def from_alignment_dataset(self,dataset, id, graphs, alignment, sequences, PDBs, Lambda=0.2, Theta=1):

        #print("defining parameters")
        self.Lambda=Lambda
        self.Theta=Theta
        self.dataset=dataset
        self.ID = id
        self.graph_list = graphs
        #print("SELF GRAPHLIST EDGES", self.graph_list[0].edges(data=True))
        self.alignment = alignment
        #print("STORING TEST SEQUENCES")
        self.test_sequences = [x[0] for x in sequences]
        self.test_positions = [x[1] for x in sequences]
        #print("DONE, now storing PDBs")
        self.PDBs = PDBs
        self.positions = sorted(list(self.nodes()))
        #print("POSITIONS",self.positions)
        
        self.folded_positions = []
        #print("done defining parameters")
        #gaps are the positions in the ss(including the * and &) at which the nucleotide must not be included in the scoring
        self.gaps = []

        #gaps_per_strand are the same, but strand by strand, not including non-nucleotide symbols
        self.gaps_per_strand = []
        
        
       # self.generate_module_attributes()
        #print("Generating module graph",flush=True)
        self.generate_module_graph()
        #print("SELF GRAPH EDGES", self.graph.edges(data=True))
        #print("ADDED GAPS:",self.gaps,self.gaps_per_strand)
        #print("setting canonical base pairs and stacking info")
        self.set_can_bp()
        #print("SETTING STACKING INFO") 
        self.set_stacking_info()
        #print("setting secondary structure")
        self.set_ss()
        #self.print_module_info()
        #print("adjusting base pairs to stackings")
        #print("STACKINGS FOUND:",self.stackings,self.canonical, "SIZE OF MODULE",self.n_nodes)
        self.adjust_bps_to_stackings()
        #print("NEW CANONICAL",self.canonical)
        #print("done")
        #print("OBTAINED SS",self.ss, self.positions)


    def from_only_graphs(RNAModule,graphs, alignement):
        print("not implemented yet")


    def generate_module_attributes(self):
        self.generate_module_graph()
        self.set_can_bp()
        self.set_ss()
        self.set_stacking_info()

    def generate_module_graph(self):
        self.graph = reindex_graph_at_0(self.graph_list[0])
        self.n_nodes = len(self.graph.nodes())
        self.n_edges = len(self.graph.edges())
        self.positions = list(self.graph.nodes())

    def print_module_info(self):
        print("===== MODULE INFO =====")
        print("ID:",self.dataset,self.ID)
        print("Info:",str(len(self.nodes)),"nodes",str(len(self.edges)),"edges")
        print("Structure",self.ss,self.stackings)
        print("=======================")

    def set_can_bp(self):
        edges = sorted(self.graph.edges(data=True))
        canon = [(x[0],x[1]) for x in edges if x[2]["label"].upper()=="CWW"]
        self.canonical = canon
        return canon

    def set_ss(self):
        self.gaps_per_strand = [[] for x in self.stackings]
        bps = self.canonical
        #print("CANONICAL BP", self.canonical)
        ss = ""
        current_strand = 0
        current_strand_position = 0
        for ind,position in enumerate(self.positions):
            if ind>0:
                if position!=(self.positions[ind-1]+1):
                    diff = abs(position-self.positions[ind-1])
                    if diff>4 and list(ss)[-1]=="(":
                        ss+="*"
                        current_strand+=1
                        current_strand_position=0
                    else:
                        for zz in range(diff-1):
                            self.gaps.append(len(ss))
                            self.gaps_per_strand[current_strand].append(current_strand_position)
                            current_strand_position += 1
                            ss+="."
                        #ss+="."
            if any(position == bp[0] for bp in bps):
                ss += "("
            elif any(position == bp[1] for bp in bps):
                ss += ")"
            else:
                ss += "."
            current_strand_position+=1
        #print("SECONDARY STRUCTURE",ss,"GAPS",self.gaps,self.gaps_per_strand)
        self.ss = ss

    def set_stacking_info(self):
        stackings = []
        current_stack = 0
        for ind,bp in enumerate(self.canonical):
            if ind==len(self.canonical)-1:
                stackings.append(current_stack)
                break

            next =self.canonical[ind+1]
            if bp[0]==next[0]-1 and bp[1]==next[1]+1:
                #need to test case (((*)))
                if ind == len(self.canonical)-2:
                    current_stack+=1
                else:
                    next_next = self.canonical[ind+2]
                    if next_next[0]-1==next[0] and next_next[1]+1 == next[1]: #case where theres multiple stackings
                        stackings.append(current_stack)
                        current_stack = 0
                    else:
                        current_stack+=1
            else:
                stackings.append(current_stack)
                current_stack=0
        self.stackings = stackings


    def adjust_bps_to_stackings(self):

        bp_list = self.canonical
        self.adjacent = []
        if len(bp_list)>1:
            non_stacking_bps = []
            bp =0
            while bp < len(bp_list):

                #if closing base pair, find last of thestacked pairs and add
                if bp == 0:
                    reached_last_stacking = False
                    counter = 0
                    while not reached_last_stacking and counter<len(bp_list)-1:
                        bp_1 = bp_list[bp]
                        bp_2 = bp_list[bp+1]
                        if (abs(bp_1[0] - bp_2[0]) == 1) and (abs(bp_1[1] - bp_2[1]) == 1):
                            if bp_1 not in self.adjacent:    
                                self.adjacent.append(bp_1)
                        elif bp_1 not in self.adjacent:
                            non_stacking_bps.append(bp_1)
                            reached_last_stacking = True
                        bp+=1
                        counter+=1
                else:
                    #if enclosed base pair, visit pairs and add first of the stacked pairs
                    if bp== len(bp_list)-1:
                        if bp_list[bp] not in self.adjacent:
                            non_stacking_bps.append(bp_list[bp])
                        bp+=1
                        continue
						
                    bp_1 = bp_list[bp]
                    bp_2 = bp_list[bp+1]
                    if (abs(bp_1[0] - bp_2[0]) == 1) and (abs(bp_1[1] - bp_2[1]) == 1):
                        if bp_1 not in self.adjacent:    
                            non_stacking_bps.append(bp_1)
                        self.adjacent.append(bp_2)
                    elif bp_1 not in self.adjacent:
                        non_stacking_bps.append(bp_1)
                    bp+=1
            #print('OBTAINING STRANDS BASE PAIRS',"bp",bp_list,'non stacking',non_stacking_bps,'stacked',self.adjacent)
            bp_list = non_stacking_bps
        self.canonical = bp_list

    def draw(self):

        g0 = self.graph
        title = "Module:"+self.dataset+str(self.ID)
        labels = {}
        elabels = {}

        pos = nx.circular_layout(g0)
        nodes = nx.draw_networkx_nodes(g0, pos, nodelist=g0.nodes(), node_color='lightgrey', node_size=500,
                                       linewidth=2,
                                       alpha=0.99)
        nodes.set_edgecolor("black")
        # for index, i in enumerate(list(g0.nodes(data=True))):
        #    labels[i] = i["nuc"]
        for index, i in enumerate(list(g0.nodes(data=True))):
            labels[i[0]] = i[0]
        nx.draw_networkx_labels(g0, pos, labels, font_size=10)

        NC_edges = []
        C_edges = []
        backbone = []
        stacking = []
        for i in g0.edges(data=True):
            if i[0] > i[1]:
                continue
            i = i[0:2]
            label = g0.get_edge_data(*i)['label'].upper()
            if label not in ["B53", "S55", "S33"]:
                if label == "CWW":
                    elabels[i] = r"●"
                elif label == "TWW":
                    elabels[i] = r"○"
                elif label == "CSS":
                    elabels[i] = r"▲"
                elif label == "TSS":
                    elabels[i] = r"∆"
                elif label == "CHH":
                    elabels[i] = r"■"
                elif label == "THH":
                    elabels[i] = r"□"

                elif label == "CHW":
                    elabels[i] = r"■●"
                elif label == "THW":
                    elabels[i] = r"□○"
                elif label == "CWH":
                    elabels[i] = r"●■"
                elif label == "TWH":
                    elabels[i] = r"○□"
                elif label == "TSH":
                    elabels[i] = r"∆□"
                elif label == "CSH":
                    elabels[i] = r"▲■"
                elif label == "THS":
                    elabels[i] = r"□∆"
                elif label == "CHS":
                    elabels[i] = r"■▲"
                elif label == "CWS":
                    elabels[i] = r"■●"
                elif label == "TWS":
                    elabels[i] = r"○∆"
                elif label == "CSW":
                    elabels[i] = r"▲●"
                elif label == "TSW":
                    elabels[i] = r"∆○"
            # elif label in ["S35", "S53", "S55", "S33"]:
            #    elabels[i] = r"s     "
            else:
                elabels[i] = r""
            if label.upper() == "CWW":
                C_edges.append(i)
            elif label.upper() == "B53":
                backbone.append(i)
            elif label.upper() in ["S55", "S33", "S35", "S53"]:
                stacking.append(i)
            elif label.upper() not in ["S55", "S33", "S35", "S53", "B53", "CWW"]:
                NC_edges.append(i)
        # print(NC_edges)
        # print(elabels)
        nx.draw_networkx_edge_labels(g0, pos, elabels, font_size=10)
        nx.draw_networkx_edges(g0, pos, edgelist=NC_edges, edge_color='purple', width=4, arrows=False)
        nx.draw_networkx_edges(g0, pos, edgelist=C_edges, edge_color='green', width=4, arrows=False)
        nx.draw_networkx_edges(g0, pos, edgelist=backbone, edge_color='black', width=3, arrowsize=25)
        nx.draw_networkx_edges(g0, pos, edgelist=stacking, edge_color='orange', width=1, arrows=False)
        plt.title(title, fontsize=24)
        # print("done")
        plt.axis("off")
        NCP = mpatches.Patch(color="purple", label="Non-canonical")
        CP = mpatches.Patch(color="green", label="Canonical")
        BP = mpatches.Patch(color="black", label="Backbone")
        SP = mpatches.Patch(color="red", label="Stacking")
        # plt.legend(handles=[NCP,CP,BP,SP],prop={'size': 16})
        plt.show()


    def eval(self, sequence):
        #print('testing seq against nodes',sequence,self.nodes())
        candidate = self.candidate_from_sequence(sequence)
        #print('testing seq against nodes',sequence,self.nodes(),candidate)
        if candidate is None:
            return [-100, None]
        # in the bayesnet object, the nucleotide are stored as number
        module_probability = len(sequence)*math.log(4)

        # initialize the array of node probabilities for returning info
        node_probabilities = []

        # iterate over nodes to compute the probability of the sequence event at each node
        for node in self.nodes():
            ancestors_of_node = sorted(self._get_ancestors_of(node))
            ancestors_of_node.remove(node)
            # print("NODE,ANCESTORS",node,ancestors_of_node)
            this_node_evidence_dict = {}
            for anc in ancestors_of_node:
                this_node_evidence_dict[anc] = candidate[anc]
                #print("evaluating",node, candidate[node], this_node_evidence_dict)
            node_probability = math.log(compute_probability(self, node, candidate[node], this_node_evidence_dict))
            node_probabilities.append(node_probability)
            module_probability = module_probability + node_probability
        #print("EVALUATED",sequence,"WITH PROB",module_probability,node_probabilities)
        return module_probability, node_probabilities
    
    def odds(self, sequence):
        candidate = self.candidate_from_sequence(sequence)
        if candidate is None:
            return [-100, None]
        # in the bayesnet object, the nucleotide are stored as number
        module_probability = len(sequence)*math.log(4)

        # initialize the array of node probabilities for returning info
        node_probabilities = []

        # iterate over nodes to compute the probability of the sequence event at each node
        for node in self.nodes():
            ancestors_of_node = sorted(self._get_ancestors_of(node))
            ancestors_of_node.remove(node)
            # print("NODE,ANCESTORS",node,ancestors_of_node)
            this_node_evidence_dict = {}
            for anc in ancestors_of_node:
                this_node_evidence_dict[anc] = candidate[anc]
            node_probability = math.log(compute_probability(self, node, candidate[node], this_node_evidence_dict))
            node_probabilities.append(node_probability)
            module_probability = module_probability + node_probability
        return module_probability, node_probabilities

    def candidate_from_sequence(self, sequence):
        #print("SCORING CANDIDATE", self.ID, self.nodes(),sequence, "GAPS",self.gaps)
        #takes as input a sequence with stars, and outputs a dictionary of node values
        candidate = {}
        values = ['A', 'C', 'G', 'U']
        if "T" in sequence:
            sequence.replace("T","U")
        sequence = list(sequence)

        #gaps seem not to be necessary?
        ##if len(self.gaps)>0:
        #    newseq = []
        #    for nuc in range(len(sequence)):
        #        if nuc not in self.gaps:
        #            newseq.append(sequence[nuc])
        #    sequence = newseq

        #TODO: include gaps
        sequence = [i for i in sequence if i.isalpha()]
        if len(sequence)!=len(self.nodes()) or "T" in sequence:
            print("ERROR: SEQUENCE DOES NOT MATCH MODULE")
            print("seq",sequence,"mod",self.ID,self.ss,self.stackings,self.nodes())
            return None
        for ind,node in enumerate(self.nodes(data=False)):
            #print(ind,node)
            candidate[node] = values.index(sequence[ind])
        return candidate

    def mapping_from_sequence(self, sequence):
        #print("SCORING CANDIDATE", self.ID, self.nodes(),sequence, "GAPS",self.gaps)
        #takes as input a sequence with stars, and outputs a dictionary of node values
        candidate = {}
        if "T" in sequence:
            sequence.replace("T","U")

        sequence = list(sequence)
        sequence = [i for i in sequence if i.isalpha()]
  
        for ind,node in enumerate(self.nodes(data=False)):
            candidate[node] = sequence[ind]
        return candidate

    def eval_constraint_folding(self,seq_score,full_seq,positions,rotated_ss, rotated_gaps, fc):
        #print('testing candidate with constraint folding',rotated_ss, rotated_ss.split("*"))


        #this is now handled at the position level
        """
        rotated_gapped_ss = ""
        if "*" not in rotated_ss:
            for ind,element in enumerate(rotated_ss):
                if ind not in rotated_gaps[0]:
                    rotated_gapped_ss += element
        else:
            for strand_ind,strand in enumerate(rotated_ss.split("*")):
                print("current strand",strand)
                for ind, element in enumerate(strand):
                    if ind not in rotated_gaps[strand_ind]:
                        rotated_gapped_ss += element
                if strand_ind<len( rotated_ss.split("*"))-1:
                    rotated_gapped_ss+="*"
        """

        self.get_constraints_from_module(full_seq,positions,rotated_ss)

        _,nocons_score,_ = fc.constraint_folding()
        _,cons_score,_ = fc.constraint_folding(c=self.constrained_ss)
        Delta = cons_score - nocons_score

        current_score = seq_score - Delta*self.Lambda

        return current_score


    def eval_aln_constraint_folding(self,seq_score,positions, rotated_ss, rotated_gaps, one_seq, fc):

        """
        rotated_gapped_ss = ""
        if "*" not in rotated_ss:
            for ind,element in enumerate(rotated_ss):
                if ind not in rotated_gaps[0]:
                    rotated_gapped_ss += element
        for strand_ind,strand in enumerate( rotated_ss.split("*")):
            for ind, element in enumerate(strand):
                if ind not in rotated_gaps[strand_ind]:
                    rotated_gapped_ss += element
            if strand_ind<len( rotated_ss.split("*"))-1:
                rotated_gapped_ss+="*"
        """


        #do the same thing as the sequence but with an alignment file, instead of rnafold, rnacofold -p and parse the same info
        #print('testing candidate with constraint folding on file',aln_file)
        self.get_constraints_from_module(one_seq,positions,rotated_ss)

        _,nocons_score,_ = fc.constraint_folding()
        _,cons_score,_ = fc.constraint_folding(c=self.constrained_ss)
        Delta = cons_score - nocons_score

        current_score = seq_score - Delta*self.Lambda

        return current_score


    def get_constraints_from_module(self,seq,positions,ss=""):
        if len(ss)<1:
            ss = list(self.ss)
        else:
            ss = list(ss)
        while "*" in ss:
            ss.remove("*")


        positions = [x for y in positions for x in y]

        constrained_ss = ["."]*len(seq)
        #print("GETTING CONSTRAINTS",self.ss,ss,"LEN OF SEQ",len(seq),positions,"MODULE POSITIONS",self.positions,self.gaps,self.gaps_per_strand)
        #print("LEN SS",len(ss),"LEN POSITIONS",len(positions))
        for ind, x in enumerate(positions):
            if ss[ind] in [".","_"]:
                    constrained_ss[x] = "x"
            else:
                constrained_ss[x] = ss[ind]
        self.constrained_ss = "".join(constrained_ss)
        #print("CONSTRAINED_SS",self.constrained_ss)
        return self.constrained_ss

def compute_probability(motif, node, value, evidence_dict):

    #print('computing probability of node',node,value,'with evidence',evidence_dict)
    #print("MOTIF CPDS",motif.cpd_order,"MOTIF NODES",motif.nodes())
    #nodes = sorted(list(motif.nodes()), key=float)
    index_of_node = list(motif.cpd_order).index(node)
    cpd_of_node = motif.get_cpds()[index_of_node]
    #print("NODE",node,"ALL NODES", motif.cpd_order,len(get_immediate_parents(motif,node)),"PARENTS","CPD",cpd_of_node.values)
    #print("NODES",motif.nodes(data=True))
    #print("CPD ORDER",motif.cpd_order)
    #print("CPDS NODE", list(motif.nodes())[0], motif.get_cpds(list(motif.nodes())[0]))
    #print("CPDS NODE", list(motif.nodes())[1], motif.get_cpds(list(motif.nodes())[1]))    

    parents = get_immediate_parents(motif,node)
    #print("NODE",node,"NODE PARENT",parents)
    number_of_parents = len(parents)
    #print("NUMBER OF PARENTS", number_of_parents)

    if number_of_parents==0:
        #print("CPD OF NODE",cpd_of_node)
        node_probability = cpd_of_node.values.item(value)
        return node_probability

    elif number_of_parents==1:

        node_probability = cpd_of_node.values.item((value,evidence_dict[parents[0]]))

    elif number_of_parents>1:
        cpd_call = [value]
        #print("CPD CALLL",cpd_call)
        for parent in sorted(parents):
            cpd_call.append(evidence_dict[parent])
        #print("CPD CALLL",cpd_call)

        final_cpd_call = tuple(cpd_call)

        #print("SHAPE OF CPD", cpd_of_node.values)
        node_probability = cpd_of_node.values.item(final_cpd_call)

    return node_probability




def get_immediate_parents(g, node):
    """
    Returns the immediate parents of a node
    :param g: bayes net object
    :param node:  node of a interest
    :return: list
    """

    parents = []
    for i in g.edges():
        if i[1] == node:
            parents.append(i[0])
    return parents






def reindex_graph_at_0(g):
    new_g = nx.DiGraph()
    old_nodes = list(g.nodes(data=True))
    old_edges = list(g.edges(data=True))
    #print("NODES OF SELF GRAPH LIST",old_nodes)
    # print(old_edges)
    # print("===========================")
    new_nodes, new_edges = [], []
    fnode = min(list(g.nodes()))
    for node in old_nodes:
        new_nodes.append((node[0] - fnode, node[1]))
    for edge in old_edges:
        new_edges.append((edge[0] - fnode, edge[1] - fnode, edge[2]))
    # print(new_nodes, new_edges)
    new_g.add_nodes_from(new_nodes)
    new_g.add_edges_from(new_edges)
    return new_g
