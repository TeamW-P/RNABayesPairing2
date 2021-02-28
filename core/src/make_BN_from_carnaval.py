import os
import sys
import pickle
from collections import OrderedDict as od
from .pgmpy.models import RNAModule
from Bio import SeqIO
from . import BN_builder
import networkx as nx
import subprocess
from matplotlib import pyplot as plt

CURRENT_DIRECTORY = os.path.dirname(__file__)

USE_RFAM_SEQS=True

def get_pdb_list(g_list):
    pdb_to_fetch = []
    for i in g_list:
        if i not in pdb_to_fetch:
            pdb_to_fetch.append(i)
    return pdb_to_fetch


def get_pdb_seqs(pdb_list):
    record_dict = SeqIO.to_dict(SeqIO.parse("fasta.txt", "fasta"))
    seqs = []
    clean_dict = {}
    for i in record_dict.keys():
        id = i[:4]
        chain = i[5]
        id = id+chain
        if id in clean_dict.keys():
            clean_dict[id] = clean_dict[id]+str(record_dict[i].seq)
        else:
            clean_dict[id]= str(record_dict[i].seq)
    for i in pdb_list:
        if '.' in i:
            name = i.split('.')[0]
            chain= i.split('.')[1]
            candidate = str(i[0])+str(i[1])
        else:
            candidate = i
        if candidate in clean_dict.keys():
            seqs.append(clean_dict[candidate])
    return seqs

def filter_extra_seqs(seqs):
    good_seqs = []
    for seq in seqs:

        if len(seq)<3:
            continue
        seq = "".join(seq)

        if "\n" in seq:
            seq = seq.replace("\n","")
        if "-" in seq:
            seq = seq.replace("-","")
        if len(seq)>3:
            good_seqs.append(seq)

    return good_seqs

def get_alignment(g,aln,extra_seqs):
    acceptable_extra_seqs = filter_extra_seqs(extra_seqs)

    #print("processed_extra seqs",acceptable_extra_seqs)
    #print(len(acceptable_extra_seqs),"SEQS TO BE ADDED",acceptable_extra_seqs,"FOR MODULE WITH NODES",aln.keys())

    aln_dict = od()
    nodes = aln.keys()
    snodes = sorted(nodes)
    #print("NODES",snodes)

    for counter,i in enumerate(snodes):
        aln_dict[i] = {}
        #for j in range(1):
        #    index = list(g[j].nodes()).index(aln[i][j])
        #    nuc = list(g[j].nodes(data=True))[index][1]['nuc']
        #    aln_dict[i][j] = nuc
        j=0
        if USE_RFAM_SEQS==True:
            for k in range(len(acceptable_extra_seqs)):
                #print(aln_dict[i])
                #print("SEQ ADDED",acceptable_extra_seqs[k],"AT POSITION",counter)
                j = j+1
                aln_dict[i][j]=acceptable_extra_seqs[k][counter]
    return(aln_dict)

##def get_addresses(aln):
#    addresses = od()
#    for i in range(len(aln.values().__iter__().__next__())):
#        addresses[i] = []
#        keys = aln.keys()
#        sortedkeys = sorted(keys)
#        for j in sortedkeys:
#            addresses[i].append(aln[j][i])
#    return(addresses)

def get_seqs (g,addresses):
    seqs = {}
    for i in range(len(g)):
        seqs[i] = ""
        for j in range(max(addresses[i])+1):
            if j in g[i].nodes():
                index = list(g[i].nodes()).index(j)
                nuc = list(g[i].nodes(data=True))[index][1]['nuc']
                seqs[i] = seqs[i] + nuc
            else:
                seqs[i] = seqs[i] + '-'
    return seqs

#currently narrowest
def get_real_size(pos):
    real_size = 0
    for ind,node in enumerate(pos):
        if ind==0:
            continue
        if 1<node-pos[ind-1]<5:
            real_size=+node-pos[ind-1]
    return real_size
    
def get_widest_graph(g,aiming_fors,indexes):
    new_nodes = {}
    selected_ind = -1
    max_size = 1000 #is actually min right now
    for ind,aiming_for in enumerate(aiming_fors):
        if ind not in indexes:
            continue
        size = get_real_size(aiming_for)
        #print(aiming_for,size)
        if size<max_size:
            max_size=size
            selected_ind = ind
            
    for indd,node in enumerate(sorted(list(g.nodes()))):
        new_nodes[node]= aiming_fors[selected_ind][indd]
    
    g2 = nx.relabel_nodes(g,new_nodes,copy=True)
    #print("WIDEST GRAPH",g.nodes(),g2.nodes(),max_size)
    return g2


def make_graph_from_carnaval(g, alignment):
    #import draw_single_graph
    #draw_single_graph.draw_graph(g[0],node_numbers=True)
    #edges_ss = []
    #edges_3d = []
    #edges_backbone = []
    #nodes3d = []
    #alignment,addresses,seqs = aln_dict
    #nodes = []
    #for i in g[0].nodes():
    #    nodes.append(int(i))
    #edges = g[0].edges()
    #correct_edges = []
    #for i in edges:
    #    if i[0]<i[1]:
    #        correct_edges.append(i)
    #edges = correct_edges
    ##print("ALL EDGES",edges)
    #for i in edges:
    #    if g[0].get_edge_data(*i)['label']=='b53' or g[0].get_edge_data(*i)['label']=='B53':
    #        edges_backbone.append((int(i[0]),int(i[1])))
    #        edges_ss.append((int(i[0]),int(i[1])))
    #    #CHANGE: added backbone edges to secondary intucture. BEFORE: edges_backbone.append((int(i[0]),int(i[1])))
    #    elif g[0].get_edge_data(*i)['label'].upper()=='CWW':
    #        edges_ss.append((int(i[0]),int(i[1])))
    #    elif g[0].get_edge_data(*i)['label'].upper() in ["S33","S35","S53","S55"]:
    #        if len(g[0].nodes())<100:
    #            edges_3d.append((int(i[0]), int(i[1])))
    #    else:
    #        edges_3d.append((int(i[0]),int(i[1])))
    #        if int(i[0]) not in nodes3d:
    #            nodes3d.append(int(i[0]))
    #        if int(i[1]) not in nodes3d:
    #            nodes3d.append(int(i[1]))#

    #print("BACKBONE EDGES",edges_backbone)
    
    #for node in nodes:
    #    node = str(node)
    #    n_node_interactions = 0
    #    needs_backbone = False
    #    for edge in edges_ss:
    #        if edge[0]==node or edge[1]==node:
    #            n_node_interactions=n_node_interactions+1
    #    for edge in edges_3d:
    ##        if edge[0]==node or edge[1]==node:
    #            n_node_interactions=n_node_interactions+1      
    #    print(n_node_interactions) 
    #    if n_node_interactions<1:
    #        needs_backbone=True#
    #
    #    if needs_backbone:
    #        print("added_backbone for node", node)
    #        for edge in edges_backbone:
    #            if edge[0]==node or edge[1]==node:
    #                print("added backbone edge")
    #                edges_ss.append((str(edge[0]),str(edge[1])))

    """
    for node in nodes:
        this_node_partners = []
        for edge in edges_ss:
            if edge[0]==node:
                this_node_partners.append(edge[1])
            if edge[1]==node:
                this_node_partners.append(edge[0])
        for edge in edges_3d:
            if edge[0]==node:
                this_node_partners.append(edge[1])
            if edge[1]==node:
                this_node_partners.append(edge[0])
        if all(int(i)>=int(node) for i in this_node_partners):
            for backbone_edge in edges_backbone:
                if backbone_edge[1]==node:
                    #print("connected backbone to rest of graph")
                    edges_ss.append((str(backbone_edge[0]),str(backbone_edge[1])))
        #if all(int(i)<=int(node) for i in this_node_partners):
        #    for backbone_edge in edges_backbone:
        #        if backbone_edge[0]==node:
        #            print("connected backbone to rest of graph")
        #            edges_ss.append((str(backbone_edge[0]),str(backbone_edge[1])))
    """


    motif = RNAModule()
    motif.add_nodes_from(g)
    motif.add_edges_from(g.edges(data=True))
    #motif.add_nodes_from(g.nodes(data=True))
    #motif.add_edges_from(g.edges(data=True))
    #motif.add_edges_from(edges_3d)

    #print("MOTIF FINAL NODES", motif.nodes())
    #print("MOTIF FINAL EDGES", motif.edges())

    #labels = {}
    #pos = nx.circular_layout(motif)
    #nodes = nx.draw_networkx_nodes(motif, pos, nodelist=motif.nodes(), node_color='lightgrey', node_size=500, linewidth=2,alpha=0.99)
    #for index, i in enumerate(list(motif.nodes(data=True))):
    #    labels[i[0]] = i[0]
    #nx.draw_networkx_labels(motif, pos, labels, font_size=10)
    #nx.draw_networkx_edges(motif, pos, edgelist=motif.edges(), edge_color='purple', width=4, arrows=False)
   # 
    #plt.show()

    return motif

def call_makeBN(mod,dataset,left_out, leave_out_seq = False, left_out_seq = "", Lambda=1, Theta=1,indexes=[],kfold=False,retrain=False):
    ok_indexes = []
    current_ID = mod
    excluded = left_out
    g_list = pickle.load(open(os.path.join(CURRENT_DIRECTORY, "../models/"+dataset + "_one_of_each_graph.cPickle"),'rb'))
    seq_list = pickle.load(open(os.path.join(CURRENT_DIRECTORY, "../models/"+dataset + "_sequences.pickle"),'rb'))
    try:
        motif_seqs  = seq_list[0][mod]
    except:
        motif_seqs = seq_list[mod]
    #print("SEQUENCES",motif_seqs)

    #motif_seqs  = pickle.load(open("../models/"+dataset + "_sequences.pickle",'rb'))[0][mod]
    #print("MOTIF SEQS", motif_seqs)

    extra_seqs=[]
    test_seqs = []
    g = g_list[current_ID][0]
    if excluded == "NONE":
        if(os.path.isfile(os.path.join(CURRENT_DIRECTORY, "../models/" + dataset + "_models.pickle" ))) and retrain==False:
            nets = pickle.load(open(os.path.join(CURRENT_DIRECTORY, "../models/" + dataset + "_models.pickle"), "rb"))
            if mod in nets:
                return nets[mod]
            else:
                aln = {}
                for n in sorted(list(g.nodes())):
                    aln[n] = n
                alignment = get_alignment(g,aln,motif_seqs)            
            
        else:
            try:
                existing_models  = pickle.load(open(os.path.join(CURRENT_DIRECTORY, "../models/"+dataset + "_models.pickle"),'rb'))
            except:
                existing_models = {}
            if current_ID in existing_models:
                return existing_models[current_ID]
            else:

                aln = {}
                for n in sorted(list(g.nodes())):
                    aln[n] = n
                alignment = get_alignment(g,aln,motif_seqs)            
            
            
            
            
    else:
        # print("Excluded is not none",excluded)
        if len(indexes)<2:
            indexes = list(range(len(g_list)))
        try:
            test_seqs  = pickle.load(open(os.path.join(CURRENT_DIRECTORY, "../models/"+dataset + "_sequences.pickle"),'rb'))[1][mod]
        except:
            test_seqs = []
        
        if kfold==True:
            excluded_indexes = []
            left_out_index = indexes.index(excluded)
            if left_out_index<=int(len(indexes)/2):
                excluded_indexes = indexes[:int(len(indexes)/2)]
            else:
                excluded_indexes = indexes[int(len(indexes)/2):]

        else:
            excluded_indexes = []
            left_out_index = indexes.index
            for ind in indexes:
                motif_seq = motif_seqs[ind]

                if excluded==ind or (leave_out_seq and motif_seq == left_out_seq):
                    excluded_indexes.append(ind)


        if 0 in excluded_indexes:
            excluded_indexes.remove(0)
        #print("EXCLUDED INDEXES",excluded_indexes)
        okay_indexes= []
        for inde in indexes:
            if inde not in excluded_indexes:
                ok_indexes.append(inde)
                extra_seqs.append(motif_seqs[inde])
        #print("DATASET",extra_seqs)
        g = g_list[current_ID][0]
        #print("AIMINGFORS",[x[1] for x in test_seqs])
        #currently changed to narrowest
        #g = get_widest_graph(g,[x[1] for x in test_seqs],ok_indexes)
        #print("SELECTED GRAPH",g.nodes())
        aln = {}
        for n in sorted(list(g.nodes())):
            aln[n] = n
        alignment = get_alignment(g,aln,extra_seqs)

    

                       

    motif = make_graph_from_carnaval(g, alignment)
    pwm = BN_builder.build_pwm(sorted(list(motif.nodes())),alignment)

    BN  = BN_builder.build_BN(motif,pwm,alignment)

    BN.from_alignment_dataset(dataset, mod, [g], alignment, test_seqs, [], Lambda, Theta)
    return BN

# Temporary placeholder function for CV for building a Bayesian Network
def make_BN(module, dataset, graphs, motif_sequences, Lambda=0.35, Theta=1):

    # Test sequences should be tuple of seq and positions, seems like it is outdated/not used
    test_seqs = []

    # Build the BN with the motif sequences
    g = graphs[module][0]
    aln = {}
    for n in sorted(list(g.nodes())):
        aln[n] = n
    alignment = get_alignment(g, aln, motif_sequences)

    motif = make_graph_from_carnaval(g, alignment)
    pwm = BN_builder.build_pwm(sorted(list(motif.nodes())), alignment)

    BN = BN_builder.build_BN(motif, pwm, alignment)

    # PDBs array is also empty, outdated?
    BN.from_alignment_dataset(dataset, module, [g], alignment, test_seqs, [], Lambda, Theta)
    return BN


if __name__ == "__main__":
    print("Please run this from parse_sequences.py")
