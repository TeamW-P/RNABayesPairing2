import pickle
import networkx as nx
import matplotlib
from matplotlib import pyplot as plt
import os
from networkx.algorithms.isomorphism import DiGraphMatcher as DGM
from collections import OrderedDict as od
import networkx.algorithms.isomorphism as iso


def get_seq_ss(PDBid):
    PDB, chain = PDBid.split(".")[0:2]
    # print(PDB)
    # print("../all_graphs_pickled/" + PDB + ".nxpickled")
    try:
        g = pickle.load(open("../models/all_graphs_pickled/" + PDB + ".nxpickled", "rb"))
    except FileNotFoundError:
        #print("PDB FILE NOT FOUND")
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

def align_multiple_graphs(g_list):
    matches = {}
    etalon = g_list[0]
    for z in etalon.nodes():
        matches[z]=[]
    #g_list.remove(g_list[0])
    for i in g_list:
        corr = align_graphs(etalon,i)
        for j in corr:
            matches[j[0]].append(j[1])
    return matches


def get_edges_at_node(g1,node):
    edgez = []
    for i in g1.edges():
        if i[0] == node or i[1] == node:
            edgez.append(i)
    return edgez


def get_edge_between_two_nodes(g1,node1,node2):
    if node2[0]==node1:
        node2 = node2[1]
    else:
        node2=node2[0]
    #print('GET EDGE BETWEEN NODES',node1,node2)
    edge = [-1,-1]
    for i in g1.edges():
        #print('THIS_EDGE=',i)
        if i==(node1,node2) or i==(node2,node1):
            edge = i

   # print('RETURNED',edge)
    return edge

def edge_match(g1,g2,edge1,edge2):
    edge_data00 = g1.get_edge_data(*edge1)
    edge_data11 = g2.get_edge_data(*edge2)

    if edge_data00 == edge_data11:
        #print("EDGES MATCHED!",edge_data00,edge_data11)
       # print(edge_data00)
      #  print(edge_data11)
        return True

    #else:
      #  print("EDGES DID NOT MATCH!", edge_data00, edge_data11)

        # print(edge_data00)
      #  print(edge_data01)
      #  print(edge_data10)
      #  print(edge_data11)
        return False

def add_to_corr(g1,g2,correspondances,matches,current_expansion=-1):
    success=False
    neighbors = g1.neighbors(correspondances[current_expansion][0])
    neighbors2 = g2.neighbors(correspondances[current_expansion][1])
    for u in neighbors:
        #print('CORRESPONDANCES:',correspondances[current_expansion][0],u, neighbors)
        if u not in [x[0] for x in correspondances]:
            #print('TESTED EDGE:',u)
            for v in matches[u]:
                #print(v,matches[u])
                if v not in [x[1] for x in correspondances]:
                   # print('TRYING PARTERN:',v)
                    if v in neighbors2:
                        edge1 = get_edge_between_two_nodes(g1,correspondances[current_expansion][0],u)
                        edge2 = get_edge_between_two_nodes(g2,correspondances[current_expansion][1],v)
                        #print('EDGEMATCH BETWEEN U AND V:',edge1,edge2)
                        edgematch = edge_match(g1,g2,edge1,edge2)
                        if edgematch==True:
                            root_match = match_node_to_node(g1,g2,matches,u)
                            if v in [x[1] for x in root_match] and u not in [x[0] for x in correspondances]:
                                correspondances.append([u,v])
                                success=True

    return correspondances,success

def fits_current_corr(g1,g2,corr,node1,node2):
    edges_to_check_1 =  []
    edges_to_check_2 =  []
    for i in g1.neighbors(node1):
        for j in corr:
            if i in j and (node1,j[0]) not in edges_to_check_1:
                matcher = j[1]
                edges_to_check_1.append((node1,j[0]))
                edges_to_check_2.append((node2,matcher))
    is_match = True

    for i in range(len(edges_to_check_1)):
        if g1.get_edge_data(*edges_to_check_1[i])!=g2.get_edge_data(*edges_to_check_2[i]):
            #print(edges_to_check_1[i],edges_to_check_2[i])
            #print(g1.get_edge_data(*edges_to_check_1[i]))
            #print(g2.get_edge_data(*edges_to_check_2[i]))
            is_match = False
    return is_match



def match_node_to_node(g1,g2,matches,corr,key_node):
    #print("MATCH NODE TO NODE")
    #print(key_node)
    root_match = []
    if key_node not in [x[0] for x in corr]:
        for i in matches[key_node]:
            #print("TESTING",str(i))
            nodes_are_aligned=False
            if i not in [x[1] for x in corr]:
                nodes_are_aligned=True
                edges1 = get_edges_at_node(g1,key_node)
                edges2 = get_edges_at_node(g2,i)
                for j in edges1:
                    node_match = False
                    edge1 = get_edge_between_two_nodes(g1,key_node,j)
                    for k in edges2:
                        edge2 = get_edge_between_two_nodes(g2,i,k)
                        #print("TRYING TO MATCH TWO EDGES:",edge1,edge2)
                        edgematch = edge_match(g1,g2,edge1,edge2)
                      #  print(edgematch,j,k)
                        if edgematch==True:
                            if fits_current_corr(g1,g2,corr,key_node,i)==True:
                                node_match=True
                    if node_match==False:
                        nodes_are_aligned=False
            if nodes_are_aligned==True:
                root_match.append(i)
        #print(root_match)
        #print(corr)
        #print("ENDING NODE TO NODE")
    return (root_match,corr)

def align_iterative(g1,g2,alignable, matches,correspondances, seed = 0):
    seed_g1,seed_g2 = -1,-1
    iterator = 0
    for i in alignable.keys():
        for j in alignable[i]:
            if iterator==seed:
                seed_g1,seed_g2 = i,j
            iterator = iterator + 1

    changed = True
    #print(alignable)
    while changed==True and len(correspondances)<len(g1.nodes()):
        #print("NEW LOOP")
        aligned1=0
        aligned2=0
        changed=False
        for i in alignable.keys():
            #print("DID FIRST PART")
            if len(alignable[i])==1:
                changed=True
                correspondances.append([i,alignable[i][0]])
                aligned1 = i
                #aligned2 = alignable[i][0]
               #for k in alignable.keys():
                 #   for j in alignable[k]:
                 #       if j == aligned2 :
                 #           alignable[k].remove(j)
               # alignable.pop(aligned1,None)
        if changed==True:
            for t in alignable.keys():
                        alignable[t],correspondances = match_node_to_node(g1, g2, matches, correspondances, t)
        #print(alignable)
        #print(correspondances)
        if changed==False and len(alignable.keys())>0:
            #print("GOT HERE")
            #print(alignable)
            topop = []
            for kk in alignable.keys():
                if len(alignable[kk])==0:
                    topop.append(kk)
            for ll in topop:
                alignable.pop(ll,None)
            if seed_g1 not in [x[0] for x in correspondances] and seed_g2 not in [x[1] for x in correspondances]:
                try:
                    correspondances.append([seed_g1,seed_g2])
                except:
                    return correspondances,alignable,False
                changed=True
                aligned1 = seed_g1
                aligned2 = seed_g2
           # for k in alignable.keys():
          #      for j in alignable[k]:
          #          if j == aligned2 :
         #               alignable[k].remove(j)
           # alignable.pop(aligned1,None)
            for t in alignable.keys():
                alignable[t],correspondances = match_node_to_node(g1, g2,matches, correspondances, t)

    return correspondances,alignable,True
def align_graphs(g1,g2):
    """
    takes as input two identical DiGraph objects
    :param g1: DiGraph 1
    :param g2: DiGraph 2
    :param aln : proto-alignment dict. List of dicts
    :return: Dictionary of alignment of nodes
    """

    correspondances = []
    matches = {}
    match = DGM(g1,g2)
    for xx in g1.nodes():
        matches[xx]=[]
        for yy in g2.nodes():
            if match.syntactic_feasibility(xx,yy) == True:
                matches[xx].append(yy)
    min_match = 100
    key_node = -1
    for i in matches.keys():
        if (len(matches[i]))<=min_match:
            key_node = i
            min_match = len(matches[i])
    #align g1 on g2
    correspondances = []
    #print(matches)
    alignable = od()
    for i in g1.nodes():
        alignable[i],correspondances = match_node_to_node(g1,g2,matches,correspondances,i)
    it = 0
    seed = 0
    base_alignable = alignable.copy()
    base_corr = []

    while len(correspondances)<len(g1.nodes()):
        if len(correspondances)+ len(alignable)<len(g1.nodes()):
            #print("NEW SEED")
            alignable = od()
            correspondances= []
            for i in g1.nodes():
                alignable[i],correspondances = match_node_to_node(g1,g2,matches,correspondances,i)
            seed = seed +1
            correspondances,alignable,completed_aln = align_iterative(g1,g2,alignable,matches,[],seed)
        it = it+1
        if it>80:
            raise ValueError('A very specific bad thing happened.')

        correspondances,alignable,completed_aln = align_iterative(g1,g2,alignable,matches,correspondances,seed)
        #print(alignable)
        #print(correspondances)
        while completed_aln==False:
            seed = seed +1
            correspondances,alignable,completed_aln = align_iterative(g1,g2,base_alignable,matches,correspondances,seed)
            #print(alignable)
            #print(correspondances)
    return correspondances

def compare_graphs(g1,g2):
    em = iso.generic_edge_match('label',["CWW","TWW","CWS","TWS","CWH","TWH","CHH","THH","THW","CHW","CHS","THS","CSW","TSW","CSH","TSH","CSS","TSS","C++","T++","C--","T--","B53"],operator.eq)

    isof = nx.is_isomorphic(g1,g2,edge_match=em)
    return isof

#t = pickle.load(open("bp2_rna3dmotif_PDB_positions.cPickle", 'rb'))

full_aln = []
PDBs = []
all_graphs = []
motif_pos = []

modules = pickle.load(open("all_modules_all_examples.cPickle","rb"))
graphs, ids, positions = modules

print(graphs[0])
print(ids[0])
print(positions[0])
print(len(graphs))

for ind, mod in enumerate(graphs):
    mod_graphs = []
    mod_positions = []
    mod_PDB_names = []
    if len(mod)<9:
        continue
    for occ_ind, occ in enumerate(mod):
        module_seq, module_ss, chain = get_seq_ss(ids[ind][occ_ind])
        if len (module_seq) < 1:
            #print("rejected for PDB not found")
            continue

        motif_seq = "".join([x[1]["nuc"] for x in sorted(graphs[ind][occ_ind].nodes(data=True))])
        real_motif_seq = ""
        if positions[ind][occ_ind][-1] < len(module_seq):
            ##print("testing position")
            real_motif_seq = "".join([module_seq[x - 1] for x in positions[ind][occ_ind]])
        if motif_seq != real_motif_seq:
            #print("rejected for wrong position")
            continue

        not_boring=False
        j = occ
        for edge in j.edges(data=True):
            if edge[2]["label"].upper() not in ["CWW","B53","S55","S35","S53","S33"]:
                not_boring=True
        if not_boring==False:
            #print("rejected for boring,", j.edges(data=True))
            continue

        size = len(list(j.nodes()))
        if size < 5 or size > 25:
            #print("rejecte for size:",size)
            continue

        last = sorted(positions[ind][occ_ind])[-1]
        first = sorted(positions[ind][occ_ind])[0]
        offset = max(0, first - 100)
        if last - first < 200:
            new_seq = module_seq[max(0, first - 100):last + 100]
            if len(new_seq) > 250 or len(new_seq) < 10:
                continue
        else:
            continue



        #print("accepted")
        mod_graphs.append(occ)
        mod_positions.append((ids[ind][occ_ind],positions[ind][occ_ind]))
        mod_PDB_names.append(ids[ind][occ_ind])
    print("SIZE OF MODULE ", ind,len(mod_graphs), "OF", len(graphs[ind]))
    if len(mod_graphs)>2:
        try:
            bobby = align_multiple_graphs(mod_graphs)
            full_aln.append(bobby)
            all_graphs.append(mod_graphs)
            PDBs.append(mod_PDB_names)
            motif_pos.append(mod_positions)
            # gg = make_align_graph(j,t,corr)
        except ValueError as e:
            print("TERRIBLE")
            continue
pickle.dump(full_aln, open("bp2_rna3dmotif_aligned_modulegraphs.cPickle", 'wb'))
pickle.dump(all_graphs, open("bp2_rna3dmotif_one_of_each_graph.cPickle", 'wb'))
pickle.dump(PDBs, open("bp2_rna3dmotif_PDB_names.cPickle", 'wb'))
pickle.dump(motif_pos, open("bp2_rna3dmotif_PDB_positions.cPickle", 'wb'))


print("N MODULES")
print(len(all_graphs))

