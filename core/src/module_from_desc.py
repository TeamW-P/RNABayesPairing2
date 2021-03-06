import os
import networkx as nx
import matplotlib
from matplotlib import pyplot as plt
import networkx.algorithms.isomorphism as iso
import operator
import os.path
import pickle
import argparse
import json
from Bio import SeqIO
import os

CURRENT_DIRECTORY = os.path.dirname(__file__)

def compare_graphs(g1,g2):
    em = iso.generic_edge_match('label',["cWW","tWW","cWS","tWS","cWH","tWH","cHH","tHH","tHW","cHW","cHS","tHS","cSW","tSW","cSH","tSH","cSS","tSS","c++","t++","c--","t--","b53"],operator.eq)

    isof = nx.is_isomorphic(g1,g2,edge_match=em)
    return isof

def cluster(desc_file, g, graphs, IDs):
    if len(graphs) == 0:
        graphs.append([g])
        IDs.append([desc_file[:-5]])
    else:
        found = False
        kk = 0
        while kk < len(graphs) and found == False:
            if compare_graphs(g, graphs[kk][0]) == True:
                graphs[kk].append(g)
                IDs[kk].append([desc_file[:-5]])
                found = True
            kk = kk + 1
        if found == False:
            graphs.append([g])
            IDs.append([desc_file[:-5]])
    return (graphs, IDs)

def make_align_graph(g1):
    corr_edges = []
    g0 = nx.DiGraph()
    g0.add_nodes_from(g1.nodes())
    g0.add_edges_from(g1.edges())
    pos = nx.spring_layout(g0)
    nx.draw_networkx_nodes(g0,pos,nodelist=g1.nodes(),node_color='red',node_size=500)
    nx.draw_networkx_edges(g0,pos,edgelist=g1.edges(),edge_color='red',width=2)
    labels={}
    elabels = {}
    for i in g0.nodes():
        labels[i] = i
    nx.draw_networkx_labels(g0,pos,labels)
    for i in g1.edges():
        elabels[i]=g1.get_edge_data(*i)['label']
    nx.draw_networkx_edge_labels(g0,pos,elabels)
    plt.show()

def make_graph(desc_file):
    g = nx.DiGraph()
    with open(desc_file, 'r') as desc:
        lines = desc.readlines()
        nodestuff = lines[1][:-1].split(" ")
        # print(nodestuff)
        for n in nodestuff:
            if len(n) > 0:
                if n[0].isdigit():
                    node_n = n[:-2]
                    nuc = n[-1]
                    p_id = 0
                    if len(g.nodes())==0:
                        g.add_node(int(node_n), nuc= nuc, part_id=0)
                    else:
                        if int(node_n) - int(list(g.nodes())[-1])<3:
                            p_id = list(g.nodes(data=True))[-1][1]['part_id']
                        else:
                            p_id = int(list(g.nodes(data=True))[-1][1]['part_id'])+1
                        g.add_node(int(node_n), nuc=nuc, part_id=p_id)
        for line in range(2, len(lines)):
            bp = lines[line].split("---")
            p1 = bp[0][4:8].replace(" ", "")
            p2 = bp[2][2:7]

            interaction = bp[1][1:4]
            orientation = bp[1][5]
            if orientation == 's':
                bond = orientation.upper() + interaction[0] + interaction[2]
                g.add_edge(int(p1), int(p2), long_range=False, label= bond)
            elif interaction == "C/C":
                g.add_edge(int(p1), int(p2), long_range = False, label = "B53")
            elif interaction == "+/+" or interaction == "-/-":
                bond = orientation.upper() + "WW"
                g.add_edge(int(p1), int(p2), long_range= False, label = bond)
            else:
                bond = orientation.upper() + interaction[0] + interaction[2]
                g.add_edge(int(p1), int(p2), long_range= False, label = bond)
    return g


                    
def make_new_graph_examples(g,fasta, pos_signature):
    SEQUENCES = []
    graphs = []
    with open(fasta,"r") as f:
        lines = f.readlines()
        for line in lines:
            if line[0]!=">":
                seq = list(line)
                del seq[-1]
                SEQUENCES.append(seq)
                seq_dict={}
                h = nx.DiGraph()
                node_correspondence = {}
                
                        
                        


                for ind, n in enumerate(pos_signature):
                    h.add_node(n,nuc=seq[ind])
                #print(list(h.nodes()))
                #print(list(g.nodes()))

                for indn,node in enumerate(list(g.nodes())):
                    node_correspondence[node] = list(h.nodes())[indn]
                for inde,edge in enumerate(list(g.edges(data=True))):
                    #print(edge[2]["label"])
                    #print(node_correspondence)
                    
                    start,end,info = edge
                    #print(start,end)
                    new_start = node_correspondence[start]
                    new_end = node_correspondence[end]
                    #print(new_start,new_end,info)
                    if new_start>new_end:
                        h.add_edges_from([(new_end,new_start,info)])
                    else:
                        h.add_edges_from([(new_start,new_end,info)])
                        
                                          

                
                graphs.append(h)
                #print("nodes",h.nodes(data=True))
    return graphs,(SEQUENCES,[])

def get_all_signatures_from_graphs(gs):
    sigs = []
    
    for g in gs:
        positions = g
        new_nodes = []
        fnode = min(positions)
        current_strand = 0
        current_strand_position = 0
        this_sig = []
        for ind,position in enumerate(positions):
            if ind>0:
                if position!=(positions[ind-1]+1):
                    diff = abs(position-positions[ind-1])
                    if diff>4:
                        current_strand+=1
                        current_strand_position=0
                    else:
                        for zz in range(diff-1):
                            current_strand_position += 1
            this_sig.append(current_strand_position+100*current_strand)
            current_strand_position+=1
        if this_sig not in sigs:
            sigs.append(this_sig)
    return sigs
    
    
def create_module(desc,fasta,model_name, list_of_nodes, PDBs=[]):
    graphs_name = os.path.join(CURRENT_DIRECTORY, "../models/" +model_name+"_one_of_each_graph.cPickle")
    aln_name = os.path.join(CURRENT_DIRECTORY, "../models/"+model_name+"_aligned_modulegraphs.cPickle")
    PDB_name = os.path.join(CURRENT_DIRECTORY, "../models/"+model_name+"_PDB_names.cPickle")
    seq_name = os.path.join(CURRENT_DIRECTORY, "../models/"+model_name+"_sequences.pickle")
    siblings =  os.path.join(CURRENT_DIRECTORY, "../models/"+model_name+"_siblings.pickle")
    repeat= os.path.join(CURRENT_DIRECTORY, "../models/repeat.pickle")
    if 'desc' in desc.lower():
        g  = make_graph(desc)
    else:
        g = pickle.load(open(os.path.join(CURRENT_DIRECTORY, desc),'rb'))
    if os.path.isfile(graphs_name):
        number=pickle.load(open(os.path.join(CURRENT_DIRECTORY, graphs_name),'rb'))
    else:
        number=[]
    moduleSignatures = get_all_signatures_from_graphs(list_of_nodes)
    print("THERE ARE",len(moduleSignatures),"MODULES SIGNATURES:")
    print(moduleSignatures)
    
    sib_sigs = list(range(len(number),len(number)+len(moduleSignatures)))
    
    if os.path.isfile(siblings):
        sibDict = pickle.load(open(siblings,'rb'))
    else:
        sibDict = {}
    for s in sib_sigs:
        if s in sibDict:
            for t in sib_sigs:
                if t!=s:
                    sibDict[s].append(t)
        else:
            sibDict[s] = [x for x in sib_sigs if x!=s]
            
    print("SIBLING MODULES DICTIONARY")
    print(sibDict)
        
    
    try:
        seqs = pickle.load(open(os.path.join(CURRENT_DIRECTORY, seq_name),'rb'))
    except:
        seqs = []

    if os.path.isfile(graphs_name):
        graphs = pickle.load(open(os.path.join(CURRENT_DIRECTORY, graphs_name),"rb"))
        aln = pickle.load(open(os.path.join(CURRENT_DIRECTORY, aln_name),"rb"))
        PDB_list = pickle.load(open(os.path.join(CURRENT_DIRECTORY, PDB_name),"rb"))
    else:
        graphs=[]
        aln = []
        PDB_list = []

    for sig in moduleSignatures:
        
        print("sig is",sig)    
        graph_list,sequences = make_new_graph_examples(g, fasta, sig)
        
       
        
                
    
        new_aln={}
        for graph in range(len(graph_list)):
            for position in list(graph_list[0].nodes()):
                #print("position",position)
                if position in new_aln:
                    new_aln[position].append(position)
                else:
                    new_aln[position]= [position]
            

        if PDBs == []:
            PDBs = ["None" for x in graph_list]

        graphs.append(graph_list)
        PDB_list.append(PDBs)
        seqs.append(sequences)
        aln.append(new_aln)
    
    pickle.dump(graphs, open(os.path.join(CURRENT_DIRECTORY, graphs_name), 'wb'))
    pickle.dump(aln, open(os.path.join(CURRENT_DIRECTORY, aln_name), 'wb'))
    pickle.dump(PDB_list, open(os.path.join(CURRENT_DIRECTORY, PDB_name), 'wb'))
    pickle.dump(seqs,open(os.path.join(CURRENT_DIRECTORY, seq_name), "wb"))
    pickle.dump(sibDict,open(os.path.join(CURRENT_DIRECTORY, siblings), "wb"))
    pickle.dump(len(moduleSignatures),open(os.path.join(CURRENT_DIRECTORY), repeat, "wb"))
    print("YOUR MODULE WAS ADDED TO DATASET ", model_name, "AND RECEIVED THE NUMBERS",len(graphs)-1-len(moduleSignatures),"TO", len(graphs)-1, " WITH ",len(sequences[0])," SEQUENCES")
    return graphs,aln,sequences

def create_module(desc, fasta, model_name, list_of_nodes, PDBs=[], pdb_info_file="", full_seq_file="", atlas_name=""):
    repeat = os.path.join(CURRENT_DIRECTORY, "../models/repeat.pickle")
    dataset_name = os.path.join(CURRENT_DIRECTORY, "../models/" + model_name + ".json")

    dataset = {}
    # load the dataset if it exists
    if os.path.isfile(dataset_name):
        with open(dataset_name) as f:
            dataset = json.load(f)

    # load in the graph
    if 'desc' in desc.lower():
        g = make_graph(desc)
    else:
        g = pickle.load(open(desc, 'rb'))

    pdb_aln = []
    pdb_seq = []
    # Try to load PDB_info
    if os.path.isfile(pdb_info_file):
        with open(pdb_info_file, 'rb') as f:
            pdb_aln, pdb_seq = pickle.load(f)

    # try to load full_seq_file
    full_records = []
    if os.path.isfile(full_seq_file):
        full_records = list(SeqIO.parse(full_seq_file, "fasta"))

    # module signatures from list of nodes
    module_signatures = get_all_signatures_from_graphs(list_of_nodes)
    print("THERE ARE", len(module_signatures), "MODULES SIGNATURES:")
    print(module_signatures)

    dataset_length = len(dataset)

    # for each signaature, add the module
    for sig in module_signatures:

        print("sig is", sig)
        graph_list, sequences = make_new_graph_examples(g, fasta, sig)

        # this gives alignment of graphs
        new_aln = {}
        for graph in range(len(graph_list)):
            for position in list(graph_list[0].nodes()):
                # print("position",position)
                if position in new_aln:
                    new_aln[position].append(position)
                else:
                    new_aln[position] = [position]

        if PDBs == []:
            PDBs = ["None" for x in graph_list]

        data = {}
        data["master_graph"] =  {"nodes":list(graph_list[0].nodes(data=True)), "edges": list(graph_list[0].edges(data=True))}

        # list comprehension could also work but we'll use a for loop for readability
        data["graphs"] = {"nodes":[], "edges": []}
        for g in graph_list:
            data["graphs"]["nodes"].append(list(g.nodes(data=True)))
            data["graphs"]["edges"].append(list(g.edges(data=True)))

        # alignment data of graphs
        data["aln"] = new_aln

        # PDB names and positions of the motif
        data["PDBs"] = dict(pdb_aln)
        # PDB names and the sequences of the motif
        data["subsequences"] = dict(pdb_seq)

        data["training_set"] = []
        # training sequences
        if len(list_of_nodes) == len(full_records):
            # adding training sequences
            for index, record in enumerate(full_records):
                data["training_set"].append({
                    "source": record.id,
                    "seq": str(record.seq),
                    "source_pos": [],
                    "seq_pos": list_of_nodes[index]
                })

        # siblings
        sib_sigs = list(range(dataset_length, dataset_length + len(module_signatures)))

        # in sib dict, just remove self
        module_num = len(dataset)
        sib_sigs.remove(module_num)
        data["siblings"] = sib_sigs

        # add atlas name
        if atlas_name:
            data["atlas_name"] = atlas_name
        else:
            data["atlas_name"] = ""

        # general motif name, we will populate this later
        data["general_name"] = ""

        dataset[len(dataset)] = data

    print("len", len(dataset))
    # dump the json
    with open(dataset_name, 'w') as out:
        json.dump(dataset, out)

    print("YOUR MODULE WAS ADDED TO DATASET ", model_name, "AND RECEIVED THE NUMBERS", len(dataset) - len(module_signatures), "TO", len(dataset) - 1, " WITH ", len(sequences[0]), " SEQUENCES")


def fuse_existing_databases(new_dataset,name1,name2,l1,l2):
    one_of_each_graph = []
    a1 = pickle.load(open(os.path.join(CURRENT_DIRECTORY, name1 + "_one_of_each_graph.cPickle"), 'rb'))
    a2 = pickle.load(open(os.path.join(CURRENT_DIRECTORY, name2 + "_one_of_each_graph.cPickle"), 'rb'))
    for ind,i in enumerate(a1):
        if ind in l1:
            one_of_each_graph.append(i)
    for ind,i in enumerate(a2):
        if ind in l2:
            one_of_each_graph.append(i)
    pickle.dump(one_of_each_graph,open(os.path.join(CURRENT_DIRECTORY, "../models/"+new_dataset+"_one_of_each_graph.cPickle"), "wb"))

    aligned_modulegraphs = []
    a1 = pickle.load(open(os.path.join(CURRENT_DIRECTORY, name1 + "_aligned_modulegraphs.cPickle"), 'rb'))
    a2 = pickle.load(open(os.path.join(CURRENT_DIRECTORY, name2 + "_aligned_modulegraphs.cPickle"), 'rb'))
    for ind,i in enumerate(a1):
        if ind in l1:
            aligned_modulegraphs.append(i)
    for ind,i in enumerate(a2):
        if ind in l2:
            aligned_modulegraphs.append(i)
    pickle.dump(aligned_modulegraphs,open(os.path.join(CURRENT_DIRECTORY, "../models/"+new_dataset+"aligned_modulegraphs.cPickle"), "wb"))


    PDB = []
    a1 = pickle.load(open(os.path.join(CURRENT_DIRECTORY, name1 + "PDB.cPickle"), 'rb'))
    a2 = pickle.load(open(os.path.join(CURRENT_DIRECTORY, name2 + "PDB.cPickle"), 'rb'))
    for ind,i in enumerate(a1):
        if ind in l1:
            PDB.append(i)
    for ind,i in enumerate(a2):
        if ind in l2:
            PDB.append(i)
    pickle.dump(PDB,open(os.path.join(CURRENT_DIRECTORY, os.path.join(CURRENT_DIRECTORY, "../models/"+new_dataset+"PDB.cPickle"), "wb")))


    PDB_pos = []
    a1 = pickle.load(open(os.path.join(CURRENT_DIRECTORY, name1 + "PDB_pos.cPickle"), 'rb'))
    a2 = pickle.load(open(os.path.join(CURRENT_DIRECTORY, name2 + "PDB_pos.cPickle"), 'rb'))
    for ind,i in enumerate(a1):
        if ind in l1:
            PDB_pos.append(i)
    for ind,i in enumerate(a2):
        if ind in l2:
            PDB_pos.append(i)
    pickle.dump(PDB_pos,open(os.path.join(CURRENT_DIRECTORY, "../models/"+new_dataset+"PDB_pos.cPickle"), "wb"))


import ast
if __name__ == "__main__":
    #desc_file = "../models/1A1T.B.1.desc"
    #fasta_file = "../models/input.fasta"
    #graphs = create_module(desc_file,fasta_file)
    #fuse_existing_databases("all_carnaval",'carnaval','jeffrey',[0, 5, 8, 26, 6, 18, 19, 21, 24, 25, 39, 53, 54, 58, 119, 122, 135, 146, 162, 168, 191, 194, 198, 216],[1,4,7,10,11,24,33,40,55,66,69,101,103,113,149,156,167])
    arguments = {}
    parser = argparse.ArgumentParser()
    #parser.add_argument("--verbose", help="increase output verbosity", action="store_true")
    parser.add_argument("-g", help="graph of the module, DESC file", required = True)
    parser.add_argument("-seq", help="sequences, FASTA format",required=True)
    parser.add_argument("-n", help="Dataset name (create new or add to existing)", required=True)
    parser.add_argument("-nodes", help="List of lsit of nodes as a string", required=True)

    parser.add_argument('-pdb', nargs='*', help='PDBs in which input is found')
    parser.add_argument('-pdb_info', help='Pickle file with tuple of PDB module positions and nucleotides for each column')
    parser.add_argument("-full_seq", help="full sequences, FASTA format")
    parser.add_argument("-atlas_name", help="name of the motif in the 3dMotifAtlas")


    args = parser.parse_args()
    create_module(args.g,args.seq, args.n, ast.literal_eval(args.nodes), args.pdb, pdb_info_file=args.pdb_info, full_seq_file=args.full_seq, atlas_name=args.atlas_name)

