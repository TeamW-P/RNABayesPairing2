import pickle
import networkx as nx
from matplotlib import pyplot as plt
from draw_single_graph import draw_graph


def bps_in_module(graph):
    contained_edges = ["None"]
    for e in graph.edges(data=True):
        if e[2]["label"] not in contained_edges:
            contained_edges.append(e[2]["label"])
    return contained_edges
def bp_intersect(bp_to_find,bp_in_module):
    #print(bp_to_find,bp_in_module)
    allin = True
    for i in bp_to_find:
        if i not in bp_in_module:
            allin=False
    #print(allin)
    return allin

def number_of_strands(nodelist):
    n_breaks = 0
    for ind,node in enumerate(nodelist):
        if ind==0:
            continue
        if node> nodelist[ind-1]+3:
            n_breaks = n_breaks + 1
    return (n_breaks+1)

def match_pos(input_pos,module_pos):
    if len(input_pos)==0:
        return True

    for i in input_pos:
        if i in module_pos:
            return True
    return False


def get_carnaval_graph(file):
    n = []
    e = []
    with open(file) as f:
        g = nx.DiGraph()
        lines = f.readlines()
        for line in lines:
            if "+++" in line:
                continue
            else:
                n1,nuc1,bp, n2, nuc2 = line.replace("\n","").split(" ")
                if (int(n1),{"nuc": nuc1 }) not in n:
                    n.append((int(n1),{"nuc": nuc1 }))
                if (int(n2),{"nuc": nuc2 }) not in n:
                    n.append((int(n2),{"nuc": nuc2 }))
                if (int(n1),int(n2),{"label" : bp}) not in e:
                    e.append((int(n1),int(n2),{"label" : bp}))
    n = sorted(n)
    e = sorted(e)
    g.add_nodes_from(n)
    g.add_edges_from(e)
    return g


def parse_dataset_for_modules(dataset,source_PDBs, n_strands=-1,bp=["None"],input_positions=[]):
    modules_in_PDB = []
    full_ids = []
    for ind,i in enumerate(dataset):
        module_graph_ex = i[0][0]
        PDBlist = list(i[1])
        PDBs = []
        for struct in PDBlist:
            if len(struct)==1:
                PDBs.append(struct[0].split(".")[0])
            else:
                PDBs.append(struct.split(".")[0])
        for pdb_name in source_PDBs:
            if pdb_name in PDBs:
                module_g = i[0][PDBs.index(pdb_name)]
                if n_strands==-1:
                    bps_in_mod = bps_in_module(module_graph_ex)
                    intersect = bp_intersect(bp,bps_in_mod)
                    positions = list(module_graph_ex.nodes(data=False))

                    if intersect==True and match_pos(input_positions,positions)==True:
                        modules_in_PDB.append(module_g)
                        full_ids.append(i[1][PDBs.index(pdb_name)])
                else:
                    n_strands_in_module = number_of_strands(list(module_g.nodes(data=False)))
                    bps_in_mod = bps_in_module(module_graph_ex)
                    intersect = bp_intersect(bp,bps_in_mod)
                    positions = list(module_graph_ex.nodes(data=False))
                    #print(intersect,n_strands_in_module)
                    if n_strands_in_module==n_strands and intersect==True and match_pos(input_positions,positions)==True:
                        modules_in_PDB.append(i[0][PDBs.index(pdb_name)])
                        full_ids.append(i[1][PDBs.index(pdb_name)])
    return modules_in_PDB,full_ids


#dataset = pickle.load(open('../models/min0_3dmotif_graphs.cPickle', 'rb'))
#PDBs = ["2GDI", "3D2V"]
#graphs, ids = parse_dataset_for_modules(dataset, PDBs, n_strands=-1,bp=["TWH"])
#for g in range(len(graphs)):
#    draw_single_graph.draw_graph(graphs[g], ids[g])