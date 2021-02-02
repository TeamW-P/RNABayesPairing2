#import sys
#sys.path.append('/home/mcb/rsarra2/anaconda3/envs/py35/lib/python3.6/site-#packages')
#sys.path.append('/home/mcb/rsarra2/usr/local/lib')
#sys.path.append('/home/mcb/rsarra2/usr/local/bin')
import treedecomp
import pickle
import sys
import networkx as nx

def node_order_map(g):
    nmap = {}
    snodes = sorted(g.nodes())
    for ind,node in enumerate(snodes):
        nmap[node]=ind
    return nmap

def get_dependencies(module, map_index, out=None):
    """
    Return deoendencies list for a given module
    """
    inv_map = {v: k for k, v in map_index.items()}

    G = nx.Graph()
    G.add_nodes_from(range(len(module.nodes)))
    G.add_edges_from([(map_index[i], map_index[j]) for i,j in module.edges()])
    td = treedecomp.TreeDecompositionFactory().create(len(G.nodes), G.edges)
    if out:
         td.writeTD(out)
    res = {}
    for idx in td.toposorted_bag_indices():
        children = td.adj[idx]
        x = td.get_bags()[idx]
        for idy in children:
            y = td.get_bags()[idy]
            key = td.diff_set(x,y)[0]
            value = td.sep_set(x,y)
            res[inv_map[key]] = [inv_map[x] for x in value]
    return res

def tree_decomposition(g):
    return get_dependencies(g,node_order_map(g))

if __name__ == "__main__":
    map_index = {1797:0, 1798:1, 1799:2, 1800:3, 1801:4, 1802:5, 1803:6}
    module = pickle.load(open("test_module", "rb"))
    dependencies = get_dependencies(module, map_index)
    print(dependencies)
    #args = sys.argv
    #sequence = sys.argv[1]
    #dataset = sys.argv[2]
    #module = int(sys.argv[3])
    #BN = makeBN.call_makeBN(module, dataset,"NONE",False,"", 10, 0.2)
    #score = jared(sequence,BN)
    #print(score)
