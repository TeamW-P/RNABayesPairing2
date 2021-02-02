from scanning.src.classes import SSE, ExactScanResult, FuzzyScanResult
from anytree.iterators.preorderiter import PreOrderIter
import itertools

## import exact module database, a dictionary, in which the key is class representation and value is the module indices list
#MODULES = {(5,): [1,2], (3,2,2):[5]}

def naive_scanning(tree, module):
    """
    For a given structure and a module, return the list of SSEs that matches with the module in the terms of structure.
    """
    scan = lambda t: SSE.is_same_form(t, SSE.from_bracket(module))
    return PreOrderIter(tree.root, filter_=scan)

def exact_matching(tree,db,fuzzy=False):
    """
    For a given RNA secondary structure in tree presentation, add exact matching module list to each node
    """
    MODULES=db
    for node in PreOrderIter(tree.root):
        node.modules = {"exact":[], "fuzzy":[]}
        ind = node.get_strands_len()
        #print("SSE",node.label,node.strands_len, node.seq, node.closing_bp)
        node.modules["exact"] = MODULES.get(ind,[])
        #with open('found_junctions_tpp.txt','a') as f:
        #    if len(node.strands_len)>=3:
        #        to_write = str(node.strands_len) + str(node.closing_bp) + "\n"
        #        f.write(to_write)
        #        #print(str(node.strands_len) + str(node.closing_bp))
        #TODO: option to not search for fuzzy

        # Fuzzy search for hairpin
        if fuzzy==True:
            #print("TRYING OUT FUZZY MATCHES")
            acceptable_strands = []
            for strand in node.strands_len:
                acceptable_strands.append((-1,0,1))
            #get all teh combinations of -1,0,1 for each strand to get all the acceptable sizes

            #note: ne pas chercher les lists dont la somme fuzz est 0
            all_acceptable_strands_combinations = list(itertools.product(*acceptable_strands))
            for fuzz in all_acceptable_strands_combinations:
                if sum(fuzz)>0:
                    result = MODULES.get(tuple([ind[x]+fuzz[x] for x in range(len(fuzz))]))
                    if result is not None:
                        node.modules["fuzzy"].append(tuple((result,fuzz)))
        #print("ALL MATCHES", node.modules)

def yield_matching(tree, fuzzy=False):
    for node in PreOrderIter(tree.root):
        for m in node.modules["exact"]:
            match, rot = m
            #if rot==0:
            #    continue
            yield ExactScanResult(tree,node,match,rot)
        if fuzzy:
            for m in node.modules["fuzzy"]:
                if m is None:
                    continue
                #print("YIELDING MATCH M",m,"from",node.modules,node.modules["fuzzy"])
                # might have a format issue here
                match, rot = m[0][0]
                fuzzform = m[1]
                yield FuzzyScanResult(tree,node,match,rot,fuzzform)



if __name__ == "__main__":
    s = "(((.((...))((...)))))"
    m = "(.((*))((*)))"
    for node in naive_scanning(s,m):
        print(node.name, node.closing_bp)
