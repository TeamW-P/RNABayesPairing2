from anytree import RenderTree
from .SSE import SSENode
from ..utils.structure import bracket_to_index


class SSETree:

    def __init__(self):
        self.root = None

    @classmethod
    def from_bracket(cls, structure, seq=""):
        tree = cls()
        tree.structure = structure
        tree.seq = seq
        # Notice that index 0 and len+1 represent the fake basepair
        index = bracket_to_index(structure)
        seq = "#"+seq+"#"
        # Recursive function to build tree decomposition and to assign sequence
        # Minus 1 for all indices, so index -1 and len for the fake basepair
        def aux(ind, parent=None):
            tmp = []
            k = ind+1
            j = max(ind,1)
            res = []
            w = []
            while True:
                pair_k = index[k]
                # If the base is non-paired, move on
                if pair_k == -1:
                    k += 1
                # If it's an open parenthese, call the recursive function
                elif pair_k > k:
                    tmp.append((k-1,pair_k-1))
                    if seq:
                        w.append(seq[j:k+1])
                    res.append(aux(k))
                    j = pair_k
                    k = pair_k + 1
                # If it's a close parenthese correspond to ind, add the sequence then break
                else:
                    if seq:
                        w.append(seq[j:min(k+1,len(seq)-1)])
                    break
            node = SSENode([(ind-1,index[ind]-1)]+tmp, root=(not ind), seq="&".join(w))
            #print("SSE creation:",node.positions)
            for s in res:
                s.parent = node
            return node
        tree.root = aux(0)
        return tree

    def print_tree(self):
        """
        A simple tree printer
        """
        for pre,_,node in RenderTree(self.root):
            treestr = u"%s%s" % (pre,node.label)
            print(treestr.ljust(0),node.get_strands_len())

    def get_sequence(self):
        if tree.seq:
            return tree.seq
        else:
            raise Exception("No assigned sequence")

if __name__ == "__main__":
    s = "(((.((...))((...)))))"
    tree = SSETree.from_bracket(s)
    tree.print_tree()
