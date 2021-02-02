from itertools import product

def is_parenthesized(inst):
    """
    Returns True if the given bracket-dot presented secondary structure is well parenthesized
    Otherwuse, returns False
    """
    n = 0
    for c in inst:
        if c == '(':
            n += 1
        elif c == ')':
            n -= 1
    return n == 0

def bracket_to_index(inst):
    """
    For a given bracket-dot presented secondary structure S, the function returns an iterger list L.
    L[i] = j if i and j are paired in S.
    """
    res = [-1]*(len(inst)+2)
    tmp = []
    for i,c in enumerate("("+inst+")"):
        if c == '(':
            tmp.append(i)
        elif c == ')':
            j = tmp.pop()
            res[i], res[j] = j, i
    return res

def decomposition(inst):
    """
    Decompose a given bracket-dot presented RNA secondary structure into several basic components 
    in tree-presentation.
    A basic component is presented by a list of its paired bases positions
    """
    index = bracket_to_index(inst)
    def aux(ind):
        """
        A recursive function decomposing a given RNA secondary structure in index list 
        from a given starting position.
        """
        tmp = []
        res = []
        k = ind+1
        while True:
            if index[k] == -1 :
                k += 1
            elif index[k] > k:
                tmp.append((k, index[k]))
                res += aux(k)
                k = index[k]+1
            else:
                break
        return [[(ind,index[ind])]+tmp]+res

    res = aux(0)
    return res

def find_minimal(stc):
    """
    Find the set of minimal structures for a given set of RNA secondary structures
    """
    lst = sorted(stc, key=len)
    res = []
    for s in lst:
        try:
            for t in res:
                assert not t in s
            res.append(s)
        except:
            pass
    return res

def gen_ext(n, allowed=3, paired=True):
    res = {}
    if paired:
        m = n-1
    else:
        m = n+1
    for ind in range(m):
        if ind >= allowed+2:
            res[ind] = ["."+s for s in res[ind-1]]
            for i in range(allowed,ind-1):
                res[ind] += ["("+s+")"+t for s,t in product(res[i],res[ind-2-i])]
        else:
            res[ind] = ["."*ind]
    if paired:
        for k in res.keys():
            res[k] = list(map(lambda t: "("+t+")", res[k]))
    return set().union(*res.values())

if __name__ == '__main__':
    s = '(((.((...))((...)))))'
    assert bracket_to_index(s) == [22,21,20,19,-1,11,10,-1,-1,-1,6,5,18,17,-1,-1,-1,13,12,3,2,1,0]
    print(decomposition(s))


