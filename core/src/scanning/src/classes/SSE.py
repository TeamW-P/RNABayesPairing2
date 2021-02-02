from collections import namedtuple
from anytree import NodeMixin
from ..utils.structure import decomposition

BasePair = namedtuple('BasePair', ('fst','snd'))

class SSE:
    """
    A class for basic components of RNA secondary structure, including Stacking Pair, Hairpin Loop,
    Internal Loop and Multiple Loop
    """

    def __init__(self, lst, root=False, seq=""):
        self.adjacent = []
        bp_list = list(map(lambda t: BasePair(*t), lst))
        nb_bp = len(bp_list)
        self.closing_bp = bp_list[0]
        self.enclosed_bp = bp_list[1:]
        self.get_strands_len()
       # print("BP LIST",self.closing_bp, self.enclosed_bp, self.strands_len)
        self.seq=seq

        # Hairpin Loop
        if nb_bp == 1:
            self.name = 'Hairpin Loop'
            self.label = 'H'
        
        # Stacking Pair and Internal Loop
        elif nb_bp == 2:
            bp_1 = self.closing_bp
            bp_2 = self.enclosed_bp[0]
            if (abs(bp_1.fst-bp_2.fst) == 1) and (abs(bp_1.snd-bp_2.snd) == 1):
                self.name = 'Stacking Pair'
                self.label = 'S'
            else:
                self.name = 'Internal Loop'
                self.label = 'I'
            
        # Multiple Loop
        elif nb_bp > 2:
            self.name = 'Multiple Loop'
            self.label = 'M'

        # The base is the root of tree-presented structure
        if root:
            self.name = 'Root'
            self.label = 'R'

    @classmethod
    def from_bracket(cls, structure):
        index = decomposition(structure)
        return cls(index[1])
    
    def get_strands_len(self):
        """
        Return the length for each strand in the secondary structure element (SSE).
        For example, two strands for the internal loop.
        """
        try:
            return self.strands_len
        except:
            self._update_strands_info()
            return self.strands_len

    def get_positions(self):
        try:
            return self.positions
        except:
            self._update_strands_info()
            return self.positions



    def _update_strands_info(self):
        tmp = BasePair(self.closing_bp.snd,self.closing_bp.fst)
        l = [tmp]+self.enclosed_bp+[tmp]
        length = []
        pos = []
        for i in range(len(l)-1):
            fst = l[i+1].fst
            snd = l[i].snd
            if snd>fst:
                continue
            length.append(fst-snd+1)
            pos.append(list(range(snd,fst+1)))
        self.positions = pos
        self.strands_len = tuple(length)

        ## Correction for stackings





    def get_sequence(self):
        try:
            return self.seq
        except:
            pass

    def is_same_form(self, sse):
        """
        Return True if the given SSE has the same form (in structure)
        More precisely, the functions return True if two SSE have same number of bases of each strands
        """
        if self.get_strands_len() == sse.get_strands_len():
            return True

        s1 = list(self.get_strands_len())
        s2 = list(sse.get_strands_len())

        total_dist = 0
        for strand in range(len(s1)):
            total_dist += (s1[strand]-s2[strand])
        #print("TOTAL DIST",total_dist)
        if total_dist<=len(s1):
            return True


        return False

    def __str__(self):
        return self.name

class SSENode(SSE, NodeMixin):

    def __init__(self, lst, root=False, parent=None, seq=""):
        super().__init__(lst, root=root, seq=seq)
        self.parent = parent

    def __str__(self):
        return self.name


