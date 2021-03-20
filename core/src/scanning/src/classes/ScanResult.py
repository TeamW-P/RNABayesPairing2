from abc import ABC, abstractmethod
import itertools

SEQ_THRESH = 5


class ScanResult(ABC):

    @abstractmethod
    def eval(self):
        pass

    @abstractmethod
    def get_position(self, detail=True):
        pass

    @abstractmethod
    def get_request_info(self):
        """
        Return (sequence, secondary structure) of scan request
        """
        pass

    def __str__(self):
        print("Position: {} Seq: {} Score: {}".format(self.get_position(),self.seq,self.eval()))



#TODO : need a way to know when the match has been rotated!
class ExactScanResult:

    def __init__(self, tree, node, module, rotation):
        #print("EXECUTING INIT OF EXACTSCANRESULT", node.positions, module.ID, rotation)
        self.tree = tree
        self.module = module
        self.node = node

        self.rotation = rotation
        self.stackings = self.rotate_stackings()
        #print("STACKINGS",self.module.stackings)

        self.pos = self.get_position()
        #print("SSE POSITIONS",self.node.positions, self.node.get_positions(), "RETRIEVED POSITIONS",self.pos)

        self.fuzzy_form= [0]*len(self.pos)
        #print("TREE SEQ",tree.seq)
        #print("SSE SEQ",self.node.seq)
        #print("STACKINGS EXACT MATCH:",self.stackings)
        #TODO : how does this behave with rotation?
        if sum(self.stackings):
            #print("RETRIEVING SEQUENCE WITH STACKING", self.pos)
            self.seq = "&".join(map(lambda t: self.tree.seq[t[0]:t[-1]+1], self.pos))
        else:
            self.seq = node.seq
        #print("FOUND SEQUENCE:",self.seq)

        #print("UNROTATED GAPS",self.module.gaps_per_strand)
        #if rotation:
            #self.rotate(rotation)
        #print("LEN POS",self.pos,"LEN SEQ",self.seq, "ROTATION",rotation)
        #print("DONE ROTATION, CURRENT SEQUENCE", self.seq)

        #TODO : need to think about adjusting POSITIONS to gaps
        self.gaps = self.module.gaps_per_strand
        #print("ROTATED GAPS", self.gaps)
        
        #print("PREFIX SEQUENCES",self.seq)
        self.reintegrate_module_gaps_to_sse_sequence_and_positions()
        #print("Reintegrated module graps, current seq", self.seq)
        self.mapping = {}



    def rotate_stackings(self):
        n = self.rotation
        lst = self.module.stackings
        return lst[n:] + lst[:n]

    def rotate_list(self, lst):
        n = self.rotation
        return lst[n:] + lst[:n]

    def rotate(self, rotation):
        #print("PRE ROTATION SEQUENCE",self.seq,"ROTATION",rotation)
        sequence_bits = self.seq.split("&")
        new_seq = ""
        bit_counter = 0
        while bit_counter<len(self.pos):
            this_bit = sequence_bits[rotation%len(self.pos)]
            if len(new_seq)>0:
                new_seq = new_seq + "&" +  this_bit
            else:
                new_seq = new_seq  + this_bit
            bit_counter+=1
            rotation+=1
        self.seq= new_seq

        #print("POST ROTATION SEQUENCE", self.seq)


    def reintegrate_module_gaps_to_sse_sequence_and_positions(self):
        #print("REGAPPING MODULE", self.module.ss, self.module.stackings, self.gaps,self.seq,self.pos)

        positions = self.rotate_list(self.pos)
        #gaps = self.rotate_list(self.gaps) this does not help, the gaps shouldnt rotate on the module
        gaps = self.gaps
        seq = self.seq

        new_seq = ""
        if "&" not in seq:
            for ind, nuc in enumerate(seq):
                if not gaps[0]:
                    new_seq+=nuc
                else:
                    if ind not in gaps[0]:
                        new_seq += nuc
        else:
            for strand,subseq in enumerate(self.rotate_list(seq.split("&"))):
            #for strand,subseq in enumerate(seq.split("&")):
                #print("STRAND",strand,"GAPS",gaps, seq, seq.split("&"))
                if not gaps[strand]:
                    new_seq+=subseq
                else:
                    for ind, nuc in enumerate(subseq):
                        if ind not in gaps[strand]:
                            new_seq += nuc
                if strand<len(seq.split("&"))-1:
                    new_seq += "&"
        self.seq = new_seq
        #print("POST REGAPPED SEQUENCE:",new_seq)

        new_pos = []
        for strand, pos_in_strand in enumerate(positions):
            new_strand = []
            if not gaps[strand]:
                new_pos.append(pos_in_strand)
            else:
                for ind, pos in enumerate(pos_in_strand):
                    if ind not in gaps[strand]:
                        new_strand.append(pos)
                new_pos.append(new_strand)
        self.module_positions = new_pos
        #print("POST REGAPPED SEQUENCE:",new_pos,new_seq, gaps)
        #still on tryout.. does it affect anything else to name pos here immediately?
        self.pos = new_pos

    #TODO : adjust aln_file to aln_sequences and following calls to fold_compound
    def eval(self,alignment=False,sequences=[], aln_sequences="", ungapped_positions = [], fold_compound = None):
        """
        Method to calculate score
        """

        #TODO: need to make sure the gaps are managed when scoring the sequence, specifically through rotation
        #print("info eval",self.rotation,self.pos,self.module.ID,self.module.n_nodes,self.seq)
        try:
            return self.score
        except:
            if len(self.seq)<(self.module.n_nodes+len(self.stackings)-1): 
                self.score = -200
                return self.score
            if alignment==False:
                #print('EVALUATION MODULE WITH :', 'MODULE', self.module.ID, self.module.nodes , 'POSITIONS',self.pos,"SEQUENCE",self.seq,'STACKINGS',self.module.stackings,'SSE',self.node.seq, self.node.positions)
                if "#" in self.node.seq or self.pos[-1][-1]>=len(self.tree.seq) or len(self.seq)<len(self.module.nodes):
                    self.score = -1
                    return -1
               # print(self.module.ID,self.node,self.pos)
                score = self.module.eval(self.seq)[0]
                self.mapping = self.module.mapping_from_sequence(self.seq)
                #print("COMPUTING SEQUENCE SCORE",score, self.pos, self.seq)

                self.score = score
            else:
                gapped_cand  = 0
                total_score = 0
                #print("evaluating alignment score for module",self.module.ID,self.module.positions)
                for seq in sequences:
                    positions = [x  for y in self.module_positions for x in y]
                    #print("SEQ BEFORE", seq)
                    seq = [seq[ungapped_pos] for ungapped_pos in ungapped_positions]
                    #print("SEQ AFTER REMOVING GAPS", seq)
                    #print("POSITIONS:",positions)
                    tmp_seq = "".join([seq[x] for x in positions])
                    if "-" in tmp_seq:
                        #print("Candidate has a gap, no scoring", tmp_seq)
                        gapped_cand+=1
                        continue
                    else:
                        #print("Candidate scoring", tmp_seq,self.module.eval(tmp_seq)[0] )
                        total_score+= self.module.eval(tmp_seq)[0]
                        if len(self.mapping)==0:
                            self.mapping = self.module.mapping_from_sequence(tmp_seq)
                mean_score = total_score/(len(sequences)-gapped_cand)
                # print("score for this alignment:",mean_score)
                self.score = mean_score
                self.pos = self.module_positions
                return mean_score

    #double check the -1
    def get_position(self, detail=True):
        try:
            return self.pos
        except:
            #print("SELF NODE GET POSITIONS BEFORE", self.node.get_positions())
            pos = [[y for y in x] for x in self.node.get_positions()]
            #print("SELF NODE GET POSITIONS AFTER", self.node.get_positions())

            #print("MODULE POSITIONS BEFORE STACKINGS ADDED",pos,self.stackings)

            if sum(self.stackings):
                for ind, n in enumerate(self.stackings):
                    #print("CURRENT STACKING",ind,n, pos[ind])
                    pos[ind] = [pos[ind][0]-1-i for i in range(n)] + pos[ind]
                    tmp = (ind-1)%len(pos)
                    #pos[tmp] += [pos[tmp][-1]+1+i for i in range(n)]
                    pos[tmp] += [pos[tmp][-1]+1+i for i in range(n)]
                    #print("POS AFTER ITERATION",ind,n,"OF ADDING STACKINGS",pos)
                    #print("SELF NODE GET POSITIONS LOOPED 3", self.node.get_positions())
            self.pos = pos
            #print("MODULE POSITIONS AFTER STACKINGS:",pos)
            #print("SELF NODE GET POSITIONS AFTER", self.node.get_positions())

            return pos

        
    def get_request_info(self):
        return (self.tree.seq,self.tree.structure)

class FuzzyScanResult:
    #idea : take as input a SSE that is a fuzzy match to a module
    #generate all possible insertion sequences (basically 2^n_fuzzy_strands
    #Score sequences and refold the one with the highest score
    #Get score based on bayespairing1 scoring.
    def __init__(self, tree, node, module, rotation, fuzz):
        self.tree = tree
        self.module = module
        self.node = node
        self.rotation = rotation
        #the fuzzy form is in the form, say, [-1, 0] for an internal loop. -1 would bean that the module is 1 smaller than the SSE at the first strand.
        self.fuzzy_form = fuzz
        self.get_position()
        self.get_sequences()
        self.stackings = self.rotate_stackings()
        if rotation>0:
            self.rotate(rotation)
        #print("LEN POS",self.pos,"LEN SEQ",self.seq, "ROTATION",rotation)




    def rotate_stackings(self):
        n = self.rotation
        lst = self.module.stackings
        return lst[n:] + lst[:n]

    def rotate_list(self, lst):
        n = self.rotation
        return lst[n:] + lst[:n]

    #edit : a function named get_position is needed; make it output the result of eval

    def get_position(self):
        self.get_potential_module_positions_from_fuzzy_match()


    def rotate(self, rotation):
        #print("PRE ROTATION SEQUENCE",self.seq,"ROTATION",rotation)
        new_seqs = []
        for seq in self.sequences:
            sequence_bits = seq.split("&")
            new_seq = ""
            bit_counter = 0
            while bit_counter<len(self.possible_positions[0]):
                this_bit = sequence_bits[rotation%len(self.possible_positions[0])]
                new_seq = new_seq + this_bit
                bit_counter+=1
                rotation+=1
            new_seqs.append(new_seq)
        self.sequences = new_seqs
        #print("POST ROTATION SEQUENCE", self.seq)

    def get_sequences(self):
        seqs = []
        for pos in self.possible_positions:
            seqs.append("&".join(map(lambda t: self.tree.seq[t[0]:t[-1]+1], pos)))
        self.sequences = seqs
        self.seq = seqs


    def get_potential_module_positions_from_fuzzy_match(self):
        potential_positions = [[] for x in self.node.positions]
        positions = self.node.positions
        for strand in range(len(positions)):
            if self.fuzzy_form[strand]==0:
                potential_positions[strand].append(positions[strand])
            elif self.fuzzy_form[strand]==-1:
                potential_positions[strand].append(positions[strand][:-1])
                potential_positions[strand].append(positions[strand][1:])
            elif self.fuzzy_form[strand]==1:
                if positions[strand][0]>0:
                    potential_positions[strand].append([positions[strand][0]-1] + positions[strand])
                if positions[strand][-1]<(len(self.tree.seq)-1):
                    potential_positions[strand].append(positions[strand] + [positions[strand][-1]+1])

        potential_positions_combinations = list(itertools.product(*potential_positions))
        self.possible_positions = potential_positions_combinations
        self.pos = self.possible_positions
        #print("POSITIONS",self.pos)

    def eval(self,alignment=False,sequences=[],aln_sequences="",ungapped_positions=[], fold_compound=None):
        """
        Method to calculate score
        """
        #print("info",self.module.stackings,self.module.ID,self.module.n_nodes,self.seq)
        try:
            return self.score
        except:
            if len(self.seq)<(len([x for y in self.possible_positions[0] for x in y])+len(self.stackings)-1):
                self.score = -104
                return self.score
		
            #print("MODULE GAPS BEFORE ROTATE",self.module.gaps_per_strand)
            SSE_module_gaps = self.rotate_list(self.module.gaps_per_strand)
            #print("MODULE GAPS AFTER ROTATE",SSE_module_gaps)
            if alignment==False:
                max_score = -100
                current_best = []
                best_seq = ""
                #enumerate
                #print("SEQUENCES TO SCORE",sequences)
                for seq in range(len(self.sequences)):
                    #print("SCORING SEQUENCE FOR ",self.sequences[seq])
                    tmp_seq,tmp_pos = self.reintegrate_module_gaps_to_sse_sequence(self.sequences[seq], SSE_module_gaps, self.possible_positions[seq])
                    #print("SEQ AFTER REMOVING GAPS", tmp_seq)
                    if len(tmp_seq)<(len([x for y in tmp_pos for x in y])+len(self.stackings)-1):
                        score = -104
                    else:
                        score = self.module.eval(tmp_seq)[0]
                        #print("COMPUTING FUZZY SEQUENCE SCORE",score, tmp_pos)

                        #print("MODULE SCORED",tmp_seq,tmp_pos,score)
                    if score>max_score:
                        current_best= tmp_pos
                        best_seq = tmp_seq
                        max_score = score
                    #consider case where two max scores
                    self.pos = current_best
                    self.seq = best_seq
                    if max_score<-50:
                        self.score=-103
                        return max_score
                    if len(self.stackings)>1:
                        if current_best[0][0]<0 or current_best[-1][-1]>=len(self.tree.seq):
                            self.score =  -100
                            return -102
                    else:
                        #in the case of the tuple of size 1
                        #print("CURRENT_BEST",current_best)
                        if current_best[0][0]<0 or current_best[0][-1]>=len(self.tree.seq):
                            self.score =  -100
                            return -101
                    self.dot_bracket = self.set_dot_bracket_ss()

                    SSE_module_gaps = self.rotate_list(self.module.gaps_per_strand)
                    if current_best in self.module.folded_positions:
                        self.score=-111
                        return self.score
                    if max_score<SEQ_THRESH:
                        #so here we don't want to fold but we should still have some form of penalty... trying a flat penalty.
                        self.score = max_score - 3
                    else:
                        self.score = self.module.eval_constraint_folding(max_score,self.tree.seq, current_best, self.dot_bracket, SSE_module_gaps, fold_compound)
                        #print("COMPUTING FUZZY SEQUENCE SCORE WITH CONSTRAINT FOLDING",self.score, tmp_pos)
                        self.module.folded_positions.append(current_best)
                   # print("DONE: ", self.pos, self.seq, self.score)
                    return self.score
            else:

                max_score = -120
                current_best = []
                for seq_ind in range(len(self.sequences)):
                    gapped_cand = 0
                    total_score = 0
                    for full_sequence in sequences:
                        this_seq_pos = self.possible_positions[seq_ind]

                        #print("SEQ BEFORE REMOVING GAPS", full_sequence)
                        full_sequence = [full_sequence[ungapped_pos] for ungapped_pos in ungapped_positions]
                        #print("SEQ AFTER REMOVING GAPS", full_sequence)
                        joined_pos = [x for y in this_seq_pos for x in y]
                        this_seq = "".join([full_sequence[x] for x in joined_pos])
                        #print("SCORING SEQUENCE FOR ",self.sequences[seq])
                        tmp_seq,tmp_pos = self.reintegrate_module_gaps_to_sse_sequence(this_seq, SSE_module_gaps, this_seq_pos)
                        if "-" in tmp_seq:
                            gapped_cand+=1
                            #print("Candidate has a gap, no scoring", this_seq)
                            continue
                        else:
                            if len(tmp_seq)<(len([x for y in tmp_pos for x in y])+len(self.stackings)-1):
                                continue
                            else:
                                score = self.module.eval(tmp_seq)[0]
                                #print("scored candidate in aln:",score, this_seq)
                                total_score+=score

                    mean_score = total_score / (len(sequences)-gapped_cand)
                    if mean_score > max_score:
                        current_best= tmp_pos
                        max_score = mean_score
                        #consider case where two max scores
                        self.pos = current_best
                        self.seq = tmp_seq
                        if max_score<-50:
                            self.score=-103
                            return max_score
                        if len(self.stackings)>1:
                            if current_best[0][0]<0 or current_best[-1][-1]>=len(self.tree.seq):
                                self.score =  -100
                                return -102
                        else:
                            #in the case of the tuple of size 1
                            #print("CURRENT_BEST",current_best)
                            if current_best[0][0]<0 or current_best[0][-1]>=len(self.tree.seq):
                                self.score =  -100
                                return -101


                        SSE_module_gaps = self.rotate_list(self.module.gaps_per_strand)
                        if current_best in self.module.folded_positions:
                            self.score=-111
                            return self.score


                self.dot_bracket = self.set_dot_bracket_ss()
                if max_score<SEQ_THRESH:
                    self.score = max_score
                else:
                    self.score = self.module.eval_aln_constraint_folding(max_score, current_best,
                                                                     self.dot_bracket, SSE_module_gaps, self.tree.seq, fold_compound)
				
                return self.score




    def set_dot_bracket_ss(self):
        #print("GENERATING DOT BRACKET GAPPED STRUCTURE")
        positions = self.pos
        stackings = self.stackings
        string_ss = ""
        ss = []
        if len(positions) == 1:
            strand = positions[0]
            this_strand_ss = ["."] * len(strand)
            this_strand_ss[0]= "("
            this_strand_ss[-1]= ")"
            #print("SINGLE STRAND SS", this_strand_ss)
            string_ss = "".join(this_strand_ss)
            return string_ss
        else:
            for ind,strand in enumerate(positions):
                this_strand_ss = ["."] * len(strand)
                if ind == 0:
                    this_strand_ss[0]= "("
                    this_strand_ss[-1]="("
                elif ind == len(positions)-1:
                    this_strand_ss[0]= ")"
                    this_strand_ss[-1]=")"
                else:
                    this_strand_ss[0]= ")"
                    this_strand_ss[-1]="("
                ss.append("".join(this_strand_ss))
        if sum(stackings)>0:
            for ind,stacking in  enumerate(stackings):
                if stacking>0:
                    if ind==0:
                        for i in range(stacking):
                            ss[0]  = "(" + ss[0]
                            ss[-1] = ss[-1] + ")"
                    else:
                        for i in range(stacking):
                            ss[ind-1]= ss[ind-1] + "("
                            ss[ind] = ")" + ss[ind]
        string_ss = "*".join(ss)
        #print("DOT BRACKET FOR POSITIONS",positions,string_ss)
        return string_ss


    def reintegrate_module_gaps_to_sse_sequence(self, seq, gaps, positions):
        #print("ADDING GAPS TO SEQUENCE",seq,gaps)
        new_seq = ""
        if "&" not in seq:
            for ind, nuc in enumerate(seq):
                if ind not in gaps[0]:
                    new_seq += nuc
        else:
            for strand,subseq in enumerate(seq.split("&")):
                for ind, nuc in enumerate(subseq):
                    if ind not in gaps[strand]:
                        new_seq += nuc
                if strand<len(seq.split("&"))-1:
                    new_seq += "&"

        new_pos = []
        for strand, pos_in_strand in enumerate(positions):
            new_strand = []
            for ind, pos in enumerate(pos_in_strand):
                if ind not in gaps[strand]:
                    new_strand.append(pos)
            new_pos.append(new_strand)
        return new_seq,new_pos
    #def get_position(self, d, detail=True):
    #    try:
    #        return self.pos
    #    except:
    #        l_node = self.node.get_strands_len()[0]
    #        l_module = self.module.get_strands_len()[0]
    #        pos = self.node.get_positions()
    #        if l_module == l_node+1:
    #            if d == 0:
    #                pos = [pos[0]-1]+pos
    #            else:
    #                pos = pos.append(pos[-1]+1)
    #        elif l_module == l_node-1:
    #            if d == 0:
    #                pos = pos[0:-1]
    #            else:
    #                pos = pos[1:]
    #        self.pos = pos
    #        return pos

    def get_request_info(self):
        return (self.tree.seq,self.tree.structure)

ScanResult.register(ExactScanResult)
ScanResult.register(FuzzyScanResult)
