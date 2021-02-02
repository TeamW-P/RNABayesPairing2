from Bio import SeqIO
import numpy as np
from numpy.random import choice


def get_consensus(module_seqs):
    consensus = []
    for i in range(len(module_seqs[0])):
        scores = [1,1,1,1]
        alph = list("ACGU")
        for j in module_seqs:
            if j[i]=="A":
                scores[0]= scores[0]+1
            if j[i]=="C":
                scores[1]= scores[1]+1
            if j[i]=="G":
                scores[2]= scores[2]+1
            if j[i]=="U":
                scores[3]= scores[3]+1
        #print(scores)
        #print(scores.index(max(scores)))
        #print(np.argmax(scores))
        tot = max(sum(scores),4)
        consensus.append([x/tot for x in scores])
    return consensus

def rfam_to_module(positions=[], family_file="",output_name="module_seqs.fasta", as_list=False, seqs=[]):
    alph = ["A","C","G","U"]


    if as_list:
        module_seqs = []
        for seq in seqs:
            module_seqs.append("".join(seq))


        cons = get_consensus(module_seqs)

        fmodule_seqs = []
        for s in module_seqs:
            this_mod = ""
            for col in range(len(s)):
                if str(s[col]) != "-":
                    this_mod = this_mod + str(s[col])
                else:
                    #print("CONSENSUS:", cons, "CURRENT COLUMN", col)
                    this_mod = this_mod + alph[choice(np.array([0,1,2,3]),p=cons[col])]

            fmodule_seqs.append(this_mod)

        with open(output_name, "w") as f:
            for j in range(len(fmodule_seqs)):
                f.write(">module_seq")
                f.write(str(j))
                f.write("\n")
                f.write(fmodule_seqs[j])
                f.write("\n")
        return

    module_seqs = []
    modules_handles = []
    for record in SeqIO.parse(family_file,"fasta"):
       # print(record.seq)
        modules_handles.append(str(record.id))
        fullseq = list(str(record.seq))
        this_mod = ""
        for col in positions:
            this_mod = this_mod + str(fullseq[col-1])
        module_seqs.append(this_mod)
    cons = get_consensus(module_seqs)
    
    fmodule_seqs = []
    for s in module_seqs:
        this_mod = ""
        for col in range(len(s)):
            if str(s[col]) != "-":
                this_mod = this_mod + str(s[col])
            else:
                this_mod = this_mod + alph[choice(np.array([0,1,2,3]),p=cons[col])]

        fmodule_seqs.append(this_mod)
    #print(module_seqs)


    #print(fmodule_seqs)

    with open(output_name,"w") as f:
        for j in range(len(modules_handles)):
            f.write(">")
            f.write(modules_handles[j])
            f.write("\n")
            f.write(fmodule_seqs[j])
            f.write("\n")



