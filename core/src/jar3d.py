import pickle
from BayesPairing import jared
from . import make_BN_from_carnaval as makeBN
import statistics
import sys
import ast
from Bio import SeqIO
#mod_seqs,full_seqs,indexes = pickle.load(open("../models/RECOMB_sequences.pickle","rb"))

tr = {"A":0,"C":1,"G":2,"U":3}
def get_weights(seql):
    weights = [0,0,0,0]
    for seq in seql:
        for c in seq:
            weights[tr[c]]+=1
    fweights=[x/sum(weights) for x in weights]
    return fweights


            
from numpy.random import choice
def gen_seq(weights,size):
    seq = ""
    nucs = ["A","C","G","U"]
    for s in range(size):
        nuc = choice(nucs,1,weights)
        #print("chosen:",nuc)
        seq= seq + nuc[0]
    return seq

def direct_score(modules,file,dataset="RECOMB",is_list=False,retrain=False,THRESHOLD=3.5,pretrained=False):
    
    
    if pretrained:
        BNs = pickle.load(open("../models/" + dataset + "_models.pickle", "rb"))
    else:
        BNs = {}
        for module in modules:
            BN = makeBN.call_makeBN(module, dataset,"NONE",False,"", retrain=retrain)
            BNs[module] = BN
        print("TRAINED BAYES NETS",BNs)
    if is_list:
        sequences = file
    else:
        seqs = []
        records = SeqIO.parse(open(file,"r"),"fasta")
        for i in records:
            seqs.append(str(i.seq))
        sequences = seqs
        
    truth_scores = {}
    mod_scores = {}
    #THRESHOLD = 4.372
    #THRESHOLD = -1.99
    positive_scores = []
    for module in modules:
        scores = []
        truths = []
        for ind,seq in enumerate(sequences):
            #print(seq)
            if "N" in seq or "D" in seq or "E" in seq or "F" in seq or "H" in seq:
                continue
            elif "*" in seq:
                seq = seq.replace("*","")
            score = jared(seq,BNs[module],THRESHOLD)
            scores.append(score[1])
            truths.append(score[0])
        mod_scores[module] = scores
        truth_scores[module] = truths.count(True)/len(truths)
    
        
    #print("SHUFFFLED",shuffled_scores)
    #print("POSITIVE",statistics.mean(mod_positive_scores))
    #pickle.dump(positive_scores,open("JAR3D_results","wb"))
    #print("NUMBER OF SEQS",len(sequences))
    return mod_scores,truth_scores





if __name__ == "__main__":
    mod = ast.literal_eval(sys.argv[1]) 
    file = sys.argv[2]
    dataset = sys.argv[3]
    pre = sys.argv[4]
    
    if pre=="1":
        print("pretrained")
        a = direct_score(mod,file,dataset,pretrained=True)
    else:
        a = direct_score(mod,file,dataset,pretrained=False)
    print(a[1])