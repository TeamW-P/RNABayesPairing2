import pickle
import numpy as np
from matplotlib import pyplot as plt
import statistics
import math
import matplotlib.pyplot as plt
#from altair import *
import matplotlib.patches as mpatches
import networkx as nx
import pandas
import seaborn as sns
from pylab import rcParams
sns.set(style="white")
import pickle
import random
top = []

INITIAL_LAMBDA = 0.35

def correct_lambda_error(mod,scores,size_dist,Lambda=0.35):
    scores2 = []
    for seq in i:
        if len(seq)<3:
            i2.append([])
            continue
        seqscore = seq[3]
        struct_score = (seq[2]-seq[3])/INITIAL_LAMBDA
        scores2.append([seq[0],seq[1],seqscore-abs(Lambda*struct_score)])
    scores2.sort(key=lambda k: (k[0],k[2]),reverse=True)

    return scores2

"""
def correct_lambda_error(mod,i,size_dist,Lambda=0.35):
    i2 = []
    for seq in i:
        #print(seq)
        if len(seq)<3:
            i2.append([])
            continue
        seqscore = seq[3]
        struct_score = (seq[2]-seq[3])/0.35
        #print(struct_score)
        i2.append([seq[0],seq[1],seqscore-abs(Lambda*struct_score)])
        #i2.append([seq[0],seq[1],seqscore-abs(Lambda*struct_score)-math.log(size_dist[mod]),seqscore])
        #print(i2[0][2])
    i2.sort(key=lambda k: (k[0],k[2]),reverse=True)

    return i2
"""
def threshold (N,Beta = 0.2):
    return Beta

def get_precision(TP,FP,TN,FN):
    return (TP)/(TP+FP)

def get_recall(TP,FP,TN,FN):
    return (TP)/(TP+FN)

def MCC(TP,FP,TN,FN):
    return (TP* TN - FP * FN)/math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))

def get_f1_score(TP,FP,TN,FN):
    p = get_precision(TP,FP,TN,FN)
    r = get_recall(TP,FP,TN,FN)
    f1 = 2*p*r/(p+r)
    return f1


graphs = pickle.load(open("../models/bp2june_one_of_each_graph.cPickle", "rb"))
excluded_modules = [194,195,196,197,198,199,200]
n_mods_of_size = {}
for m in range(len(graphs)):
    if m not in excluded_modules:
        size = len(graphs[m][0].nodes())
        if size in n_mods_of_size:
            n_mods_of_size[size] = n_mods_of_size[size]+1
        else:
            n_mods_of_size[size]=1
SIZE_DIST= {}
for m in range(len(graphs)):
    if m not in excluded_modules:
        size = len(graphs[m][0].nodes())
        SIZE_DIST[m] = n_mods_of_size[size]

positives= []
negatives = []
pos_sizes = []
neg_sizes = []

bayespairing_prediction = []
actual_module = []
moip_motifs =  []
to_sort = []


a = pickle.load(open("../src/BP2_JUNE2020_66MODULES.pickle","rb"))
#a = pickle.load(open("../src/2fold_RECOMB.pickle","rb"))
for motif in a:
    if motif in excluded_modules:
        continue
    motif_size = len(graphs[motif][0].nodes)
    motif_scores = []
    
    for sequence in a[motif]:
        if len(sequence)==0:
            continue
       
        #sequence=eval_mod(motif,sequence,SIZE_DIST)
        sequence = correct_lambda_error(motif,sequence,SIZE_DIST,0.5)
        seq_scores = []
        for result in sequence:
            #print(motif,result)
            score = result[2]
            seq_scores.append(score)
            actual_module.append(score)
            if result[0]==1:
                bayespairing_prediction.append(result[2])
                positives.append(result[2])
                pos_sizes.append(motif_size)
            else:
                negatives.append(result[2])
                neg_sizes.append(motif_size)
        if len(seq_scores) == 0:
            negatives.append(-20)
            neg_sizes.append(motif_size)
            seq_scores = [0]
        motif_scores.append(max(max(seq_scores),0.00))
    if len(motif_scores) > 0:
        motif_mean = np.mean(motif_scores)
    else:
        motif_mean = 0
    if motif_mean>-1 and len(motif_scores)>0:
        top.append(motif_mean)
        moip_motifs.append(motif)
        to_sort.append((len(motif_scores), motif_mean, motif))
final_scores1 = sorted(to_sort,key=lambda tup: tup[0], reverse=True)
print("N POSITIVES",len(positives), "N NEGATIVES", len(negatives))
positives = random.sample(positives,1000)
negatives = random.sample(negatives,len(positives))

scores = []
classes = []
for i in positives:
    scores.append(i)
    classes.append("True")
for j in negatives:
    scores.append(j)
    classes.append("False")
df = pandas.DataFrame({"scores":scores, "Ground truth":classes})
df['zero'] = ""
f, ax = plt.subplots(figsize=(7,5.75))
violin = sns.violinplot(x="scores",y="zero",hue='Ground truth', split=True, data=df, palette="Set1", inner="quart", bw =.2, cut=0, linewidth=1.8)
#leg=violin.collections
plt.setp(violin.get_legend().get_texts(),fontsize="18")
plt.setp(violin.get_legend().get_title(),fontsize="18")

sns.despine(left=True)
f.suptitle('BayesPairing2 log odds score on true and false hits', fontsize=19)
ax.set_ylabel("")
ax.set_yticklabels([])
sns.set_context("poster")
ax.set_xlabel("log odds score",size = 19,alpha=0.99)
ax.xaxis.set_ticks_position('none')
violin.tick_params(axis='x', labelsize=14)
plt.savefig("violin_null_final_bigger.pdf",format="pdf")
plt.show()

Betas = np.linspace(-10,10,200)

optimal_mcc = -1
optimal_cutoff = -1111
optimal_line = {}


ROC_values = []
ROC_keys = []
FDR = []
for Beta in Betas:
    true_pos = 0
    false_pos = 0
    true_neg=0
    false_neg =0
    for ind,neg in enumerate(negatives):
        if neg>threshold(neg_sizes[ind],Beta):
            false_pos = false_pos + 1
        if neg<threshold(neg_sizes[ind],Beta):
            true_neg = true_neg +1
    for ind,pos in enumerate(positives):
        if pos>threshold(pos_sizes[ind],Beta):
            true_pos+=1
        if pos<threshold(pos_sizes[ind],Beta):
            false_neg+=1
    sens = true_pos/(true_pos + false_neg)
    fpr = (false_pos/(false_pos+true_neg))
    tFDR= false_pos/(false_pos+true_pos)
    f1 = get_f1_score(true_pos,false_pos,true_neg,false_neg)
    try:
        mcece = MCC(true_pos,false_pos,true_neg,false_neg)
    except:
        mcece=-1
    ROC_values.append(sens)
    ROC_keys.append(fpr)
    FDR.append(tFDR)
    if mcece>optimal_mcc:
        optimal_mcc = mcece
        optimal_cutoff = Beta
        optimal_line = {"cutoff":Beta,"MCC":mcece,"FPR":fpr,"FDR":tFDR,"F1":f1}
        #print(Beta,mcece,tFDR)
        
print("CUTOFF",round(optimal_line["cutoff"],3),"MCC",round(optimal_line["MCC"],3),"FPR",round(optimal_line["FPR"],3), "FDR", round(optimal_line["FDR"],3),"F1",round(optimal_line["F1"],3))

plt.clf()
plt.plot(ROC_keys, ROC_values, color="Blue", label = "ROC curve")
#plt.ylim(0.2,1.1)
#plt.xlim(-0.001,0.14)
plt.xlabel("FPR")
plt.ylabel("TPR")

#plt.axvline(x=optimal_cutoff,color="red", label = "Selected Beta value")
plt.plot([0],[0], color="green", label = "Beta values")

count=0
for i,j in zip(ROC_keys,ROC_values):
    if count%18==0:
        plt.annotate(round(Betas[ROC_keys.index(i)],1),xy=(i,j+0.03), color="green", label = "Cutoff")
    count+=1
plt.legend()
plt.savefig("FDR.pdf",format="pdf")
plt.show()
