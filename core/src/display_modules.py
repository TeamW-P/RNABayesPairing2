import pickle
import make_logo
import argparse
import networkx as nx
import os
import numpy as np
import matplotlib.pyplot as plt
#from altair import *
import matplotlib
import matplotlib.patches as mpatches
import networkx as nx
import seaborn as sns
#sns.set_style("white")

USE_RFAM=False
CURRENT_DIRECTORY = os.path.dirname(__file__)
GRAPHS_DIRECTORY = os.path.join(CURRENT_DIRECTORY, '../Graphs')

def plot_module(graphs,name,latex=False):
    if latex==True:
        print('using Latex')
        matplotlib.rcParams['text.usetex'] = True

        params= {'text.latex.preamble' : [r'\usepackage{fdsymbol}\usepackage{xspace}']}
        plt.rcParams.update(params)

        for m in range(len(graphs)):
            plt.figure(m)
            labels = {}
            elabels = {}
            elabels = {}
            g0 = graphs[m][0]
            pos = nx.circular_layout(g0)
            nodes = nx.draw_networkx_nodes(g0, pos, nodelist=g0.nodes(), node_color='lightgrey', node_size=500, linewidths=2, alpha=0.99)
            nodes.set_edgecolor("black")
            #print(g0.nodes(data=True))
            #print(g0.edges(data=True))
            for index,i in enumerate(list(g0.nodes())):

                #labels[i] = i
                labels[i] = list(g0.nodes(data=True))[index][1]["nuc"].upper()
            #print(labels)
            nx.draw_networkx_labels(g0,pos,labels,font_size=10,font_weight="bold")

            NC_edges = []
            C_edges = []
            backbone = []
            stacking = []
            for i in g0.edges(data=True):
                if i[0]>i[1]:
                    continue
                i=i[0:2]
                label = g0.get_edge_data(*i)['label'].upper()
                if label not in ["B53","S55","S33"]:
                    if label=="CWW":
                        elabels[i]= r"$\medblackcircle$\xspace"
                    if label=="TWW":
                       elabels[i] = r"$\medcircle$\xspace"
                    if label=="CSS":
                        elabels[i]= r"$\medblacktriangleright$\xspace"
                    if label=="TSS":
                       elabels[i] = r"$\medtriangleright$\xspace"
                    if label=="CHH":
                        elabels[i]= r"$\medblacksquare$\xspace"
                    if label=="THH":
                       elabels[i] = r"$\medsquare$\xspace"


                    if label=="THW":
                       elabels[i] = r"$\medsquare$\xspace $\medcircle$"
                    if label=="CHW":
                       elabels[i] = r"$\medblacksquare$\xspace $\medblackcircle$"
                    if label=="TWH":
                        elabels[i]=r"$\medcircle$\xspace $\medsquare$"
                    if label=="CWH":
                        elabels[i]=r"$\medblackcircle$\xspace $\medblacksquare$"
                    if label=="TSH":
                        elabels[i]=r"$\medtriangleright$\xspace $\medsquare$"
                    if label=="CSH":
                        elabels[i]=r"$\medblacktriangleright$\xspace $\medblacksquare$"
                    if label=="THS":
                        elabels[i]=r"$\medsquare$\xspace $\medtriangleright$"
                    if label=="CHS":
                        elabels[i]=r"$\medblacksquare$\xspace $\medblacktriangleright$"
                    if label=="CWS":
                        elabels[i]=r"$\medblackcircle$\xspace $\medblacktriangleright$"
                    if label=="TWS":
                        elabels[i]=r"$\medcircle$\xspace $\medtriangleright$"
                    if label=="CSW":
                        elabels[i]=r"$\medblacktriangleright$\xspace$\medblackcircle$\xspace"
                    if label=="TSW":
                        elabels[i]=r"$\medtriangleright$\xspace $\medcircle$"

                else:
                    elabels[i]=r""
                if label.upper()=="CWW":
                    C_edges.append(i)
                if label.upper()=="B53":
                    backbone.append(i)
                if label.upper() in ["S55","S33","S35","S53"]:
                    stacking.append(i)
                if label.upper() not in ["S55","S33","S35","S53","B53","CWW"]:
                    NC_edges.append(i)
            #print(NC_edges)
            #print(elabels)
            nx.draw_networkx_edge_labels(g0, pos, elabels, font_size=10)
            nx.draw_networkx_edges(g0, pos, edgelist=NC_edges, edge_color='purple', width=4, arrows=False)
            nx.draw_networkx_edges(g0, pos, edgelist=C_edges, edge_color='green', width=4, arrows=False)
            nx.draw_networkx_edges(g0, pos, edgelist=backbone, edge_color='black', width=3, arrowsize=18)
            nx.draw_networkx_edges(g0, pos, edgelist=stacking, edge_color='orange', width=1, arrows=False)
            plt.title("Module "+str(m),fontsize= 36)
            #print("done")
            plt.axis("off")
            NCP = mpatches.Patch(color="purple", label="Non-canonical" )
            CP = mpatches.Patch(color="green", label="Canonical")
            BP = mpatches.Patch(color="black", label="Backbone" )
            SP = mpatches.Patch(color="red", label="Stacking")
            #plt.legend(handles=[NCP,CP,BP,SP],prop={'size': 16})
            #plt.show()
            plt.savefig(os.path.join(CURRENT_DIRECTORY, "../Graphs/"+name+"_graph"+str(m)+".png"),format="png")
    else:
        print('plotting without latex; for optimal results use latex')
        matplotlib.rcParams['text.usetex'] = False
        for m in range(len(graphs)):
            plt.figure(m)
            labels = {}
            elabels = {}
            g0 = graphs[m][0]
            pos = nx.circular_layout(g0)
            nodes = nx.draw_networkx_nodes(g0, pos, nodelist=g0.nodes(), node_color='lightgrey', node_size=500,
                                           linewidths=2, alpha=0.99)
            nodes.set_edgecolor("black")
            #print(g0.nodes(data=True))
            #print(g0.edges(data=True))
            for index, i in enumerate(list(g0.nodes())):
                # labels[i] = i
                labels[i] = list(g0.nodes(data=True))[index][1]["nuc"].upper()
            #print(labels)
            # nx.draw_networkx_labels(g0,pos,labels,font_size=25,font_weight="bold")
            nx.draw_networkx_labels(g0, pos, labels, font_size=20)

            NC_edges = []
            C_edges = []
            backbone = []
            stacking = []
            for i in g0.edges(data=True):
                if i[0] > i[1]:
                    continue
                i = i[0:2]
                label = g0.get_edge_data(*i)['label'].upper()
                if label not in ["B53", "S55", "S33"]:
                    if label == "CWW":
                        elabels[i] = r"●"
                    if label == "TWW":
                        elabels[i] = r"○"
                    if label == "CSS":
                        elabels[i] = r"▲"
                    if label == "TSS":
                        elabels[i] = r"∆"
                    if label == "CHH":
                        elabels[i] = r"■"
                    if label == "THH":
                        elabels[i] = r"□"

                    if label == "CHW":
                        elabels[i] = r"■●"
                    if label == "THW":
                        elabels[i] = r"□○"
                    if label == "CWH":
                        elabels[i] = r"●■"
                    if label == "TWH":
                        elabels[i] = r"○□"
                    if label == "TSH":
                        elabels[i] = r"∆□"
                    if label == "CSH":
                        elabels[i] = r"▲■"
                    if label == "THS":
                        elabels[i] = r"□∆"
                    if label == "CHS":
                        elabels[i] = r"■▲"
                    if label == "CWS":
                        elabels[i] = r"■●"
                    if label == "TWS":
                        elabels[i] = r"○∆"
                    if label == "CSW":
                        elabels[i] = r"▲●"
                    if label == "TSW":
                        elabels[i] = r"∆○"
                #elif label in ["S35","S53", "S55", "S33"]:
                #    elabels[i]=r"s     "

                else:
                    elabels[i] = r""
                if label.upper() == "CWW":
                    C_edges.append(i)
                if label.upper() == "B53":
                    backbone.append(i)
                if label.upper() in ["S55", "S33", "S35", "S53"]:
                    stacking.append(i)
                if label.upper() not in ["S55", "S33", "S35", "S53", "B53", "CWW"]:
                    NC_edges.append(i)
            #print(NC_edges)
            #print(elabels)
            nx.draw_networkx_edge_labels(g0, pos, elabels, font_size=10)
            nx.draw_networkx_edges(g0, pos, edgelist=NC_edges, edge_color='purple', width=4, arrows=False)
            nx.draw_networkx_edges(g0, pos, edgelist=C_edges, edge_color='lightseagreen', width=4, arrows=False)
            nx.draw_networkx_edges(g0, pos, edgelist=backbone, edge_color='black', width=3, arrowsize=18)
            nx.draw_networkx_edges(g0, pos, edgelist=stacking, edge_color='orange', width=1, arrows=False)
            plt.title("Module " + str(m), fontsize=24)
            #print("done")
            plt.axis("off")
            NCP = mpatches.Patch(color="purple", label="Non-canonical")
            CP = mpatches.Patch(color="green", label="Canonical")
            BP = mpatches.Patch(color="black", label="Backbone")
            SP = mpatches.Patch(color="red", label="Stacking")
            # plt.legend(handles=[NCP,CP,BP,SP],prop={'size': 16})
            #plt.show()
            plt.savefig(os.path.join(GRAPHS_DIRECTORY, name + "_graph" + str(m) + ".png"), format="png")
            plt.clf()

"""
def logo_module(aln,name):
    extra_seqs = {}
    #extra_seqs  = pickle.load(open("../models/"+name + "_rfam.cPickle",'rb'))
    for i,gs in enumerate(aln):
        if not USE_RFAM:
            extra_seqs[i]= []
        #print("MODULE NUMBER",i)
        with(open("../test/seq"+str(i)+".txt", "w")) as f:
            for g in gs[1:]:
                n = sorted(g.nodes(data=True))
                pos = [z[0] for z in n]
                seq ="".join([z[1]['nuc'] for z in n])
                length = len(seq)
                #print(pos)
                #print(seq)
                if len(seq)>0:
                    f.write(seq)
                    if "\n" not in seq:
                        f.write('\n')
            for seq in extra_seqs[i]:
                print(seq)
                if "-" not in seq :
                    if len(seq)>length:
                        seq=seq[len(seq)-length:]
                    f.write("".join(seq))
                    f.write("\n")
        make_logo.make_logo("../test/seq"+str(i)+".txt", "../Graphs/"+name+"_logo"+str(i))
        #os.remove("../test/seq" + str(i) + ".txt")
"""
def logo_module(name):
    seqs = pickle.load(open(os.path.join(CURRENT_DIRECTORY, "../models/"+name + "_sequences.pickle"),'rb'))[0]
    for i,all_seqs in enumerate(seqs):
        with(open(os.path.join(CURRENT_DIRECTORY,"../test/seq"+str(i)+".txt"), "w")) as f:
            for seq in all_seqs:
                if "-" not in seq :
                    f.write(seq)
                    f.write("\n")
    
        make_logo.make_logo(os.path.join(CURRENT_DIRECTORY, "../test/seq"+str(i)+".txt"), os.path.join(GRAPHS_DIRECTORY, name+"_logo"+str(i)))

#deprecated, need to update
def tikz_module(data,name):
    data = [x[0] for x in data]
    print("LEN DATA")
    print(len(data))
    data = ncDraw.remove_already_there(data, GRAPHS_DIRECTORY)
    lb, ub = 0, len(data)
    tot = sum(1 for x in os.listdir(os.path.join(CURRENT_DIRECTORY, '../Graphs')) if (not x.startswith('.') and x.endswith('nxpickle')))

    # lb,ub = 1,5
    for i in range(lb, ub):
        g = data[i]
        i += tot
        #        print "NCM",(i+1)
        nx.write_gpickle(g, GRAPHS_DIRECTORY + '/' + name+'_module-%s.nxpickle' % (i))
        print('preparing for strands')
        strands = ncDraw.buildStrands(g)
        print('ordering strands')
        orderStrands = ncDraw.placeStrands(g, strands)
        print('preparing coords')
        coords = ncDraw.shiftHeight(g, strands, orderStrands)
        print('creating file')
        ncDraw.createTikz(i, g, coords,name)
        # collate(lb, ub)




if __name__ == "__main__":
    arguments = {}
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", help="Dataset name (create new or add to existing)", required=True)
    parser.add_argument("-latex", help="Use latex)")
    args = parser.parse_args()
    model_name = args.n
    if args.latex != None:
        latex=True
    else:
        latex=False
    aln = pickle.load(open(os.path.join(CURRENT_DIRECTORY, "../models/"+model_name+"_one_of_each_graph.cPickle"),"rb"))
    #for i in aln[0]:
    #   print(len(list(i.nodes)))

    if not os.path.exists(GRAPHS_DIRECTORY):
        os.makedirs(GRAPHS_DIRECTORY)

    logo_module(model_name)

    plot_module(aln,model_name,latex=latex)
    print('done')
