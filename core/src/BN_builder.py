__author__ = 'Roman'
import itertools
import matplotlib.pyplot as plt
import networkx as nx
import pgmpy
from pgmpy.models import BayesianModel
from pgmpy.models import RNAModule
from module_td import tree_decomposition


def get_cons_seq_from_column(col, aln):
    """
    Takes as input the full alignment information, and the column of interest
    """
    #aln, addresses, seqs = aln_dict
    all_nucs = []
    # print(aln.keys())
    tot = len(aln[col].values())
    for i in aln[col].values():
        all_nucs.append(i)
    # print(all_nucs)
    scores_dict = {}
    for i in "ACGU":
        z = all_nucs.count(i)
        but = all_nucs.count('-')
        scores_dict[i] = float(z + 0.1) / float(tot + 0.4 - but)
    return scores_dict


def get_key_interaction_nodes(cols, n, type, addresses):
    """
    Get the important nodes, that have a frequency over some threshold, that will be used to build the bayes net,
    and outputs them in the format of a list
    """
    nodes = []
    if type == '3d':
        threshold = 0.15
    else:
        threshold = 0.35
    for i in range(len(cols)):
        if cols[i][1] != None:
            if (int(cols[i][1]) / float(n)) > threshold:
                nodes.append(addresses[i])
    return nodes


def build_graph(cols_ss, cols_3d, nodesSS, nodes3d, details, aln_dict, motif_zone):
    """
    Builds a graph with all the nucleotides and nodes that are close to important 3D nodes

    :param cols_ss: all most common secondary structure interactions
    :param cols_3d: all most common 3D interactions
    :param nodesSS:  Important secondary structure nodes
    :param nodes3d:  Important 3D structure nodes
    :param details:  leontis westhof details for the 3D interactions
    :param aln_dict: all alignment information
    :param motif_zone: exact alignment positions of the motifs if given (nodes)
    :return:
    """
    motif = BayesianModel()
    edge_3d = []
    edge_ss = []
    allnodes = []
    keynodes = []
    nearkeynodes = []
    backbone = []

    # nodes that have a significant 3D interaction inside the motif region are keynodes
    for i in nodes3d:
        for j in cols_3d:
            if j[0] == i and int(i) in motif_zone:
                if j[0] == i and int(i) not in keynodes:
                    keynodes.append(int(i))
                if j[0] == i and i not in allnodes:
                    allnodes.append(int(i))
    for zz in motif_zone:
        if int(zz) not in keynodes:
            keynodes.append(int(zz))
    for j in keynodes:
        if j not in nearkeynodes:
            nearkeynodes.append(int(i))

    # add secondary intucture interactions only if they are between two nodes in the immediate vicinity of a keynode
    for i in range(len(cols_ss)):
        if cols_ss[i][0] != None:
            #print(int(motif_zone[i]), cols_ss[i][0])
            if (int(motif_zone[i]) in nearkeynodes and int(cols_ss[i][0]) in nearkeynodes) and (
            int(motif_zone[i]), int(cols_ss[i][0])) not in edge_ss and (
            int(cols_ss[i][0]), int(motif_zone[i])) not in edge_ss:
                if i < cols_ss[i][0]:
                    edge_ss.append((int(motif_zone[i]), int(int(cols_ss[i][0]))))
                else:
                    edge_ss.append((int(int(cols_ss[i][0])), int(motif_zone[i])))

    for i in keynodes:
        allnodes.append(i)
    for i in nearkeynodes:
        if i not in allnodes:
            allnodes.append(i)

    finalnodes = sorted(allnodes)
    weightss = []
    for i in finalnodes:
        weightss.append(1)
    motif.add_nodes_from(finalnodes, weightss)

    # add secondary intucture nodes
    motif.add_edges_from(edge_ss)
    new_edgess = []
    for i in edge_ss:
        new_edgess.append((int(i[0]), int(i[1])))

    # add 3D interaction nodes
    for i in range(len(cols_3d)):
        if cols_3d[i][0] != None:
            if (int(motif_zone[i]) in keynodes and int(cols_3d[i][0]) in keynodes) and motif_zone[i] != cols_3d[i][
                0] and (int(motif_zone[i]), int(cols_3d[i][0])) not in edge_3d and (
            int(cols_3d[i][0]), int(motif_zone[i])) not in edge_3d:
                if i < cols_3d[i][0]:
                    edge_3d.append((int(motif_zone[i]), int(cols_3d[i][0])))
                else:
                    edge_3d.append((int(cols_3d[i][0]), int(motif_zone[i])))

    for i in edge_3d:
        i = (int(i[0]), int(i[1]))
    motif.add_edges_from(edge_3d)

    othernodes = []
    for node in motif.nodes():
        if node not in keynodes:
            othernodes.append(node)

    # builds a secondary graph that copies the bayes net, but with added backbone.
    plotter = nx.Graph()
    plotter.add_nodes_from(motif.nodes())
    plotter.add_edges_from(motif.edges())
    sorted_nodes = (sorted(plotter.nodes()))
    for i in range(1, len(sorted_nodes)):
        plotter.add_edge(sorted_nodes[i - 1], sorted_nodes[i])
        backbone.append((sorted_nodes[i - 1], sorted_nodes[i]))
    pos = nx.spring_layout(plotter)
    nx.draw_networkx_nodes(plotter, pos, node_color='g', node_size=400)
    nx.draw_networkx_nodes(plotter, pos, nodelist=othernodes, node_color='r', node_size=400)
    nx.draw_networkx_edges(plotter, pos, edgelist=edge_ss, width=1.6, edge_color='r')
    nx.draw_networkx_edges(plotter, pos, edgelist=edge_3d, width=2.0, edge_color='g')
    nx.draw_networkx_edges(plotter, pos, edgelist=backbone, width=0.2, edge_color='b')

    labels = {}
    for i in motif.nodes():
        labels[i] = int(int(i) % (motif_zone[0]) + 1)
    nx.draw_networkx_labels(motif, pos, labels, font_size=12)
    plt.show()
    return motif


def build_pwm(positions, aln_dict):
    """
    This function gets the nucleotide statistics for each column of the alignment and returns them as a dictionary
    :param positions: the positions of 
    the alignment, input as a list
    :param aln_dict: all alignment information, needed to call get_cons_seq_from_column
    :return: dictionary of the nuc statistics for each column.
    """
    #print("building position weight matrix")
    scores = {}
    #regex = r""
    for i in positions:
        scores[i] = get_cons_seq_from_column(int(i), aln_dict)
    return scores


def get_immediate_parents(g, node):
    """
    helper function that returns a list of the parents of a node
    :param g: graph of the motif
    :param node: the node we want parents for
    :return: a list of parents
    """
    #parent_dict = tree_decomposition(g)
    #print("TREE DECOMPOSED",parent_dict)
    #parents = parent_dict[node]
    parents = []
    for i in g.edges():
        if i[1] == node:
            parents.append(i[0])
    if node in parents:
        parents.remove(node)
    return parents


def get_partial_prob_data(aln_dict, node, partners, nuc, more_than_one=False):
    #print("GETTING PARTIAL PROB OF",node)
    """
    Extract data for conservatory mutations.
    :param aln_dict: all nucleotide information
    :param node: node of interest
    :param partners: interacting node of interest
    :param nuc: nucleotide of the partner column
    :return: a list of scores, but only including sequences including the nucleotide nuc
    at the partner position
    """
    #print(node,partners,nuc)
    if more_than_one:
        oneD_partners=[]
        for i in partners:
            for j in i:
                oneD_partners.append(int(j))
    else:
        oneD_partners=partners
    aln = aln_dict
    #print("ALL PARTNERS IN ONE LIST:",oneD_partners)
    if len(oneD_partners) == 1:
        partner = int(partners[0])

        example_key = list(aln.keys())[0]
        nucs = []
        scores_dict = {}
        for i in aln[example_key].keys():
            if aln[partner][i] == str(nuc):
                nucs.append(aln[int(node)][i])
        for i in "ACGU":
            z = nucs.count(i)
            gaps = nucs.count('-')
            scores_dict[i] = float(z + 1) / float(len(nucs) - gaps + 4)
    else:
        #print("DONE")
        partners=oneD_partners
        #print("partners",partners)
        nucs = []
        scores_dict = {}
        example_key = list(aln.keys())[0]

        for i in aln[example_key].keys():
            match = True
            for ind,j in enumerate(partners):
                #print("ALIGNMENT",aln)
                #print("nuc",nuc,ind)
                if aln[j][i] != str(nuc[ind]):
                    match = False
            if match == True:
                nucs.append(aln[int(node)][i])
        #print("NUCLEOTIDES OBSERVED",nucs)
        for i in "ACGU":
            z = nucs.count(i)
            gaps = nucs.count('-')
            scores_dict[i] = float(z + 0.1) / float(len(nucs) - gaps + 0.4)

    return scores_dict

def flatten_parents_and(g,node,parents):
    if int(node) in parents:
        parents.remove(node)
    #print("FLATTENING PARENTS of",node)
    clustered_parents= []
    for p in parents:
        anc = list(g._get_ancestors_of(int(p)))
        anc.remove(int(p))
        #print("ANCESTRY OF PARENT:",p,"IS :",anc)
        #anc.remove(node)
        if len(anc)==0:
            #print("no anc")
            clustered_parents.append([int(p)])
            #print(clustered_parents)
        elif len(anc)==1:
           # print("one anc")
            clustered_parents.append([int(anc[0]),int(p)])
           # print(clustered_parents)
        else:
            #print("mutltiple anc")
           # print("ancestors",anc)
            #imm_par= get_immediate_parents(g, p)
            #print("parents of current parent",imm_par)

            #exit(0)
            anc.append(int(p))
            clustered_parents.append(sorted(anc))
            #print(clustered_parents)
    #print("INITIAL:",parents,"CLUSTERED:",clustered_parents)
    return clustered_parents

def parents_in_single_list(flattened_parents):
    oneD_partners=[]
    for i in flattened_parents:
        for j in i:
            oneD_partners.append(int(j))
    return oneD_partners

def get_priors(g, node, seqfreq, priors, aln_dict):
    """
    Get priors for this node  (building the bayes net)
    :param g: graph of the motif (bayes net)
    :param node: node of interest
    :param seqfreq: nucleotide frequency dictionary for each position of the motif
    :param priors: dictionary of all priors (updated recursively)
    :param aln_dict: all alignment information
    :return: updated prior dictionary
    """
    #print("CURRENTLY BUILDING BAYES NET AT NODE :", node)
    if node not in priors.keys():
        parents = get_immediate_parents(g,node)
        #print('NODE',node,"parents",parents, flush=True)

        # if the node has no parents, then it's easy, the conditional probability distribution is a PWM
        if len(parents) == 0:
            values = []
            priors[node] = [[seqfreq[node]['A']], [seqfreq[node]['C']], [seqfreq[node]['G']], [seqfreq[node]['U']]]
        else:
            for j in parents:
                if j not in priors.keys():
                    priors = get_priors(g, j, seqfreq, priors, aln_dict)
        cpd_list = []

        # If the node is not yet in the prior dict, and only has one parent, the size of the CPD is 16 :
        if node not in priors.keys() and len(parents) == 1:
            for parent in parents:
                x = priors[parent]
                for y in range(4):
                    nuclist = ['A', 'C', 'G', 'U']
                    exp_score = get_partial_prob_data(aln_dict, node, [int(parent)], nuclist[y % 4])
                    cpd_list.append([float(exp_score['A']) * float(x[y][0])])
                    cpd_list.append([float(exp_score['C']) * float(x[y][0])])
                    cpd_list.append([float(exp_score['G']) * float(x[y][0])])
                    cpd_list.append([float(exp_score['U']) * float(x[y][0])])
            priors[node] = cpd_list

        # Here if more than one parent, we solve recursively
        elif node not in priors.keys() and len(parents)> 1:
            #print("NODE HAS MORE THAN 1 PARENT")
            cpd_list = []
            cpd_list_with_node_data = []
            prob_parents = []
            for i in parents:
                temp = [0,1,2,3]

                prob_parents.append(temp)
           # print("EXAMPLE ARRAY PROB PARENTS", flush=True)
            #print("ARRAY SHAPE FOR PARENTS PROB",prob_parents)
            tuple_iter = itertools.product(*prob_parents)
            #print("PARENTS PROB TUPLES", prob_parents)

            # generate the list of tuples necessary to fit the structure of the pgmpy CPD format
            # it is the number of combinations between all ancestors = size of the CPD
            tuples = [item for item in tuple_iter]
            #print(tuples)
            # We generated the format of the list, now we fill it
            #all_parents = flatten_parents_and(g,int(node),parents)

            #print("SIZE OF TUPLE SPACE",len(tuples), flush=True)
            for j in tuples:

                combined_prob = 1
                combined_nuc = ""
                nuclist = ['A', 'C', 'G', 'U']
                # print(j)
                for i in range(len(j)):
                    combined_prob = combined_prob * priors[int(parents[i])][j[i]][0]
                    combined_nuc = combined_nuc + nuclist[j[i] % 4]
                #if (tuples.index(j))%10000==0:
                    #print("step",tuples.index(j),"of",len(tuples), flush=True)
                    #print("NUMBER OF COMBINED NUCLEOTIDES:",len(combined_nuc))
                    #print(combined_nuc, combined_prob)
                exp_score = get_partial_prob_data(aln_dict, node, [int(x) for x in parents], combined_nuc)
                #print("current_probs:",exp_score)
                cpd_list_with_node_data.append([float(exp_score['A'])])
                cpd_list_with_node_data.append([float(exp_score['C'])])
                cpd_list_with_node_data.append([float(exp_score['G'])])
                cpd_list_with_node_data.append([float(exp_score['U'])])
                #print("CARDINALITY OF PRIOR ",int(node),int(len(cpd_list_with_node_data)))
            priors[node] = cpd_list_with_node_data
    return priors


def add_cpd(g, i, priors, cpds):
    """
    Once we have all the prior information, we input the cpd in the pygmpy BayesianModel object
    :param g: Bayes Net
    :param i: node of interest
    :param priors: All prior information in terms of parents, nucleotide statistics and probabilities
    :param cpds: list of node CPDs, updated recursively
    :return: updated graph, updated list of CPD
    """
    #print("============================== ADDING CPD",i,"PARENTS",get_immediate_parents(g,i),"====================================")
    #print("===CPDS SO FAR===")
    #print(cpds)
    #print("ADDING CPD FOR NODE",i,"WITH PRIORS",priors[i])
    #print("ADDING CPD OF LENGTH",len(priors[i]), flush=True)
    #print("parents added",cpds)

    if i in cpds:
        #print("nothing to do here")
        return g, cpds
    
    
    parents_card = []
    pars = get_immediate_parents(g, i)
    if len(pars) == 0:
        #print("CREATING CPD",i,priors[i],pars,parents_card,"CPDS SO FAR",cpds)
        #cpd = pgmpy.factors.discrete.CPD.TabularCPD(i, len(priors[i]), priors[i])
        cpd = pgmpy.factors.discrete.CPD.TabularCPD(i, len(priors[i]), priors[i])
        #print("ADDING CPD",i,type(cpd),cpd.variable,cpd.variable_card,cpd.cardinality,cpd.values,"CPDS SO FAR",cpds)
        g.add_cpds(i,cpd)
       
        cpds.append(i)
        #print("CPD SUCCESSFULLY ADDED",cpds)
    else:
        for j in pars:
            if j not in cpds:
                #print("recursively calling add_cpd",j)
                g, cpds = add_cpd(g, j, priors, cpds)
        for j in pars:
            parents_card.append(4)

        # add the CPD to the model
        #print("CREATING CPD",i,priors[i],pars,parents_card,"CPDS SO FAR",cpds)
        #this_cpd = pgmpy.factors.discrete.CPD.TabularCPD(i, 4, priors[i], pars, parents_card)
        this_cpd = pgmpy.factors.discrete.CPD.TabularCPD(i, 4, priors[i], pars, parents_card)
        #print("CPD CREATED", this_cpd.values)
        #print()
        cpd=this_cpd
        #print("ADDING CPD",i,type(cpd),cpd.variable,cpd.variable_card,cpd.cardinality,cpd.values,"CPDS SO FAR",cpds)
        if i not in cpds and i not in g.cpd_order:
            #try:
            g.add_cpds(i,this_cpd)
            cpds.append(i)
            #print("CPD SUCCESSFULLY ADDED",cpds)
            #except:
            #    print("attempted to add int as cpd")
            #    return g, cpds
    return g, cpds


def adjust_order(probs):
    """
    This is a bit technical but the CPDs have to be ordered in a very specific way, so we have to transpose the
    matrix from the previous functions
    :param probs: the CPDs from the graph
    :return: same, but transposed
    """
    new_order = [probs[0]]
    for i in range(1, len(probs)):
        if i % 4 == 0:
            new_order.append(probs[i])
    for i in range(1, len(probs)):
        if i % 4 == 1:
            new_order.append(probs[i])
    for i in range(1, len(probs)):
        if i % 4 == 2:
            new_order.append(probs[i])
    for i in range(1, len(probs)):
        if i % 4 == 3:
            new_order.append(probs[i])
    return (new_order)


def check_sum_1(probs):
    """
    Verifies all CPDs sum up to 1, and correct them if necessary.
    :param probs: the CPDs of the bayes net
    :return: same, but corrected
    """
    i = 0
    return_list = []
    while i < len(probs) - 3:
        templist = [probs[i], probs[i + 1], probs[i + 2], probs[i + 3]]
        sum = 0.0
        for j in templist:
            sum = sum + j[0]
        corrector = []
        newsum = 0
        for t in templist:
            corrector.append(float(t[0]) / sum)
        for t in corrector:
            newsum = newsum + t
        diff = 1.0 - newsum
        templist[0][0] = templist[0][0] + diff
        return_list.append([corrector[0]])
        return_list.append([corrector[1]])
        return_list.append([corrector[2]])
        return_list.append([corrector[3]])
        i = i + 4
    return adjust_order(return_list)

def rebuild_graph(g0, parent_dict):
    #print("REBUILDING GRAPH")
    g = RNAModule()
    g.add_nodes_from(g0)

    for n in parent_dict:
        for parent in parent_dict[n]:
            g.add_edge(parent,n)
    return g

def build_BN(g, seqfreq, aln_dict):
    #print("building BN",flush=True)
    """
    Actually builds the bayes net using the above functions
    :param g: Bayes net object
    :param seqfreq: PWM as a dictionary
    :param aln_dict: all the information about the alignment
    :return: the BayesianModel object
    """
    #print("OLD GRAPH",g.nodes(data=True),g.edges(data=True))
    parent_dict = tree_decomposition(g)
    #print("TREE DECOMP",parent_dict)
    motif = rebuild_graph(g,parent_dict)
    g = motif
    labels = {}
    #pos = nx.circular_layout(motif)
    #nodes = nx.draw_networkx_nodes(motif, pos, nodelist=motif.nodes(), node_color='lightgrey', node_size=500, linewidth=2,alpha=0.99)
    #for index, i in enumerate(list(motif.nodes(data=True))):
    #    labels[i[0]] = i[0]
    #nx.draw_networkx_labels(motif, pos, labels, font_size=10)
    #nx.draw_networkx_edges(motif, pos, edgelist=motif.edges(), edge_color='purple', width=2, arrowsize=25)
    #plt.axis("off")
    #plt.savefig("graphBN.png",format = "png")


    cpds = []
    priors = {}

    # first, get the prior information for each nodei
    #print("ALL NODES:", g.nodes(data=True),g.edges(data=True))
    for i in sorted(list(g.nodes())):
        #print("currently setting priors of node",i, flush=True)
        ancestors = get_immediate_parents(g,i)
        parents = list(ancestors)
        #print("there are",len(parents),"parents",parents,flush=True)
        priors = get_priors(g, i, seqfreq, priors, aln_dict)
        correct_edges =  []
    # verifies the priors are in the correct format
    #print("CHECKING SUM TO 1")
    for i in priors.keys():
        priors[i] = check_sum_1(priors[i])
    #print("ADDING CPDs TO MODEL", flush=True)
    # Adds the CPDs to the model
    for i in sorted(list(g.nodes()),key=float):
        #print("ADDING CPD FOR NODE:",i,priors,cpds, flush=True)
        if i not in cpds:
            g, cpds = add_cpd(g, i, priors, cpds)
    #print("GETTING CPDs")
    #probas = g.get_cpds()
    #for i in g.get_cpds():
        #print(i)
        #print(i.values)

    ##g_copy= nx.DiGraph()
    #g_copy.add_nodes_from(g.nodes())
    #g_copy.add_edges_from(g.edges())
    #bayes_dict = {}
    #for ind,node in enumerate(sorted(list(g.nodes()),key=float)):
    #    bayes_dict[node] = probas[ind].values
    #import pickle
    #pickle.dump(g,open("bayes_net_hairpin.cPickle","wb"))
    return g


