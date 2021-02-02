# -*- coding: utf-8 -*-
import subprocess
import numpy as np
import os
from subprocess import run, PIPE
def map_position_to_nodes(nodelist,positions):
    ordered_nodes = []
    sortednodes = sorted(nodelist)
    sortedpositions=positions.copy()
    sortedpositions=positions
    for i in positions:
        corr_pos = sortedpositions.index(i)
        corr_node= sortednodes[corr_pos]
        ordered_nodes.append(corr_node)
    return ordered_nodes


import copy
def get_constraints_from_BN(graph, positions):
    sortedpos = copy.deepcopy(positions)
    sortedpos = sorted(positions)
    real_positions = []
    #print(graph.nodes(data=True))
    #print(graph.edges(data=True))
    #positions = [x for y in positions for x in y]
    constraints = []
    bps = []
    for i in list(graph.edges(data=True)):
        if (i[2]['label'].upper()=="CWW") and i[0]<i[1]:
            bps.append(i)
    nodes = []
    for i in graph.nodes():
        nodes.append(int(i))
    sortednodes = map_position_to_nodes(nodes,positions)
    #sortednodes = nodes
    matched = False
    visited = []
    #print("sortednodes",sortednodes)
    #print("positions",positions)
    #print("sorted positions", sortedpos)
    #print("bps",bps)
    for j in range(len(sortednodes)):
        #print("currently doing",j)
        #if j>0:
        #    dist = sortednodes[j]-sortednodes[j-1]
        #    if 1<dist<5:
        #        for dd in range(dist):
        #            constraints.append("x")
        #            real_positions.append(real_positions[-1]+dd+1)
                
        n = sortednodes[j]
        matched = False

        for bp in bps:
            (a,b) = bp[0],bp[1]
            #print("BP",bp[0],bp[1],"IN NODES INDEX", nodes.index(bp[0]),nodes.index(bp[1]),"IN POSITIONS", positions[nodes.index(bp[0])],positions[nodes.index(bp[1])])
            if (n==a and b in visited) or (n==b and a in visited):
                
                if positions[nodes.index(bp[0])]>positions[nodes.index(bp[1])]:
                    #print(bp,"we are reversed, we completed bp, adding (")
                    constraints.append("(")
                else:                    
                    #print(bp,"we are straight, we completed bp, adding )")
                    constraints.append(')')
                matched = True
                break
            elif (n==a and b not in visited) or (n==b and a not in visited):
                if positions[nodes.index(bp[0])]>positions[nodes.index(bp[1])]:
                    #print(bp,"we are reversed, starting bp, adding )")
                    constraints.append(")")
                else: 
                    #print(bp,"we are straight, starting bp, adding (")
                    constraints.append('(')
                #constraints.append('(')
                matched = True
                break
        if matched == False:
            constraints.append('x')
        visited.append(n)
        real_positions.append(sortednodes[j])
    return constraints


def call_rnasubopt(seq):
    #out = run(['RNAsubopt', '-e 5', '--noLP'], input=seq, stdout=PIPE, universal_newlines=True)
    #lines = out.stdout.splitlines()
    #print(lines)
    output = []
    #for line in lines:
    #    output.append(line.split(" ")[0])
    lines = []
    while len(lines)<100:
        out = run(['RNAsubopt', '--stochBT_en=5000', '--nonRedundant'], input=seq, stdout=PIPE, universal_newlines=True)
        #out = run(['RNAsubopt', '−−stochBT_en 100000','--nonRedundant'], input=seq, stdout=PIPE, universal_newlines=True,shell=True)
        #out = run(['RNAsubopt', '−−stochBT_en 100000','--nonRedundant'], input=seq, stdout=PIPE, universal_newlines=True,shell=True)
        lines = out.stdout.splitlines()
    #print("SUBOPT OUTPUT",lines)
    for line in lines[1:]:
        #print(line)
        line = [x for x in (line.split(" ")) if x]
        try:


            output.append((line[0],float(line[2])))
        except:
            print("rnasubopt error")
            print(line)
    return output

def call_alifold(file):
    #for line in lines:
    output = []
    print("RUNNING ALIFOLD ON FILE",file,os.getcwd())
    f=open(file,"r")
    out = run(['RNAalifold', "--input-format=S",'−−stochBT_en=5000'], stdin=f, stdout=PIPE, universal_newlines=True, )
    f.close()
    lines = out.stdout.splitlines()
    for line in lines[1:]:
        #print(line.split(" "))
        #print(line)
        try:
            #print((line.split(" ")[0], float(line.split(" ")[1])))
            output.append((line.split(" ")[0], float(line.split(" ")[1])))
        except:
            print("rnaalifold error")
    return output


def call_rnafold(seq,cons="",just_ss=False):
   
    E = 0
    p = 0
    Z = 0

    if cons=="":
        out = run(['RNAfold', '-p', '--noPS'], input=seq, stdout=PIPE, universal_newlines=True)
    else:
        out = run(['RNAfold', '-p', '--noPS', '-C', '--enforceConstraint'], input='{}\n{}'.format(seq, cons), stdout=PIPE, universal_newlines=True)



    lines = out.stdout.splitlines()
    ss = lines[1].split(" ")[0]
    if just_ss:
        return ss
    energy = (lines[1].split(" ")[-1]).strip()
    eval = energy[1:]
    final_energy = eval[:-1]
    E = float(final_energy)
    weight = (lines[-1].split(" ")[7]).strip()
    p = float(weight[:-1])

    #formula :   Z = 1/p * exp(E/ RT))
    T = 274.5
    R =  0.00198717
    first_half = (1/p)
    inside_par = (E/(R*T))
    second_half = np.exp(inside_par)
    Z = first_half * second_half
    #print(Z)  
    return Z


def call_rnaalifold(fname, cons="", just_ss=False):
    E = 0
    p = 0
    Z = 0
    f = open(fname,"r")

    if cons == "":
        #print("RUNNING",'RNAalifold', "--input-format","S", '-p', '--noPS')
        out = run(['RNAalifold', "--input-format","S", '-p', '--noPS'],stdin=f, stdout=PIPE, universal_newlines=True)
    else:
        #print("RUNNING",'RNAalifold', "--input-format","S", '-p', '--noPS', '-C',fname )

        out = run(['RNAalifold',  '-p', '--noPS', '-C',fname], input=cons, stdout=PIPE, universal_newlines=True)
        #out = run(['RNAalifold', "--input-format=S", '-p', '--noPS', '-C=fname'], input='{}\n{}'.format(f, cons), stdout=PIPE, universal_newlines=True)


    lines = out.stdout.splitlines()
   # print("RNAalifold STDOUT lines")
    #for i in range(len(lines)):
    #    print(i,lines[i])
   # print("2ND LINE OF ALIFOLD OUT",lines[2])
    #print("LAST LINE OF ALIFOLD OUT",lines[-1])

    ss = lines[1].split(" ")[0]
    if just_ss:
        return ss
    energy = lines[2].split("=")[0].split("(")[-1].strip()
    #print("MFE ENERGY", energy)

    E = float(energy)
    weight = lines[-1].split(";")[0].split(" ")[-1].strip()
    #print("MFE WEIGHT", weight)

    p = float(weight[:-1])

    # formula :   Z = 1/p * exp(E/ RT))
    T = 274.5
    R = 0.00198717
    first_half = (1 / p)
    inside_par = (E / (R * T))
    second_half = np.exp(inside_par)
    Z = first_half * second_half
    # print(Z)
    return Z


if __name__ == "__main__":
    a = call_rnasubopt("CAACGGCGGGGGUAACUAUGACCCUCUUAAGGUAGCGUAGUACCUUGCCGCAUCAGUAGCGGCUUGCAUGAAUGGAUUAACCAGAGCUUCACUGUCCCAACGUUGGGCCCGGUGAACUGU")
    print(a)
