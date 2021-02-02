import pickle
from Bio import SeqIO
from Bio import AlignIO
import networkx as nx
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align.Applications import ClustalwCommandline

from Bio.Align.Applications import ClustalOmegaCommandline
import os
from bs4 import BeautifulSoup
from bs4.element import Comment
import urllib
from urllib.request import urlopen
import time
import pickle
import requests

rfam = pickle.load(open("../models/RFAM_databse.cPickle", "rb"))


def get_from_rfam(fam):
    url = "http://rfam.org/family/"+fam+"/alignment/stockholm"
    r  = requests.get(url)
    open(fam+'.stockholm.txt', 'wb').write(r.content)


def get_rfam_alignment(family):
    n = "RF"+format(family,"05d")
    if (n+".stockholm.txt") not in os.listdir("."):
        get_from_rfam(n)
    try:
        align = AlignIO.read(open(n+".stockholm.txt"), "stockholm")
        
    except:
        align = AlignIO.read(open(n+".stockholm.txt",encoding="utf-16"), "stockholm")
        
    AlignIO.write(align,open("RF"+format(family,"05d")+".clustal","w"), "clustal")
    return align

def get_PDB_sequence(PDBid):
    #print("PDB")
    print(PDBid)
    PDB, chain = PDBid.split(".")
    # print(PDB)
    # print("../all_graphs_pickled/" + PDB + ".nxpickled")
    try:
        g = pickle.load(open("../Graphs/all_graphs_pickled/" + PDB + ".nxpickled", "rb"))
    except FileNotFoundError:
        print("PDB FILE NOT FOUND")
        return ("", 0)
    seq = ""
    nodes = []
    for node in g.nodes(data=True):
        # print(node)
        # print(node[0][0],chain)
        if node[0][0] == chain:
            nodecode = (node[0][1])
            nodes.append((int(nodecode), node[1]["nt"]))
    sortednodes = sorted(list(nodes))
    nuc_by_node = {}
    missing_nuc = False
    # print("NODES")
    numericals = [x[0] for x in sortednodes]
    decalage = 0
    if 1 not in numericals:
        decalage = decalage + 1
        sortednodes.append((1, "N"))
        sortednodes = sorted(sortednodes)
        # missing_nuc=True
        # decalage = decalage +1
        newnodes = []
        # for i in sortednodes:
        #    newnodes.append((i[0],i[1]))
        # sortednodes = sorted(list(newnodes))
        numericals = [x[0] for x in sortednodes]
        # print("MISSING 1", PDBid)
        # print(numericals)
        #
        # sortednodes=sorted(sortednodes)
        # numericals = [x[0]-1 for x in sortednodes]
        # numericals.insert(0,0)
        # else:
        # print("NOT MISSING", PDBid)
    for i in sortednodes:
        nuc_by_node[i[0]] = i[1]
    # print(sortednodes)
    for i in range(1, int(sortednodes[-1][0]) + 1):
        if i not in numericals:
            "NOT IN NODES"
            seq = seq + "-"
        else:
            seq = seq + nuc_by_node[i]
    ss = g.graph['ss']
    print(seq)
    # print("MISSING_NUC",PDBid,missing_nuc)
    return (seq, 0)

def match_rfam_to_PDB(PDB, family, range="",written=[],newfam=False):
    if newfam==True:
        written=["empty"]
    aln = get_rfam_alignment(family)
    #print("PDB_rfam")
    #print(PDB)
    seq,decalage= get_PDB_sequence(PDB)
    gapped_seq = ""
    #print("FAMILY")
    #print(rfam[family])
    for struct in rfam[family][PDB.split(".")[0]]:
        #print("SEQ IN FAM ", family)
       # print(struct)
       # print(PDB)
        #print(rfam[family][PDB.split(".")[0]])
        if rfam[family][PDB.split(".")[0]][0]==PDB.split(".")[1]:
           # print("VALID")
            positions = rfam[family][PDB.split(".")[0]][1]

    #seq_range= [int(i) for i in positions.split(" - ")]
    #print("SEQ RANGE",seq_range)
    #seq=seq[seq_range[0]-1:seq_range[1]-1]
    #print(seq)

    #print(written)
    if len(seq)>5 and PDB not in written:
        #print("writing")
        with open("seq2.fasta","a") as f:
            f.write(">seq"+str(PDB)+"\n")
            f.write(str(seq)+"\n")
        #written.append(PDB)
    #else:
       # print("NOT WRITTEN")
        #SeqIO.write(seqObj,open("seq.clustal","w"),"clustal")
    #for record in aln:
    #    current_seq = str(record.seq).replace("-","")
    #    print(seq)
    #    print("----------")
    #    print(format_alignment(*(pairwise2.align.globalxx(seq,current_seq))[0]))
    #    if seq in current_seq:
    #        gapped_seq=record.seq
    #        break
    #if len(gapped_seq)==0:
    #    exit("NO MATCH FOUND")
    return gapped_seq,written, 0,0


def get_matches(rfam,rfam_by_PDB, PDBs):
    i =PDBs
    matches = {"None":0}
    positions= {}
    for j in i:
        PDB_name = (j.split(".")[0], j.split(".")[1])
        PDB = ".".join(PDB_name)
        ind=i.index(j)
        if (PDB_name in rfam_by_PDB) and (PDB_name is not None):
            print("PDB :",PDB_name)
            for fam in rfam_by_PDB[PDB_name]:
                if fam in matches:
                    print(matches[fam][1])
                    new_list = list(matches[fam][1])
                    new_list.append(PDB_name)
                    matches[fam]=(matches[fam][0]+1,new_list)
                    #print("AFTER",matches[k][fam])
                    print("test",rfam[fam][PDB_name[0]][0])
                    if rfam[fam][PDB_name[0]][0]==PDB_name[1]:
                        if fam in positions:
                            #print(list(graphs[k][ind].nodes(data=True)),fam,rfam[fam][PDB_name[0]][1])
                            positions[fam].append(rfam[fam][PDB_name[0]][1])
                        else:
                            positions[fam]=[rfam[fam][PDB_name[0]][1]]

                else:
                    matches[fam]=(1,[PDB_name])
        else:
            matches["None"]+=1
    return matches

def gapped_column_from_ungapped_nuc_pos(gapped,ungapped,ungapped_pos):

    # --01-2

    gapped_counter=-1
    ungapped_counter = -1
    for i in range(len(gapped)):
        if gapped[i]=="-":
            gapped_counter=gapped_counter+1
        else:
            gapped_counter=gapped_counter+1
            ungapped_counter = ungapped_counter+1


        if ungapped_counter==ungapped_pos:
            #print("GAPPED",gapped_counter-1)
            return gapped_counter

    return -1000

def pos_from_fred(fred_list):
    return [x[0] for x in fred_list[0]]
def nucs_from_fred(fred_list):
    return [x[1]["nuc"] for x in fred_list[0]]

    #KDE for figures

def get_tot_counts(counts):
    tot = 0
    for i in counts:
        tot = tot + counts[i]
    return tot
def define_acceptable_columns(counts):
    #print("ACCEPTABLE COLUMNS")
    #print(counts)
    acceptable_scores = []
    acceptable_cols = []
    noise = []
    threshold = 0+0.00 * get_tot_counts(counts)
    for i in counts:
        if counts[i]>threshold and i!=-1000:
            acceptable_cols.append(i)
            acceptable_scores.append(counts[i])
    for i in range(len(acceptable_cols)):
        for j in range(i+1,len(acceptable_cols)):
            #print(abs(acceptable_cols[i]-acceptable_cols[j])<10)
            if abs(acceptable_cols[i]-acceptable_cols[j])<10:
                if acceptable_scores[i]> acceptable_scores[j]:
                    noise.append(acceptable_cols[j])
                else:
                    noise.append(acceptable_cols[i])
    noise = set(noise)
    #print("NOISE", noise)
    for i in noise:
        acceptable_cols.remove(i)
    return acceptable_cols

def get_bps_from_ss(ss, columns):
    starts = []
    ends = []
    ss=list(ss)
    for i in range(len(list(ss))):
        if ss[i]=="(" or ss[i]=="<" or ss[i]=="{":
            starts.append(i)
        if ss[i]==")" or ss[i]==">" or ss[i]=="}":
            ends.append(i)
    ends.reverse()
    #struct = [(columns[starts[i]],columns[ends[i]]) for i in range(min(len(starts),len(ends)))]
    if min(len(starts),len(ends))==0:
        return ()
    else:
        struct = [(starts[i],ends[i]) for i in range(min(len(starts),len(ends)))]
        return struct
def verify_columns_ss(columns,family,module, graphs):
    rfam_ss=get_ss_from_rfam_family(family,columns)
    rfam_bps = get_bps_from_ss(rfam_ss,columns)
    graph_ss=get_ss_from_graph(graphs[module][0])
    print("COMPARING :",rfam_bps,graph_ss)
    compatible = True
    for bp in graph_ss:
        if bp not in rfam_bps:
            compatible=False

    if compatible:
        return True
    else:
        return False


def find_good_columns(module_columns_by_family, family):
    print('FINDING GOOD COLUMNS',module_columns_by_family)
    good_cols=[]
    cols = []
    module_length = -1
    for col in module_columns_by_family[family]:
        module_length = len(module_columns_by_family[family][col])
        break
    # print(module_length)
    for node_number in range(module_length):
        good_cols_for_this_node = []
        counts={}
        total = 0
        for example in module_columns_by_family[family]:
            this_col = module_columns_by_family[family][example][node_number]
            #print(this_col)
            if this_col in counts:
                counts[this_col]= counts[this_col]+1
            else:
                counts[this_col]=1
            total = total+1
        print("COUNTS", counts)
        acceptable_cols = define_acceptable_columns(counts)
        for i in counts:
            if i in acceptable_cols and i not in good_cols_for_this_node:
                good_cols_for_this_node.append(i)
        good_cols.append(good_cols_for_this_node)
        #print(counts)
    return good_cols

def check_columns_ss(good_cols,family,module,graphs):
    final_cols = [[] for x in good_cols]
    for i in range(len(good_cols[0])):
        this_cols = [x[i] for x in good_cols]
        print("GETTING SS FOR COLUMNS:", this_cols)
        it_fits = verify_columns_ss(this_cols, family, module,graphs)
        if it_fits==True:
            for col in range(len(final_cols)):
                final_cols[col].append(this_cols[col])
    good_cols = final_cols
    return good_cols
def most_common(lst):
    return max(set(lst), key=lst.count)

def get_ss_from_graph(graph):
    #print(graph.nodes())
    #print(graph.edges(data=True))
    ss = []

    edges = list(graph.edges(data=True))
    for i in edges:
        if i[2]["label"]=="cWW" or i[2]["label"]=="CWW":
            if i[0]<i[1]:
                a = list(graph.nodes()).index(i[0])
                b = list(graph.nodes()).index(i[1])
                ss.append((a,b))
    return ss

def is_complementary(a,b):
    if a=="C":
        if b=="G":
            return True
    if a=="G":
        if b=="C":
            return True
    if a=="A":
        if b=="T":
            return True
    if a=="T":
        if b=="T":
            return True
    return False

def rule_out_impossible_seqs(graph,sequences):
    ss = get_ss_from_graph(graph)
    final_sequences = sequences

    for seq in sequences:
        if seq in final_sequences:
            for bp in ss:
                a,b= bp
                if is_complementary(seq[a],seq[b])==False:
                    if seq in final_sequences:
                        final_sequences.remove(seq)
    return final_sequences


def sanity_check(struct_seqs,exp_seqs):
    consensus_struct = []
    consensus_exp = []
    struct_cols = zip(*struct_seqs)
    for col in struct_cols:
        consensus_struct.append(most_common(col))
    exp_cols = zip(*exp_seqs)
    for col in exp_cols:
        consensus_exp.append(most_common(col))

    if consensus_exp==consensus_struct:
        return 0
    elif consensus_exp[1:]==consensus_struct[:-1]:
        return -1
    elif consensus_exp[2:]==consensus_struct[:-2]:
        return -2
    elif consensus_exp[1:]==consensus_struct[:-1]:
        return -3

    elif consensus_exp[:-1]==consensus_struct[1:]:
        return 1
    elif consensus_exp[:-2]==consensus_struct[2:]:
        return 2
    elif consensus_exp[:-3]==consensus_struct[3:]:
        return 3

    return -1000

def get_ss_from_rfam_family(family,cols):
    ss = ""
    fname = "RF"+format(family,"05d")+".stockholm.txt"
    with open(fname,"r") as f:
        lines = f.readlines()
        full_ss=lines[-3].split()[-1]
    print("GENERATING SS")
    for i in cols:
        print(full_ss[i])
        ss = ss+full_ss[i]
    print('SS for columns : ')
    print(cols,ss)
    return ss


def alignment_column_from_seq_positions(positions_for_this_pdb, seq):
    print("GETTING ALIGNMENT COLUMNS")
    gapped_seq = seq
    s2 = str(seq)
    ungapped_seq = s2.replace("-","")

    print(gapped_seq,ungapped_seq)


    aln_pos = []

    for i in positions_for_this_pdb:

        j = gapped_column_from_ungapped_nuc_pos(gapped_seq,ungapped_seq,i)
        aln_pos.append(j)
    return(aln_pos)


def get_module_seqs_for_family(CURRENTLY_STUDIED_FAMILY,matches,module_PDB_addresses):
    newfam=True
    written = []
    seq_ranges = {}
    decales = {}
    profile = "RF"+format(CURRENTLY_STUDIED_FAMILY,"05d")+".clustal"
    bugged = ["3BBX"]
    rfam_module_seqs = []
    #print(matches)
    for family_match in matches:
        #print("WRITTEN",written)
        written = []
        if family_match=="None":
            continue
        #print(matches[i][family_match])
        if matches[family_match][0]>0 and family_match==CURRENTLY_STUDIED_FAMILY:
            for k in set(matches[family_match][1]):
                if k[0] not in bugged:
                    if newfam==True:
                        print("WRITING TO PDB")
                        match,written,seq_range,decalage = match_rfam_to_PDB(k[0]+"."+k[1],family_match, newfam=True)
                        newfam=False
                    else:
                        match,written,seq_range,decalage = match_rfam_to_PDB(k[0]+"."+k[1],family_match, written)
                    seq_ranges[".".join(k)]=seq_range
                    if decalage>0:
                        decales[".".join(k)]=decalage
                    #print(match)
    #print(positions[i])
    #print("-----------")
    written = []
    if os.stat("seq2.fasta").st_size == 0:
        print("sequence file is empty")
        return -1

    nseq = len(open("seq2.fasta").readlines(  ))
    

    if nseq>2:
        a =ClustalwCommandline("../src/clustalo", infile="seq2.fasta", outfile="seq2.clustal")
        print(a)
        a()

        print(module_PDB_addresses)
        if os.path.isfile("new_aln3.fasta"):
            os.remove("new_aln3.fasta")
        b =ClustalOmegaCommandline("../src/clustalo",profile1=profile,profile2="seq2.clustal",outfile="new_aln3.fasta",force=True)
        #./clustalw2 -PROFILE1=RF02540.clustal -SEQUENCES  -PROFILE2=seq.clustal -OUTPUT=new_alignment.clustal
        print(b)
        b()
    else:
        print(module_PDB_addresses)
        if os.path.isfile("new_aln3.fasta"):
            os.remove("new_aln3.fasta")
        b =ClustalOmegaCommandline("../src/clustalo",profile1=profile,profile2="seq2.fasta",outfile="new_aln3.fasta",force=True)
        #./clustalw2 -PROFILE1=RF02540.clustal -SEQUENCES  -PROFILE2=seq.clustal -OUTPUT=new_alignment.clustal
        print(b)
        b()        
    
    if os.path.isfile("seq2.fasta"):
        os.remove("seq2.fasta")
    if os.path.isfile("seq2.clustal"):
        os.remove("seq2.clustal") 
        
    module_columns_by_family={}
    module_columns_by_family[CURRENTLY_STUDIED_FAMILY]={}
    aln =  AlignIO.read(open("new_aln3.fasta"), "fasta")
    appearances_of_this_module_in_this_family = []
    for record in aln:
        current_PDB= record.id[3:]
        print('PDB:',current_PDB)
        if "seq" not in record.id:
            continue
        appearances_of_this_module_in_this_family.append(current_PDB)
        gapped_seq = str(record.seq)
        ungapped_seq=str(record.seq).replace("-","")
        this_example_columns = []
        if current_PDB in module_PDB_addresses:
            positions_for_this_pdb = module_PDB_addresses[current_PDB]
            print("positions",positions_for_this_pdb)
            columns_for_this_pdb = alignment_column_from_seq_positions(positions_for_this_pdb, record.seq)
            print("COLUMNS FOR THIS PDB",columns_for_this_pdb)
            module_columns_by_family[CURRENTLY_STUDIED_FAMILY][current_PDB]=columns_for_this_pdb

    #print("---------")
    print("Module columns positions by family",module_columns_by_family)
    good_cols= find_good_columns(module_columns_by_family,CURRENTLY_STUDIED_FAMILY)
    print('FOUND COLUMNS:',good_cols)
    if [-1000]==good_cols[0]:
        return []
    print("GOOD COLUMNS",good_cols)
    struct_seqs = []
    exp_seqs = []
    for PDB in module_columns_by_family[CURRENTLY_STUDIED_FAMILY]:
        columns = module_columns_by_family[CURRENTLY_STUDIED_FAMILY][PDB][0]
        sequence = module_columns_by_family[CURRENTLY_STUDIED_FAMILY][PDB][1]
        experiment_seq = []
        #print(columns)
        #print(good_cols)

        if columns in good_cols:
            #print(sequence)

            aln =  AlignIO.read(open("new_aln3.fasta"), "fasta")
            for record in aln:
                #print("PDB :", PDB)
                if PDB in record.id and -1000 not in columns:
                    #print(record.id)
                    seq = list(str(record.seq))
                    for col in columns:
                        experiment_seq.append(seq[col-0])
            #print("STRUCTURE SEQUENCE: ",sequence)
            struct_seqs.append(sequence)
            #print("RFAM SEQUENCE : ", experiment_seq)
            exp_seqs.append(experiment_seq)
    alignment_is_right = sanity_check(struct_seqs,exp_seqs)
    if alignment_is_right!=0:
        module_columns_by_family[CURRENTLY_STUDIED_FAMILY]={}
        for record in aln:
                current_PDB= record.id[3:]
                print("CURRENT PDB",current_PDB)
                if "seq" not in record.id:
                    continue
                appearances_of_this_module_in_this_family.append(current_PDB)
                gapped_seq = str(record.seq)
                ungapped_seq=str(record.seq).replace("-","")
                this_example_columns = []
                if current_PDB in module_PDB_addresses:
                    positions_for_this_pdb = module_PDB_addresses[current_PDB]
                    for pos in pos_from_fred(positions_for_this_pdb):
                        start_pos = seq_ranges[current_PDB][0]
                        this_example_columns.append(gapped_column_from_ungapped_nuc_pos(gapped_seq,ungapped_seq,pos-(start_pos)-alignment_is_right))
                    module_columns_by_family[CURRENTLY_STUDIED_FAMILY][current_PDB]=(this_example_columns,nucs_from_fred(positions_for_this_pdb))
        good_cols= find_good_columns(module_columns_by_family,CURRENTLY_STUDIED_FAMILY)
    if good_cols==[] or [-1000] in good_cols:
        return -1
    aln =  AlignIO.read(open("new_aln3.fasta"), "fasta")
    for PDB in module_columns_by_family[CURRENTLY_STUDIED_FAMILY]:
        columns = module_columns_by_family[CURRENTLY_STUDIED_FAMILY][PDB]
        sequence = module_columns_by_family[CURRENTLY_STUDIED_FAMILY][PDB][1]
        experiment_seq = []
        #print(columns)
        #print(good_cols)
        if -1000 not in columns:
            for record in aln:
                if PDB in record.id:
                    seq = list(str(record.seq))
                    for col in columns:
                        experiment_seq.append(seq[col-0])
        print("STRUCTURE SEQUENCE: ",sequence)
        print("RFAM SEQUENCE : ", experiment_seq)
    for record in aln:
        #print('FINDING PURE RFAM SEQ:',record)
        if "seq" not in record.id:
           # print("COLUMNS",good_cols)
            seq = list(str(record.seq))
            for i in range(len(good_cols[0])):
                experiment_seq = []
                columns = [col[i] for col in good_cols]
                if -1000 in columns or -1 in columns:
                    break
                for col in columns:
                    experiment_seq.append(seq[col-0])
                print("PURE RFAM SEQUENCE : ", experiment_seq)
                rfam_module_seqs.append(experiment_seq)
    return rfam_module_seqs

def sanity_check_existing_seqs():
    interesting = [2, 8, 9, 14, 20, 28, 36, 59, 113, 127, 133, 150, 162, 194, 195]
    fixed_seqs = []
    seqs = pickle.load(open("test2_rfam.cPickle",'rb'))
    graphs = pickle.load(open("test2_one_of_each_graph.cPickle",'rb'))
    for module in range(len(seqs)):
        print("MODULE:", str(interesting[module]))
        print("BEFORE REFINING:",str(len(seqs[module])))
        new_seqs = rule_out_impossible_seqs(graphs[module][0],seqs[module])
        fixed_seqs.append(new_seqs)
        print("AFTER REFINING:",str(len(fixed_seqs[module])))
        print("--------------------")
    pickle.dump(fixed_seqs,open("test3_ram.cPickle","wb"))





def PDB_to_Rfam(PDB_code, PDB_positions):
    rfam_by_PDB = {}
    number = 0
    for i in rfam:
        if len(rfam[i])>0:
            for j in rfam[i]:
                PDB_id = (j,rfam[i][j][0])
                if PDB_id not in rfam_by_PDB:
                    rfam_by_PDB[PDB_id]= [i]
                else:

                    rfam_by_PDB[PDB_id].append(i)


    rfam_extra_sequences = []
    format_aligned_modulegraphs = []
    format_PDB_positions = []
    format_graphs = []
    format_PDB_names = []


    matches={}
    positions = {}
    module_PDB_addresses = {}
    rfam_module_seqs = []
    #module_PDB = ['3D2V.B.2',"2GDI.Y.1"]
    #module_positions= {'3D2V.B' : [6,7,39,40,41,42,43,44,45,70,71,72,73],"2GDI.Y": [14,15,51,52,53,54,55,56,82,83,84,85]}
    module_PDB = PDB_code
    print(module_PDB,PDB_positions)

    module_positions = {}
    for i in range(len(module_PDB)):
        module_positions[module_PDB[i]] = [x-1 for x in PDB_positions[i]]
    #print("Rfam by PDB",rfam_by_PDB)
    matches= get_matches(rfam,rfam_by_PDB, module_PDB)
    print("matches",matches)
    #print("KEYS",list(matches.keys()))
    if os.path.isfile("seq2.fasta"):
        os.remove("seq2.fasta")
    with open("seq2.fasta", "w") as f:
        f.write("")
    time.sleep(1)
    #print("KEYS", list(matches.keys()))
    for family_match in matches.keys():
        print("CURRENTKEY",matches[family_match])
        if family_match=="None":
            continue
        if matches[family_match][0] > 0:
            # if family_match==2541:
            # print("ZIPCODE",module_addresses)
            module_seqs = get_module_seqs_for_family(family_match, matches, module_positions)
            if module_seqs == -1:
                print("FAILED TO FIND NEW SEQUENCES FOR FAMILY ", family_match)
                continue
            for seq in module_seqs:
                rfam_module_seqs.append(seq)
    print(rfam_module_seqs)
    return(rfam_module_seqs)


if __name__ == "__main__":
    PDB_pos = [6, 7, 39, 40, 41, 42, 43, 44, 45, 70, 71, 72, 73]
    a = PDB_to_Rfam(["3D2V.B"], [PDB_pos])
    print(a)