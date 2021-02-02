import pickle
import os
from Bio import SeqIO
import glob

def leave_seq_out(in_file,out_file,to_remove):
    with open(out_file,"w") as outf:
        for row,record in enumerate(SeqIO.parse(in_file, "fasta")):
            id  = str(record.id)
            seq = str(record.seq)
            if row!=to_remove:
                outf.write(">"+id+"\n")
                outf.write(seq)


SEQS_TO_LEARN_2GDI = "../test/tpp_seqs_2GDIb.fasta"
SEQS_TO_LEARN_3D2V = "../test/tpp_seqs_3D2Vb.fasta"
SEQS_TO_LEARN_TLOOP = "../test/t_loop_seqs.fasta"
SEQ_TO_SEARCH = "../test/RF00059.fasta.txt"

SKIPPED_INDEX = 0
CVED_SEQS_2GDI = "../test/tpp_seqs_2GDI_loo.fasta"
CVED_SEQS_3D2V = "../test/tpp_seqs_3D2Vb_loo.fasta"
CVED_SEQS_TLOOP = "../test/t_loop_seqs_loo.fasta"

for row,record in enumerate(SeqIO.parse(SEQ_TO_SEARCH, "fasta")):
    this_seq = str(record.seq)
    leave_seq_out(SEQS_TO_LEARN_2GDI, CVED_SEQS_2GDI)
    leave_seq_out(SEQS_TO_LEARN_3D2V, CVED_SEQS_3D2V)
    leave_seq_out(SEQS_TO_LEARN_TLOOP, CVED_SEQS_TLOOP)

    files = glob.glob("../models/*tppcv*")
    if len(files)>0:
        for file in files:
            os.remove(file)

    %run ../src/module_from_desc.py -g ../test/2GDI.X.10.desc -seq ../test/tpp_seqs_2GDIb_loo.fasta -n tppcv
    %run ../src/module_from_desc.py -g ../test/3D2V.B.20.desc -seq ../test/tpp_seqs_3D2Vb_loo.fasta -n tppcv
    %run ../src/module_from_desc.py -g ../test/2GDI.T.1.desc.txt -seq ../test/t_loop_seqs_loo.fasta -n tppcv

    %run ../src/parse_sequences.py -d tppcv -seq ../test/RF00059.fasta.txt -t 14 -p 100