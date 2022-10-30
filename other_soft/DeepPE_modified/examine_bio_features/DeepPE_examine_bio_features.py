from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import RNA
import numpy as np
import pandas as pd
import re

def reverse_complement(seq_arg):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    rc = "".join(complement.get(base, base) for base in reversed(seq_arg))
    return rc

def complement(seq_arg):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    c = "".join(complement.get(base, base) for base in seq_arg)
    return c

def get_mfe(x):
    """ Calculates secondary structure approximated by Vienna Fold Energy. """
    vf = RNA.fold(x)
    if vf[1] == np.nan:
        return 0
    else:
        return vf[1]

def get_tm(string_arg):
    """ Calculates the melting temperature of x based on Biopython's TM_NN."""
    myseq = Seq(string_arg)
    return mt.Tm_NN(myseq, nn_table=mt.R_DNA_NN1)

def get_tm2(string_arg):
    """ Calculates the melting temperature of x based on Biopython's TM_NN."""
    myseq = Seq(string_arg)
    return mt.Tm_NN(myseq)

def get_tm3(string_arg, c_seq_arg):
    """ Calculates the melting temperature of x based on Biopython's TM_NN."""
    myseq = Seq(string_arg)
    c_seq = Seq(c_seq_arg)
    return mt.Tm_NN(myseq, c_seq = c_seq)

def get_GCcount(string_arg):
    """ Calculates GC count."""
    myseq = Seq(string_arg)
    GCcount = myseq.count("C") + myseq.count("G") 
    return GCcount

def get_GCcontent(string_arg):
    """ Calculates GC content (percentage)."""
    myseq = Seq(string_arg)
    GCcount = myseq.count("C") + myseq.count("G")
    seq_len = len(myseq)
    GCcontent = GCcount/seq_len*100
    return GCcontent
    
#def get_deep_sp_cas9(seq_30nt_arg):
#    """Gets DeepSpCas9 score"""


#arguments-----------------------------------------------------------------------------------------------------------------------------------
left_neighbor = 4
spacer_len = 20
1
GtoC_pos = 5
scaffold = "gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgc"
filepath = "./DeepPE_test_input.txt"
test_input = pd.read_csv(filepath, sep="\t")

#Get sequences correspoonding to "wide target sequence" and "prime edited sequence" inputs in the DeepPE model-------------------------------
for row in range(test_input.shape[0]):
    line = test_input.iloc[row]
    wide,edit,pbs_len,rtt_len, pbs_rtt_len = line[0],line[1],int(line[2]),int(line[3]),int(line[4])
    spacer = wide[left_neighbor:(left_neighbor+spacer_len)]
    pbs_rtt_cor = re.sub("X","",edit)
    pbs_cor = pbs_rtt_cor[0:pbs_len]
    pbs = reverse_complement(pbs_cor)
    rtt_cor = pbs_rtt_cor[pbs_len:pbs_rtt_len]
    rtt = reverse_complement(rtt_cor)
    rtt_cor_woGtoC = rtt_cor[0:GtoC_pos-1] + "G" + rtt_cor[GtoC_pos:pbs_rtt_len]
    PAMop_woGtoC = complement(rtt_cor_woGtoC)
    PAMop_wGtoC = complement(rtt_cor)

    #Calculate 20 biological features------------------------------------------------------------------------------------------------------------
    #1-3 (lengths)
    # pbs_len = len(pbs_cor)
    # RTT_len = len(RTT_cor)
    # PBS_RTT_len = len(PBS_RTT_cor)
    # print("PBS_len:"+str(PBS_len)+", RTT_len:"+str(RTT_len)+", PBS_RTT_len :"+str(PBS_RTT_len ))

    #4-8 (melting temperatres)
    Tm1 = round(get_tm(pbs),4)
    Tm2 = get_tm2(rtt_cor_woGtoC)
    Tm3 = get_tm3(rtt_cor_woGtoC, complement(rtt_cor))
    Tm4 = get_tm(rtt)
    deltaTm = Tm3 - Tm2
    #print("Tm1:"+str(Tm1)+", Tm2:"+str(Tm2)+", Tm3:"+str(Tm3)+", Tm4:"+str(Tm4)+", deltaTm:"+str(deltaTm))
    test_input.iloc[row,5] = round(Tm1,4)
    test_input.iloc[row,6] = round(Tm2,4)
    test_input.iloc[row,7] = round(Tm3,4)
    test_input.iloc[row,8] = round(Tm4,4)
    test_input.iloc[row,9] = round(deltaTm,4)

    #9-14 (GCs)
    GCcount1 = get_GCcount(pbs_cor)
    GCcount2 = get_GCcount(rtt_cor)
    GCcount3 = get_GCcount(pbs_rtt_cor)
    GCcontent1 = get_GCcontent(pbs_cor)
    GCcontent2 = get_GCcontent(rtt_cor)
    GCcontent3 = get_GCcontent(pbs_rtt_cor)
    #print("GCcount1:"+str(GCcount1)+", GCcount2:"+str(GCcount2)+", GCcount3:"+str(GCcount3)+", GCcontent1:"+str(GCcontent1)+", GCcontent2:"+str(GCcontent2)+", GCcontent3:"+str(GCcontent3))
    test_input.iloc[row,10] = GCcount1
    test_input.iloc[row,11] = GCcount2
    test_input.iloc[row,12] = GCcount3
    test_input.iloc[row,13] = round(GCcontent1,8)
    test_input.iloc[row,14] = round(GCcontent2,8)
    test_input.iloc[row,15] = round(GCcontent3,8)

    #15-19 (MFEs)
    mfe_1 = get_mfe("G" + spacer[1:20] + scaffold + rtt + pbs + "TTTTTTT")
    mfe_2 = get_mfe(scaffold + rtt + pbs + "TTTTTTT")
    mfe_3 = get_mfe(rtt + pbs + "TTTTTTT")
    mfe_4 = get_mfe("G" + spacer[1:20])
    mfe_5 = get_mfe("G" + spacer[1:20] + scaffold)
    #print("MFE_1:"+str(MFE_1)+", MFE_2:"+str(MFE_2)+", MFE_3:"+str(MFE_3)+", MFE_4:"+str(MFE_4)+", MFE_5:"+str(MFE_5))
    test_input.iloc[row,16] = round(mfe_1,1)
    test_input.iloc[row,17] = round(mfe_2,1)
    test_input.iloc[row,18] = round(mfe_3,1)
    test_input.iloc[row,19] = round(mfe_4,1)
    test_input.iloc[row,20] = round(mfe_5,1)

    #20 (DeepSpCas9)
    #20 (DeepSpCas9)
    #deep_sp_cas9 = DeepSpCas9.main(wide[0:30])

test_input.to_csv(filepath, sep="\t",index=False)

