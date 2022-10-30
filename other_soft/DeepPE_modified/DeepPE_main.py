import warnings
warnings.filterwarnings('ignore')
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import RNA
import numpy as np
import pandas as pd
import os
import sys
PACKAGE_DIRECTORY = os.path.dirname(os.path.realpath(__file__))
sys.path.append(PACKAGE_DIRECTORY)
import DeepPE_model
parent = os.path.dirname(PACKAGE_DIRECTORY)
sys.path.append(parent+"/DeepSpCas9_modified/")
import DeepSpCas9

## Paths ##
#PACKAGE_DIRECTORY = os.path.dirname(os.path.realpath(__file__))
#MODEL_DIR = PACKAGE_DIRECTORY + '/DeepPE_weights/'
#best_model_path_list = [MODEL_DIR]
TEST_NUM_SET         = [0] # List can be expanded in case of multiple test parameters

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


def main(seq47_arg=None, pbs_arg=None, rtt_arg=None, scaf_arg="gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgc"):
    #arguments-----------------------------------------------------------------------------------------------------------------------------------
    #parser = argparse.ArgumentParser(description='DeepPE')
    #optional = parser.add_argument_group('optional arguments')
    #optional.add_argument('-seq47', '--seq_47nt', dest = 'seq47', type=str, help ='4-bp left neighbor, 20-bp protospacer, 3-bpPAM,and 20-bp right neighbor sequences', required=False)
    #optional.add_argument('-pbs', '--pbs', dest = 'pbs', type=str, help ='PBS sequence', required=False)
    #optional.add_argument('-rtt', '--rtt', dest = 'rtt', type=str, help ='RTT sequence', required=False)
    #optional.add_argument('-scaf', '--scaffold', dest = 'scaf', type = str, help ='Scaffold for PE2', default="gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgc",required=False)
    #args = parser.parse_args()
    #if args.seq47 is not None:
    #    seq47_arg = args.seq47
    #if args.pbs is not None:
    #    pbs_arg = args.pbs
    #if args.rtt is not None:
    #    rtt_arg = args.rtt
    #if args.scaf is not None:
    #    scaf_arg = args.scaf
    
    scaffold = scaf_arg
    wide = seq47_arg
    pbs = pbs_arg
    rtt = rtt_arg
    pbs_len = len(pbs)
    rtt_len = len(rtt)
    pbs_rtt_len = pbs_len + rtt_len
    pbs_cor = reverse_complement(pbs)
    rtt_cor = reverse_complement(rtt)
    pbs_rtt_cor = pbs_cor + rtt_cor
    PBS_SPACE = 21 #refer to DeepPE second input ("edit" sequence, 47nt, filled with "X")
    RTT_SPACE = 26 #refer to DeepPE second input ("edit" sequence, 47nt, filled with "X")
    edit = "".join("X" for i in range(PBS_SPACE - pbs_len)) + pbs_rtt_cor + "".join("X" for i in range(RTT_SPACE - rtt_len)) #refer to DeepPE second input ("edit" sequence, 47nt, filled with "X")
    
    left_neighbor = 4
    spacer_len = 20
    spacer = wide[left_neighbor:(left_neighbor+spacer_len)]
    
    #GtoC_pos = 5 
    #rtt_cor_woGtoC = rtt_cor[0:GtoC_pos-1] + "G" + rtt_cor[GtoC_pos:pbs_rtt_len]

    #Calculate 20 biological features------------------------------------------------------------------------------------------------------------
    #1-3 (lengths)
    # pbs_len = len(pbs_cor)
    # RTT_len = len(RTT_cor)
    # PBS_RTT_len = len(PBS_RTT_cor)
    # print("PBS_len:"+str(PBS_len)+", RTT_len:"+str(RTT_len)+", PBS_RTT_len :"+str(PBS_RTT_len ))

    #4-8 (melting temperatres)
    Tm1 = get_tm(pbs)
    #Tm2 = get_tm2(rtt_cor_woGtoC)
    #Tm3 = get_tm3(rtt_cor_woGtoC, complement(rtt_cor))
    Tm2 = get_tm2(rtt_cor)
    Tm3 = get_tm3(rtt_cor, complement(rtt_cor))
    Tm4 = get_tm(rtt)
    deltaTm = Tm3 - Tm2
    

    #9-14 (GCs)
    GCcount1 = get_GCcount(pbs_cor)
    GCcount2 = get_GCcount(rtt_cor)
    GCcount3 = get_GCcount(pbs_rtt_cor)
    GCcontent1 = get_GCcontent(pbs_cor)
    GCcontent2 = get_GCcontent(rtt_cor)
    GCcontent3 = get_GCcontent(pbs_rtt_cor)
    #print("GCcount1:"+str(GCcount1)+", GCcount2:"+str(GCcount2)+", GCcount3:"+str(GCcount3)+", GCcontent1:"+str(GCcontent1)+", GCcontent2:"+str(GCcontent2)+", GCcontent3:"+str(GCcontent3))

    #15-19 (MFEs)
    mfe_1 = get_mfe("G" + spacer[1:20] + scaffold + rtt + pbs + "TTTTTTT")
    mfe_2 = get_mfe(scaffold + rtt + pbs + "TTTTTTT")
    mfe_3 = get_mfe(rtt + pbs + "TTTTTTT")
    mfe_4 = get_mfe("G" + spacer[1:20])
    mfe_5 = get_mfe("G" + spacer[1:20] + scaffold)
    #print("MFE_1:"+str(MFE_1)+", MFE_2:"+str(MFE_2)+", MFE_3:"+str(MFE_3)+", MFE_4:"+str(MFE_4)+", MFE_5:"+str(MFE_5))

    #20 (DeepSpCas9)
    deep_sp_cas9 = DeepSpCas9.main(wide[0:30])
        
    #DeepPE
    score = DeepPE_model.main(wide, edit, pbs_len, rtt_len, pbs_rtt_len, round(Tm1,4),round(Tm2,4),round(Tm3,4),round(Tm4,4),round(deltaTm,4),GCcount1,GCcount2,GCcount3,round(GCcontent1,8),round(GCcontent2,8),round(GCcontent3,8),round(mfe_1,1), round(mfe_2,1),round(mfe_3,1),round(mfe_4,1),round(mfe_5,1),round(deep_sp_cas9,3))
    #print(score)
    #print(wide, edit, pbs_len, rtt_len, pbs_rtt_len, round(Tm1,4),round(Tm2,4),round(Tm3,4),round(Tm4,4),round(deltaTm,4),GCcount1,GCcount2,GCcount3,round(GCcontent1,8),round(GCcontent2,8),round(GCcontent3,8),round(mfe_1,1), round(mfe_2,1),round(mfe_3,1),round(mfe_4,1),round(mfe_5,1),round(deep_sp_cas9,3))
    print(score)

    

if __name__ == '__main__':
    main()