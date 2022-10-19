import os
import subprocess
import logging
import pandas as pd
import numpy as np
import other_soft.DeepSpCas9_modified.DeepSpCas9 as DeepSpCas9
import other_soft.DeepPE_modified.DeepPE_main as DeepPE
from Bio.Seq import Seq
#import concurrent.futures
import itertools
import multiprocessing as mp

def shell(command_arg, out=False):
    """Runs the given shell command(s). Returns the output if out=True."""
    run = subprocess.Popen(command_arg, shell=True, stdout=subprocess.PIPE)
    run.communicate()
    if out==True:
        run_out = subprocess.check_output(command_arg, shell=True)
        return run_out

#good, dont worry about it

def validate_smin(smin_arg, seql_arg, gap_arg, spacer_len_arg):
    """Validates of the -smin/--start_min argument or sets the default values when the argument is not passed"""
    sminr = 1
    smin = 1
    if smin_arg is not None: 
        if smin_arg > (seql_arg - spacer_len_arg *2):
            raise Exception("Argument -smin/--start_min too large in relation to the given sequence length, which is " + str(seql_arg) + ".")
        elif smin_arg < 1:
            raise Exception("Argument -smin/--start_min less than 1. Index starts at 1.")
        else: 
            sminr = smin_arg + gap_arg - 1 - spacer_len_arg + 1
            smin = smin_arg  
    return sminr, smin

def validate_smax(smax_arg, seql_arg, gap_arg, spacer_len_arg):
    """Validates of the argument -smax/--start_max, or sets the default values when the argument is not passed"""
    smaxr = seql_arg
    smax = seql_arg
    if smax_arg is not None:
        if smax_arg > seql_arg :
            raise Exception("Argument -smax/--start_max too large in relation to the given sequence length, which is " + str(seql_arg) + ".")
        elif smax_arg < spacer_len_arg*2 :
            raise Exception("Argument -smin/--start_min too small.")
        else: 
            smaxr = smax_arg + gap_arg - 1
            smax = smax_arg
    return smaxr, smax

def validate_emin(emin_arg, seql_arg, gap_arg, context_nt_arg):
    """Validates of the argument -emin/--end_min, or sets the default values when the argument is not passed"""
    eminr = 1
    emin = 1
    if emin_arg is not None:
        if emin_arg > (seql_arg - context_nt_arg):
            raise Exception("Argument -emin/--end_min too large in relation to the given sequence length, which is " + str(seql_arg) + ".")
        elif emin_arg < 1 : #or context
            raise Exception("Argument -emin/--end_min too small.")
        else: 
            eminr = emin_arg - gap_arg + 1
            emin = emin_arg
    return eminr, emin

def validate_emax(emax_arg, seql_arg, gap_arg, spacer_len_arg,context_nt_arg):
    """Validates the argument -emax/--end_max, or sets the default values when the argument is not passed"""
    emaxr = seql_arg
    emax = seql_arg
    if emax_arg is not None:
        if emax_arg > seql_arg:
            raise Exception("Error: Argument -emax/--end_max larger than the given sequence length, which is " + str(seql_arg) + ". Check if the indexing is right or do not specify this argument.")
        elif emax_arg < context_nt_arg*2 :
            raise Exception("Error: Argument -emax/--end_max too small.")
        else: 
            emaxr = emax_arg - gap_arg + 1 + spacer_len_arg - 1 
            emax = emax_arg
    return emaxr, emax

def validate_start_range(smin_arg, smax_arg, threshold_arg):
   if (smax_arg - smin_arg) < threshold_arg :
        raise TypeError("The distance between -smin/--start_min and -smax/--start_max is too short. If you have only specified one of them, either -smin/--start_min needs to be decresed or -smax/--start_max needs to be increased.") 

def validate_end_range(emin_arg, emax_arg, threshold_arg):
    if (emax_arg - emin_arg) < threshold_arg :
        raise TypeError("The distance between -emin/--end_min and -emax/--end_max is too short. If you have only specified one of them, either -emin/--end_min needs to be decresed or -emax/--end_max needs to be increased.")

def validate_search_range(minr_arg, maxr_arg, threshold_arg):
    if (maxr_arg - minr_arg) < threshold_arg :
        raise TypeError("Unable to search pairs fulfilling all the arguments. Check and try expanding the ranges of the deletion length, deletion start position, or deletion end position")

def validate_lmin(lmin_arg, seql_arg):
    if (lmin_arg is not None): 
        if lmin_arg > seql_arg:
            raise Exception("Error: Argument -lmin/--length_min larger than the length of given sequence")
    return lmin_arg

def validate_lmax(lmax_arg, seql_arg):
    if lmax_arg > seql_arg:
        raise Exception("Error: Argument -lmax/--length_max is larger than the length of given sequence. Check if you have put the right sequence/argument. Sequence length is " + str(seql_arg))
    return lmax_arg

def validate_lenr(lenr_arg):
    """Validates the -range/--ranges_steps"""
    lenr_list = lenr_arg
    for l in range(len(lenr_list)):
        if len(lenr_list[l])!=3:
            raise Exception("Error: the length of each list in Argument -range/--ranges_steps needs to be 3. Please refer to the manual page.")
        #if lenr_list[0] >= lenr_list[1]: #prepare for errors
    return lenr_list

# def validate_npair(npair_arg):
#     if npair_arg is not None: 
#         if npair_arg < 1:
#             raise Exception("Argument -n/--num_of_pairs less than 1.")

def narrow_search_range(seql_arg, sminr_arg, smaxr_arg, eminr_arg, emaxr_arg, lmin_arg, lmax_arg, value_arg):
    #improve
    sminr = sminr_arg
    smaxr = smaxr_arg
    eminr = eminr_arg
    emaxr = emaxr_arg
    if (lmin_arg is not None):   
        if (((eminr_arg - sminr_arg) + value_arg ) < lmin_arg): 
            eminr = sminr_arg + lmin_arg - value_arg
        if (((emaxr_arg - smaxr_arg) + value_arg*2 ) < lmin_arg): 
            smaxr = emaxr_arg - lmin_arg + value_arg
        lmin = lmin_arg
    else:
        lmin = 1
    if (lmax_arg is not None):
        if (((eminr - sminr_arg) + value_arg ) > lmax_arg): 
            sminr = eminr - lmax_arg - value_arg
        if (((emaxr - smaxr) + value_arg ) > lmax_arg): 
            emaxr = smaxr + lmax_arg + value_arg
        lmax = lmax_arg
    else:
        lmax = seql_arg
    if (sminr_arg > eminr_arg):
        eminr = sminr_arg
    if (smaxr_arg > emaxr_arg):
        smaxr = emaxr_arg

    return sminr, smaxr, eminr, emaxr, lmin, lmax

def forw_in_range(matrix_arg, sminr_arg, smin_arg,  smax_arg, gap_arg):
    forw = matrix_arg[matrix_arg['orientation'].str.contains('FWD')].copy()
    del_start = forw['stop'].astype(int) - gap_arg + 1 + sminr_arg - 1
    range = (del_start <= smax_arg) & (del_start >= smin_arg )
    forw_range = forw[range] 
    forw_range = forw_range[ forw_range['context'] != "NONE" ]
    return forw_range

def rev_in_range(matrix_arg, eminr_arg, emin_arg, emax_arg, gap_arg):
    rev = matrix_arg[matrix_arg['orientation'].str.contains('RVS')].copy()
    del_end = rev['start'].astype(int) + gap_arg - 1 + eminr_arg
    range = (del_end <= emax_arg ) &(del_end >= emin_arg )
    rev_range = rev[range]  
    rev_range = rev_range[ rev_range['context'] != "NONE" ]
    return rev_range

def read_fasta(fasta_arg):
    """Reads and returns a fasta-formatted file. Gives an error when it contains more than one sequence."""
    fasta = pd.read_csv(fasta_arg, header=None)
    if (len(np.where(fasta.iloc[:,0].str.contains(">"))[0]) > 1):
        raise Exception("fasta file should not contain more than one sequence")
    return fasta

def get_outelement(name_arg, fasta_arg, design_arg):
    """Returns the output file name element. When not specified by the user, taken from the name of the fasta file."""
    if name_arg is not None:
        outelement = name_arg
    else: 
        if "." in fasta_arg:
            outelement = '.'.join((fasta_arg.split("."))[:-1])
        else:
            outelement = fasta_arg
        if "/" in outelement:
            outelement = (outelement.split("/"))[-1] 
        outelement = outelement + "_" + design_arg
    return outelement

def get_seq_len(fasta_arg, threshold_arg):
    """Returns the sequence and the length of the given sequence."""
    namerow = np.where(fasta_arg.iloc[:,0].str.contains(">"))[0][0]
    seq = ""
    for i in range((namerow+1), fasta_arg.shape[0]):
        seq += fasta_arg.iloc[i,0]
    seql = len(seq)
    if (seql < threshold_arg):
        raise Exception("Error: sequence too short. If the sequence length is larger than "+ str(threshold_arg) +", check whether the input file is in FASTA format correctly")
    return seq, seql

def output_path(outdir_arg, outelement_arg):
    """Returns the output path and makes the directory if necessary."""
    if outdir_arg is not None: 
        outdir = outdir_arg
    else:
        outdir = '_'.join([outelement_arg,'SoftwareName','out'])
    if os.path.isdir(outdir)==False:
        os.mkdir(outdir)
    return outdir

def tmp_output_path(outdir_arg):
    """Makes a directory under the output directory for intermediate files that are to be eventually deleted."""
    if os.path.isdir(outdir_arg+"/tmp")==False:
       os.mkdir(outdir_arg+"/tmp")
       tmp_outdir = outdir_arg+"/tmp" 
    else:
        flag=0
        num=1
        while flag==0:
            if os.path.isdir(outdir_arg+"/tmp"+str(num))==True:
                num=num+1
            else:
                os.mkdir(outdir_arg+"/tmp"+str(num))
                flag=1
                tmp_outdir = outdir_arg+"/tmp"+str(num) 
    return tmp_outdir

def flashfry_discover(fastapath_arg, ref_arg):
    """Finds (possible spacers and) off-target sites using FlashFry"""
    command = "java -Xmx4g -jar other_soft/FlashFry_package/FlashFry-assembly-1.12.jar \
    discover \
    --database " + ref_arg + " \
    --fasta " + fastapath_arg + ".fasta \
    --output " + fastapath_arg + ".FlashFry.output.txt" 
    shell(command)
    path = fastapath_arg + ".FlashFry.output.txt"  
    return path

def flashfry_score(metrics_arg, spacers_path_woex_arg, ref_arg):
    """Scores spacers by defined metrics by FlashFry"""
    command = "java -Xmx4g -jar other_soft/FlashFry_package/FlashFry-assembly-1.12.jar \
    score \
    --input " + spacers_path_woex_arg + ".txt \
    --output " + spacers_path_woex_arg + ".scored.txt \
    --scoringMetrics " + metrics_arg + "\
    --database " + ref_arg
    shell(command)
    path = spacers_path_woex_arg + ".scored.txt"
    return path

def trim_fasta(seq_arg, start_arg, end_arg,tmp_outdir_arg,outelement_arg, add=""):
    header = str(">" + outelement_arg + add)
    filename = tmp_outdir_arg + "/" + outelement_arg + add + ".fasta"
    path_wo_ex =  tmp_outdir_arg + "/" + outelement_arg + add
    with open(filename, "w") as ofile:
        ofile.write(header + "\n" + seq_arg[start_arg-1:end_arg] + "\n")
    return path_wo_ex

def trim_seq(seq_arg, start_arg, end_arg):
    trimmed_seq = seq_arg[start_arg-1:end_arg]
    return trimmed_seq 

def reverse_complement(seq_arg):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    rc = "".join(complement.get(base, base) for base in reversed(seq_arg))
    return rc

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

def context_smin(sminr_arg, context_nt_arg):
    """Gives context to the minimum starting site for FlashFry. Refer to FlashFry manpage for details."""
    if (sminr_arg >= context_nt_arg):
        sminr_arg = sminr_arg - context_nt_arg
    if (sminr_arg < context_nt_arg):
        sminr_arg = 1
    return int(sminr_arg)

def context_smax(smaxr_arg, context_nt_arg, seql_arg):
    """Gives context to the maximum starting site for FlashFry. Refer to FlashFry manpage for details."""
    if (smaxr_arg + context_nt_arg <= seql_arg):
        smaxr_arg = smaxr_arg + context_nt_arg
    else:
        smaxr_arg = seql_arg
    return int(smaxr_arg)

def context_emin(eminr_arg, context_nt_arg):
    """Gives context to the minimum ending site for FlashFry. Refer to FlashFry manpage for details."""
    if (eminr_arg - context_nt_arg < 1):
        eminr_arg = 1
    else:
        eminr_arg = eminr_arg - context_nt_arg
    return int(eminr_arg)

def context_emax(emaxr_arg, context_nt_arg, seql_arg):
    """Gives context to the maximum ending site for FlashFry. Refer to FlashFry manpage for details."""
    if ((seql_arg - emaxr_arg) >= context_nt_arg):
        emaxr_arg = emaxr_arg + context_nt_arg
    if ((seql_arg - emaxr_arg) < context_nt_arg):
        emaxr_arg = seql_arg
    return int(emaxr_arg)

def remove_tmp_outdir(tmp_outdir_arg):
    for f in os.listdir(tmp_outdir_arg):
        os.remove(os.path.join(tmp_outdir_arg, f))
    os.rmdir(tmp_outdir_arg)

def get_log(logfile_path_arg, args_arg):
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",datefmt="%Y-%m-%d %H:%M:%S")
    logger = logging.getLogger(__name__)
    logfile = logging.FileHandler(filename=logfile_path_arg, mode='w')
    logger.addHandler(logfile)
    logfile.setLevel(logging.INFO)
    logger.info("=============New session=============")
    logger.info("===Timestamp in YYYY-MM-DD format===")
    fileformat = logging.Formatter("%(asctime)s %(message)s",datefmt="%Y-%m-%d %H:%M:%S")
    logfile.setFormatter(fileformat)
    logger.info("User-given arguments (default arguments applied when arguments are not provided): "+ str(vars(args_arg)))
    return logger

def pair_narrowed_fasta(sminr_arg, smaxr_arg, emaxr_arg, eminr_arg, seq_arg, fork_value_arg, tmp_outdir_arg, outelement_arg):
    if ((eminr_arg - sminr_arg) < fork_value_arg) & ((emaxr_arg - smaxr_arg) < fork_value_arg):
        both_fasta = trim_fasta(seq_arg, sminr_arg, emaxr_arg, tmp_outdir_arg, outelement_arg, add="_start_end")
        path_wo_ex = [both_fasta] 
    else:
        start_fasta = trim_fasta(seq_arg, sminr_arg, smaxr_arg, tmp_outdir_arg, outelement_arg, add="_start")
        end_fasta = trim_fasta(seq_arg, eminr_arg, emaxr_arg, tmp_outdir_arg, outelement_arg, add="_end")
        path_wo_ex = [start_fasta, end_fasta]
    return path_wo_ex

def validate_endby(endby_arg,seql_arg,fixw_arg): #improve
    if endby_arg is not None:
        endby = endby_arg
    else:
        if fixw_arg=="start":
            endby = seql_arg
        else:
            endby = 1
    return endby

def validate_metrics(metrics_args, filterby_args, rankby_each_args, rankby_pair_args,design_args):
    metrics_output = metrics_args
    filterby = filterby_args
    metrics_filterby = [metric[0] for metric in filterby]
    if design_args == "PRIME-Del":
        if "DeepPE" in metrics_filterby: raise Exception("-filter/--filter_by cannot include DeepPE when PRIME-Del is chosen. Please refer to the documentation.")
    rankby_each = rankby_each_args
    if design_args == "twinPE":
        if rankby_each is None: rankby_each = ["DeepPE","CFDscore"]
    if design_args == "PRIME-Del":
        if rankby_each is None: rankby_each = ["DeepSpCas9","CFDscore"]
        if "DeepPE" in rankby_each: raise Exception("-erankby/--rankby_each cannot be DeepPE when PRIME-Del is chosen. Please refer to the documentation.")
    rankby_pair = rankby_pair_args
    metrics_rankby_pair = [metric[0] for metric in rankby_pair]
    metrics_all_ = metrics_output + metrics_filterby + rankby_each + metrics_rankby_pair
    metrics_all = list(dict.fromkeys(metrics_all_)) #remove duplicates
    return metrics_output, filterby, rankby_each, rankby_pair, metrics_all

# def primedel(seq_arg, out_arg, pbs_arg, rtt_arg):
#     forw_pbs =  []
#     forw_rtt = []
#     rev_pbs = []
#     rev_rtt = []
#     for i in range(out_arg.shape[0]):
#         forw_pbs = np.append(forw_pbs, reverse_complement( trim_seq(seq_arg, out_arg.iloc[i,0] - pbs_arg, out_arg.iloc[i,0]-1 + 1)  ))
#         forw_rtt = np.append(forw_rtt, reverse_complement( trim_seq(seq_arg, out_arg.iloc[i,2] + 1, out_arg.iloc[i,2] + rtt_arg + 1)  )) 
#         rev_pbs = np.append(rev_pbs, reverse_complement( trim_seq(seq_arg, out_arg.iloc[i,2] + 1, out_arg.iloc[i,2] + pbs_arg + 1)  ))
#         rev_rtt = np.append(rev_rtt, reverse_complement( trim_seq(seq_arg, out_arg.iloc[i,0] - rtt_arg, out_arg.iloc[i,0] -1 + 1)  )) 
#     return forw_pbs, forw_rtt,rev_pbs,rev_rtt

# def twinpe(seq_arg, out_arg, pbs_arg):
#     forw_pbs =  []
#     rev_pbs = []
#     for i in range(out_arg.shape[0]):
#         forw_pbs = np.append(forw_pbs, reverse_complement( trim_seq(seq_arg, out_arg.iloc[i,0] - pbs_arg, out_arg.iloc[i,0]-1 + 1)  ))
#         rev_pbs = np.append(rev_pbs, reverse_complement( trim_seq(seq_arg, out_arg.iloc[i,2] + 1, out_arg.iloc[i,2] + pbs_arg + 1)  ))
#     return forw_pbs, rev_pbs 

# def primedel_forw(seq_arg, out_arg, pbs_arg, rtt_arg):
#     forw_pbs =  []
#     forw_rtt = []
#     for i in range(out_arg.shape[0]):
#         del_start = out_arg.loc[i,'del_start']
#         del_end = out_arg.loc[i,'del_end'] 
#         forw_pbs = np.append(forw_pbs, reverse_complement( trim_seq(seq_arg, del_start - pbs_arg, del_start-1 + 1)  ))
#         forw_rtt = np.append(forw_rtt, reverse_complement( trim_seq(seq_arg, del_end + 1, del_end + rtt_arg + 1)  )) 
#     return forw_pbs, forw_rtt

# def primedel_rev(seq_arg, out_arg, pbs_arg, rtt_arg):
#     rev_pbs = []
#     rev_rtt = []
#     for i in range(out_arg.shape[0]):
#         del_start = out_arg.loc[i,'del_start']
#         del_end = out_arg.loc[i,'del_end'] 
#         rev_pbs = np.append(rev_pbs, reverse_complement( trim_seq(seq_arg, del_end + 1, del_end + pbs_arg + 1)  ))
#         rev_rtt = np.append(rev_rtt, reverse_complement( trim_seq(seq_arg, del_start - rtt_arg, del_start -1 + 1)  )) 
#     return rev_pbs, rev_rtt

def primedel_both(forw_arg, rev_arg, seq_arg, pbs_len_arg,rtt_len_arg):
    del_start = forw_arg['del_start']
    del_end = rev_arg['del_end']
    forw_pbs = reverse_complement( trim_seq(seq_arg, del_start - pbs_len_arg , del_start-1 )  )
    forw_rtt = reverse_complement( trim_seq(seq_arg, del_end + 1, del_end + rtt_len_arg )  )
    rev_pbs = reverse_complement( trim_seq(seq_arg, del_end + 1, del_end + pbs_len_arg)  )
    rev_rtt = reverse_complement( trim_seq(seq_arg, del_start - rtt_len_arg, del_start -1)  )
    return forw_pbs,forw_rtt,rev_pbs,rev_rtt

def twinpe_forw(seq_arg, out_arg, pbs_arg):
    forw_pbs =  []
    for i in range(out_arg.shape[0]):
        del_start = out_arg.loc[i,'del_start'] 
        forw_pbs = np.append(forw_pbs, reverse_complement( trim_seq(seq_arg, del_start - pbs_arg, del_start-1)  ))
    return forw_pbs

def twinpe_forw_each(seq_arg, out_arg,i_arg):
    del_start = out_arg.loc[i_arg,'del_start'] 
    pbs_len = out_arg.loc[i_arg,'PBS_len'] 
    forw_pbs = reverse_complement( trim_seq(seq_arg, del_start - pbs_len, del_start-1)  )
    return pbs_len,forw_pbs

def twinpe_rev(seq_arg, out_arg, pbs_arg): 
    rev_pbs = []
    for i in range(out_arg.shape[0]):
        del_end = out_arg.loc[i,'del_end'] 
        rev_pbs = np.append(rev_pbs, reverse_complement( trim_seq(seq_arg, del_end + 1, del_end + pbs_arg)  ))
    return rev_pbs

def twinpe_rev_each(seq_arg, out_arg,i_arg):
    del_end = out_arg.loc[i_arg,'del_end'] 
    pbs_len = out_arg.loc[i_arg,'PBS_len'] 
    forw_pbs = reverse_complement( trim_seq(seq_arg, del_end + 1, del_end + pbs_len )  )
    return pbs_len,forw_pbs

def get_DeepSpCas9(seq_arg, spacers_matrix_arg):
    NEIGHBOR_LEFT = 4
    NEIGHBOR_RIGHT = 25
    score = []
    for spacer in range(spacers_matrix_arg.shape[0]):
        fwd_or_rvs = spacers_matrix_arg.loc[spacer, 'orientation']
        if fwd_or_rvs == 'FWD':
            start = spacers_matrix_arg.loc[spacer, 'start']
            seq30 = seq_arg[start-NEIGHBOR_LEFT-1: start+NEIGHBOR_RIGHT] 
        else:
            start = spacers_matrix_arg.loc[spacer, 'stop']
            seq30 = reverse_complement(seq_arg[start-NEIGHBOR_RIGHT-1 : start+NEIGHBOR_LEFT])
        score_each = DeepSpCas9.main(seq30)
        score.append(score_each)
    return score

def get_DeepPE_twin(seq_arg, spacers_matrix_arg, scaf_arg, pbs_len_arg, twin_rtt_arg):
    NEIGHBOR_LEFT = 4
    NEIGHBOR_RIGHT = 42
    PBS_SPACE = 21
    RTT_SPACE = 26
    score = []
    rtt = twin_rtt_arg
    for spacer in range(spacers_matrix_arg.shape[0]):
        fwd_or_rvs = spacers_matrix_arg.loc[spacer, 'orientation']
        if fwd_or_rvs == 'FWD':
            start = spacers_matrix_arg.loc[spacer,'start']
            seq47 = seq_arg[start-NEIGHBOR_LEFT-1 : start+NEIGHBOR_RIGHT] 
            pbs = twinpe_forw(seq_arg, spacers_matrix_arg, pbs_len_arg)[spacer]
        else:
            start = spacers_matrix_arg.loc[spacer,'stop']
            seq47 = reverse_complement(seq_arg[start-NEIGHBOR_RIGHT-1 : start+NEIGHBOR_LEFT])
            pbs = twinpe_rev(seq_arg, spacers_matrix_arg, pbs_len_arg)[spacer]
            rtt = reverse_complement(rtt)
        if len(pbs) > PBS_SPACE: pbs = pbs[0:PBS_SPACE]
        if len(rtt) > RTT_SPACE: rtt = rtt[0:RTT_SPACE]
        score_each = DeepPE.main(seq47_arg=seq47, pbs_arg=pbs, rtt_arg=rtt, scaf_arg=scaf_arg)
        score.append(score_each)
    return score

def get_DeepPE_primedel_each(seq_arg, spacer_vector_arg,pbs_arg, rtt_arg, scaf_arg):
    NEIGHBOR_LEFT = 4
    NEIGHBOR_RIGHT = 42
    PBS_SPACE = 21
    RTT_SPACE = 26
    pbs = pbs_arg
    rtt = rtt_arg
    fwd_or_rvs =  spacer_vector_arg['orientation']
    if fwd_or_rvs == 'FWD':
        start = spacer_vector_arg['start']
        seq47 = seq_arg[start-NEIGHBOR_LEFT-1 : start+NEIGHBOR_RIGHT] 
    else:
        start = spacer_vector_arg['stop']
        seq47 = reverse_complement(seq_arg[start-NEIGHBOR_RIGHT-1 : start+NEIGHBOR_LEFT])
    if len(pbs) > PBS_SPACE: pbs = pbs[0:PBS_SPACE]
    if len(rtt) > RTT_SPACE: rtt = rtt[0:RTT_SPACE]
    score = DeepPE.main(seq47_arg=seq47, pbs_arg=pbs, rtt_arg=rtt, scaf_arg=scaf_arg)
    return score

def score_each(seq_arg, spacers_path_woex_arg,ref_arg,metrics_arg, scaf_arg, pbs_len_arg, design_arg, gap_arg, wminr_arg, rtt_len_arg=None,twin_rtt_arg = None):
    flashfry_metrics = "minot"
    if ("CRISPRscan" in metrics_arg): flashfry_metrics += ",moreno2015"
    if ("RuleSet1" in metrics_arg): flashfry_metrics += ",doench2014ontarget"
    if ("CFDscore" in metrics_arg): flashfry_metrics += ",doench2016cfd"
    if ("MITscore" in metrics_arg): flashfry_metrics += ",hsu2013"
    scored_path = flashfry_score(flashfry_metrics, spacers_path_woex_arg, ref_arg)
    scored = pd.read_csv(scored_path, sep="\t")
    print("until setting paths and metrics")
    if scored.shape[0] != 0 : #== 0 happens when FlashFry "overflow"s - https://github.com/mckennalab/FlashFry/wiki/Discovery-options
        fwd_or_rvs = scored.loc[0, 'orientation']
        if fwd_or_rvs == 'FWD':
            scored['del_start'] = scored['stop'].astype(int) - gap_arg + 1 + wminr_arg - 1
        else:
            scored['del_end'] = scored['start'].astype(int) + gap_arg - 1 + wminr_arg 
        scored['PBS_len'] = pbs_len_arg
        if rtt_len_arg is not None: scored['RTT_len'] = rtt_len_arg
        if twin_rtt_arg is not None: scored['RTT_len'] = len(twin_rtt_arg)
        if ("DeepSpCas9" in metrics_arg):
            scored['DeepSpCas9'] = get_DeepSpCas9(seq_arg, scored)
        if ("DeepPE" in metrics_arg):
            if design_arg == 'twinPE':
                scored['DeepPE'] = get_DeepPE_twin(seq_arg, scored, scaf_arg, pbs_len_arg, twin_rtt_arg)
    return scored

def filter(scored_spacers_matrix_arg, filterby_arg):
    correspond = {"DeepPE":"DeepPE","DeepSpCas9":"DeepSpCas9","CRISPRscan":"Moreno2015","RuleSet1":"Doench2014OnTarget","CFDscore":"DoenchCFD_specificityscore","MITscore":"Hsu2013"}
    passed = scored_spacers_matrix_arg

    n_filter = len(filterby_arg)
    for_log = []
    for n in range(n_filter):
        before = passed.shape[0]
        filter = filterby_arg[n][0]
        threshold = float(filterby_arg[n][1])
        if filter=="mismatch_hit":
            filter_ = "0-1-2-3-4_mismatch"
            mismatch_threshold = int(threshold)
            hit_threshold = int(filterby_arg[n][2])
            hit = []
            for i in range(passed.shape[0]):
                hit.append(int(passed.loc[:,filter_].values[i].split(",")[int(mismatch_threshold)]))
            passed = passed[np.array(hit) < hit_threshold]
            after = passed.shape[0]
            for_log.append(filter+" (mismatch"+str(mismatch_threshold)+" < "+str(hit_threshold)+"): "+str(before-after)+" spacers filtered out for having "+" >= "+str(hit_threshold)+" off-target sites with "+str(mismatch_threshold)+" mismatch(es)")
        else:
            #print(filter)
            filter_ = correspond.get(filter)
            #print('filter:',filter_)
            passed = passed[passed.loc[:,filter_].values >= threshold]
            after = passed.shape[0]
            for_log.append(filter + "( >= " + str(threshold) + "): " + str(before-after) + " spacers filtered out for its score being " + " < "+str(threshold))
        passed = passed.reset_index(drop=True)
    return passed,for_log

def rank_each(passed_spacers_matrix_arg, rankby_each_arg):
    correspond = {"DeepPE":"DeepPE","DeepSpCas9":"DeepSpCas9","CRISPRscan":"Moreno2015","RuleSet1":"Doench2014OnTarget","CFDscore":"DoenchCFD_specificityscore","MITscore":"Hsu2013"}
    passed = passed_spacers_matrix_arg
    rankby_each = rankby_each_arg
    rankby_each_ = []
    for metric in rankby_each:
        metric_ = correspond.get(metric)
        rankby_each_.append(metric_)
    ascending = [False]*len(rankby_each_) #all scores higher the better
    passed_ranked = passed.sort_values(rankby_each_,ascending=ascending) 
    passed_ranked = passed_ranked.reset_index(drop=True)
    return passed_ranked

def polyT_flag(spacer_arg, pbs_arg, rtt_arg): #examine
    flag = "OK"
    extension = pbs_arg + rtt_arg
    if ("TTTT" in spacer_arg) or ("AAAA" in spacer_arg):
        flag = "FLAG"
    if ("TTTT" in extension) | ("AAAA" in extension):
        flag = "FLAG"
    return flag

def spacerGC_flag(seq_arg, spacerGC_low_arg, spacerGC_high_arg):
    flag = "OK"
    GC = get_GCcontent(seq_arg)
    if (GC <= spacerGC_low_arg) or (GC >= spacerGC_high_arg):
        flag = "FLAG"
    return flag

def pbsGC_flag(seq_arg, pbsGC_low_arg, pbsGC_high_arg):
    flag = "OK"
    GC = get_GCcontent(seq_arg)
    if (GC <= pbsGC_low_arg) or (GC >= pbsGC_high_arg):
        flag = "FLAG"
    return flag

def pbs_rtt_GCrich_flag(pbs_arg, rtt_arg):
    flag = "OK"
    WINDOW = 5
    extension = pbs_arg + rtt_arg
    for i in range(len(extension)-WINDOW+1):
        GC = get_GCcount(extension[i:(i+WINDOW)])
        if GC >=3:
            flag = "FLAG"
    return flag


##################OPTIMIZATION##################

#fasta = read_fasta("example_files/sequence/MAPK1.fna")
#seq, seql = get_seq_len(fasta, 20*2)
#design = "PRIME-Del"
#metrics_output, filterby, rankby_each, rankby_pair, metrics_all = validate_metrics("DeepPE DeepSpCas9 CFDscore MITscore mismatch_hit",[["mismatch_hit",0,2],["MITscore",50],["CFDscore",0.8]],"DeepSpCas9 CFDscore", [["DeepPE","product"],["CFDscore","sum"]])
#forw_arg = pd.read_csv("~/Downloads/Sanger/01.Project/gitrepos/my_repositories/FlashFry/example_files/output_base/opt_forw_MAPK1_PRIME-Del_example_single_SoftwareName.csv", sep=",")  
#rev_arg = pd.read_csv("~/Downloads/Sanger/01.Project/gitrepos/my_repositories/FlashFry/example_files/output_base/opt_rev_MAPK1_PRIME-Del_example_single_SoftwareName.csv", sep=",") 

def make_pair(forw_arg, rev_arg, metrics_all_arg,design_arg,lmin_arg,lmax_arg,seq_arg,pbs_len_arg,scaf_arg,rtt_len_arg=None,twin_rtt_arg=None):
    default_columns = ["del_start", "del_end","del_length","FWD_spacer","RVS_spacer","FWD_PBS","FWD_PBS_length","RVS_PBS","RVS_PBS_length","FWD_RTT","FWD_RTT_length","RVS_RTT","RVS_RTT_length"]
    flags_columns_raw = ["poly_T_4","SpacerGC_25_75","PBS_GC_30_60","PBS_RTT_GCrich"]
    spacerGC_low = 25
    spacerGC_high = 75
    pbsGC_low = 30
    pbsGC_high = 60
    flags_columns = []
    for flag in flags_columns_raw: 
        flags_columns.append(str("FWD_"+flag))
        flags_columns.append(str("RVS_"+flag))
    correspond = {"DeepPE":"DeepPE","DeepSpCas9":"DeepSpCas9","CRISPRscan":"Moreno2015","RuleSet1":"Doench2014OnTarget","CFDscore":"DoenchCFD_specificityscore","MITscore":"Hsu2013","mismatch_hit":"0-1-2-3-4_mismatch"}

    metrics_columns_raw = metrics_all_arg
    metrics_columns = []
    metrics_ = []
    for metric in metrics_columns_raw: 
        metrics_columns.append(str("FWD_"+metric))
        metrics_columns.append(str("RVS_"+ metric))
        metrics_.append(correspond.get(metric))

    column_names = default_columns + metrics_columns + flags_columns 
    pair = pd.DataFrame(columns = column_names)
    if design_arg == "PRIME-Del":
        if "DeepPE" in metrics_all_arg:
            forw_arg["DeepPE"] = None
            rev_arg["DeepPE"] = None
    if design_arg == "twinPE":
        forw_rtt = twin_rtt_arg
        rev_rtt = reverse_complement(twin_rtt_arg)

    pool = mp.Pool(mp.cpu_count())

    
    for x, y in itertools.product(range(forw_arg.shape[0]), range(rev_arg.shape[0])):
        size = int(rev_arg.loc[y, 'del_end']) -  int(forw_arg.loc[x, 'del_start']) + 1
        forw_spacer = forw_arg.loc[x,'target']
        rev_spacer = rev_arg.loc[y,'target']
        
        if design_arg == "PRIME-Del":
            forw_pbs,forw_rtt,rev_pbs,rev_rtt = primedel_both(forw_arg.loc[x,:], rev_arg.loc[y,:], seq_arg, pbs_len_arg,rtt_len_arg)
            if "DeepPE" in metrics_all_arg: 
                forw_arg.loc[x,"DeepPE"] = get_DeepPE_primedel_each(seq_arg, forw_arg.loc[x,:],forw_pbs, forw_rtt, scaf_arg)
                rev_arg.loc[y,"DeepPE"] = get_DeepPE_primedel_each(seq_arg, rev_arg.loc[y,:],rev_pbs, rev_rtt, scaf_arg)
        
        if design_arg == "twinPE":
            forw_pbs = twinpe_forw_each(seq_arg, forw_arg, x)[1]
            rev_pbs = twinpe_rev_each(seq_arg, rev_arg, y)[1]
        
        if (size <= lmax_arg ) & (size >= lmin_arg ):
            default_part = [ forw_arg.loc[x,'del_start'], 
            rev_arg.loc[y,'del_end'],
            size,
            forw_spacer, 
            rev_spacer,
            forw_pbs,
            len(forw_pbs),
            rev_pbs,
            len(rev_pbs),
            forw_rtt,
            len(forw_rtt),
            rev_rtt,
            len(rev_rtt)]
            metrics_part = []
            for metric in metrics_:
                metrics_part.append(forw_arg.loc[x,metric])
                metrics_part.append(rev_arg.loc[y,metric])
            print("until metrics metrics")
            flags_part = [ polyT_flag(forw_spacer, forw_pbs, forw_rtt),
            polyT_flag(rev_spacer,rev_pbs,rev_rtt),
            spacerGC_flag(forw_spacer, spacerGC_low, spacerGC_high),
            spacerGC_flag(rev_spacer, spacerGC_low, spacerGC_high),
            pbsGC_flag(forw_pbs, pbsGC_low, pbsGC_high),
            pbsGC_flag(rev_pbs, pbsGC_low, pbsGC_high),
            pbs_rtt_GCrich_flag(forw_pbs, forw_rtt),
            pbs_rtt_GCrich_flag(rev_pbs, rev_rtt)]

            pair.loc[pair.shape[0],:] = default_part + metrics_part + flags_part
        pool.close() 
    return pair


#make_pair(forw_arg, rev_arg, ["DeepPE","DeepSpCas9","CFDscore","MITscore","mismatch_hit"],design,23500,26500,seq,11,"gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgc",rtt_len_arg=11,twin_rtt_arg=30)



def rank_pair(passed_pairs_matrix_arg, rankby_pair_arg):
    DEF_COLUMN_NUM = 13 # the number of default columns in make_pair() function above
    passed = passed_pairs_matrix_arg
    rankby_pair = []
    passed_added = passed
    for metric in range(len(rankby_pair_arg)):
        m1 = rankby_pair_arg[metric][0]
        m2 = rankby_pair_arg[metric][1]
        col_name = "PAIR_"+m2+"_"+m1
        fwd = passed.loc[:,"FWD_"+m1]
        rvs = passed.loc[:,"RVS_"+m1]
        if m2=="sum": pmetric = pd.DataFrame({col_name: list(fwd+rvs)})
        if m2=="product":pmetric = pd.DataFrame({col_name: list(fwd*rvs)})
        if m2=="min":pmetric = pd.DataFrame({col_name: list(np.minimm(fwd,rvs))})
        passed_added = pd.concat([passed_added.iloc[:, 0:DEF_COLUMN_NUM+metric],pmetric,passed_added.iloc[:, DEF_COLUMN_NUM+metric:passed_added.shape[1]+1]],axis=1)
        rankby_pair.append(col_name)
    ascending = [False]*len(rankby_pair) #all scores higher the better
    passed_ranked = passed_added.sort_values(rankby_pair,ascending=ascending) 
    return passed_ranked

def read_extend_input(input_path_arg):
    list = pd.read_csv(input_path_arg, header=0)
    query_names = list.iloc[:,0].values.tolist()
    index = 1
    for i in np.where(pd.isnull(query_names).tolist())[0].tolist():
        query_names[i] = str("NoQueryName"+str(index))
        index = index + 1
    len_list = list.iloc[:,1].values.tolist()
    error_list = list.iloc[:,2].values.tolist()
    for i in np.where(np.isnan(error_list))[0].tolist():
        error_list[i] = round(len_list[i]*0.05)
    return query_names, len_list, error_list

#def compare_extend_sets(set_names_arg, outelement_arg, pair_rankby_arg):