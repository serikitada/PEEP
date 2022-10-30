import argparse
import functions as fun
import pandas as pd

spacer_len = 20
context_nt = 30
gap = 6 #the gap between the positions where actual deletion take place and the ending/starting edge of forward/reverse spacer+PAM sequences
fork_value = 300 #the threshold value for the difference in the search ranges of the forward and reverse spacers, to decide whether to run flashfry once for both or twice respectively considering the computational cost.
top_value = 50 #the cutoff number to be used to narrow down the number of forward/reverse spacer candidates before making pairs, when too many candidate spacers exist 

#loading user-defined pegRNA properties --------------------------
parser = argparse.ArgumentParser(description='SoftwareName demo')
args = parser.parse_args("")

args.fasta = "example_files/ALDH2.txt"
args.ref = "example_files/ref_genome/chr22_cas9ngg_database"
args.outdir = "example_files/output_base"
args.smin = 700
args.smax = 800
args.emin = 1100
args.emax = 1300
args.lmin = 500
args.lmax = 700
args.metrics = ["DeepPE","DeepSpCas9","CFDscore", "MITscore", "mismatch_hit"]
args.filterby = [["mismatch_hit",0,1],["MITscore",50]]
args.rankby_each = ["DeepPE"]
args.rankby_pair = [["DeepPE","product"],["CFDscore","sum"]]
args.design = "TwinPE"
args.pbs_len = 11
args.twin_rtt = "ATAACTTCGTATAATGTATGCTATACGAAGTTATGGGAT"
args.scaf = "gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgc"
args.name = "ALDH2"
args.rtt_len = 30



#validating and reflecting arguments--------------------------------------------------------
outelement = fun.get_outelement(args.name, '{args.fasta}'.format(args=args))
outdir = fun.output_path(args.outdir, outelement)
tmp_outdir = fun.tmp_output_path(outdir)
logger = fun.get_log(outdir+"/"+outelement+"_base_SoftwareName.log", args) #making a log file
logger.info("Validating arguments...")
fasta = fun.read_fasta(r'{args.fasta}'.format(args=args))
seq, seql = fun.get_seq_len(fasta, spacer_len*2)
ref = args.ref
scaf = args.scaf
pbs_len = args.pbs_len
rtt_len = args.rtt_len
design = args.design
twin_rtt = args.twin_rtt
metrics_output, filterby, rankby_each, rankby_pair, metrics_all = fun.validate_metrics(args.metrics, args.filterby,args.rankby_each, args.rankby_pair, design)
#args.metrics = ["DeepPE","DeepSpCas9","CFDscore", "MITscore", "mismatch_hit"]

sminr, smin = fun.validate_smin(args.smin, seql, gap, spacer_len)
smaxr, smax = fun.validate_smax(args.smax, seql, gap, spacer_len)
eminr, emin = fun.validate_emin(args.emin, seql, gap, context_nt)
emaxr, emax = fun.validate_emax(args.emax, seql, gap,spacer_len,  context_nt)
lmin = fun.validate_lmin(args.lmin, seql)
lmax = fun.validate_lmax(args.lmax, seql)
fun.validate_start_range(sminr, smaxr, spacer_len)
fun.validate_end_range(eminr, emaxr, spacer_len)

#roughly narrow down the searching ranges for forward and reverse pegRNAs---------------
logger.info("Determining(narrowing down if applicable) the searching ranges for forward and reverse spacers...")
sminr, smaxr, eminr, emaxr, lmin,lmax = fun.narrow_search_range(seql, sminr, smaxr, eminr, emaxr, lmin, lmax, gap*2)

#abort if the search ranges are too short to possibly find spacers within------------------------
fun.validate_search_range(sminr, smaxr, spacer_len)
fun.validate_search_range(eminr, emaxr, spacer_len)

#giving context
sminr = fun.context_smin(sminr, context_nt)
smaxr = fun.context_smax(smaxr, context_nt, seql)
eminr = fun.context_emin(eminr, context_nt)
emaxr = fun.context_emax(emaxr, context_nt, seql)

#making FASTA files of the narrowed down search ranges for flashfry------------------------
flashfry_input = fun.pair_narrowed_fasta(sminr, smaxr, emaxr, eminr, seq, fork_value, tmp_outdir, outelement)

#flashfry discover module, finding potential spacers and off-target sites------------------------
logger.info("Searching for forward and reverse spacers and potential off-target sites...")
if len(flashfry_input) == 1:
    flashfry_output = [fun.flashfry_discover(flashfry_input[0], ref)]
    res_forw = pd.read_csv(flashfry_output[0], sep="\t")  
    res_rev = pd.read_csv(flashfry_output[0], sep="\t")  
else:
    flashfry_output = [fun.flashfry_discover(flashfry_input[0], ref), fun.flashfry_discover(flashfry_input[1], ref)]
    res_forw = pd.read_csv(flashfry_output[0], sep="\t")
    res_rev = pd.read_csv(flashfry_output[1], sep="\t")   

#removes forward/reverse spacers without context or not within the ranges------------------------
forw_range = fun.forw_in_range(res_forw, sminr, smin, smax, gap)
logger.info(str(forw_range.shape[0]) + " forward spacers found in range...")
forw_spacers_path_woex = tmp_outdir + "/" + outelement + "_forw_spacers"
forw_range.to_csv(forw_spacers_path_woex + ".txt", sep="\t",index=False)

rev_range = fun.rev_in_range(res_rev, eminr, emin, emax, gap)
logger.info(str(rev_range.shape[0]) + " reverse spacers found in range...")
rev_spacers_path_woex = tmp_outdir + "/" + outelement + "_rev_spacers"
rev_range.to_csv(rev_spacers_path_woex + ".txt", sep="\t",index=False)

#score forward/reverse spacers------------------------
logger.info("Scoring forward spacers...")
forw_scored = fun.score_each(seq, forw_spacers_path_woex,ref,metrics_all, scaf, pbs_len, design, gap, sminr, eminr, rtt_len,twin_rtt)
logger.info("Scoring reverse spacers...")
rev_scored = fun.score_each(seq, rev_spacers_path_woex,ref,metrics_all, scaf, pbs_len, design, gap, sminr, eminr, rtt_len,twin_rtt)

#filters bad spacers------------------------
logger.info("Filtering out forward spacers based on metrics and thresholds...")
forw_pass,forw_for_log = fun.filter(forw_scored, filterby)
for log in range(len(forw_for_log)): logger.info("---"+forw_for_log[log])
logger.info("Filtering out reverse spacers based on metrics and thresholds...")
rev_pass,rev_for_log = fun.filter(rev_scored, filterby)
for log in range(len(rev_for_log)): logger.info("---"+rev_for_log[log])

#rank spacers----------------------------------------------
logger.info("Ranking forward spacers...")
forw_ranked = fun.rank_each(forw_pass, rankby_each)
logger.info("Ranking reverse spacers...")
rev_ranked = fun.rank_each(rev_pass, rankby_each)

#make pairs and output-----------------------------------------------------------------------------
#outputs following columns: default columns, metrics-related columns, and flag-related columns (see documentation).
#default_columns: ["del_start", "del_end","del_size","FWD_spacer","RVS_spacer","FWD_PBS","FWD_PBS_length","RVS_PBS","RVS_PBS_length","FWD_RTT","FWD_RTT_length","RVS_RTT","RVS_RTT_length"]
#metrics-related columns: default or user-defined metrics for forward and reverse spacers, respectively
#flag-related columns: ["poly_T_4","SpacerGC_25_75","PBS_GC_30_60","PBS_RTT_GCrich"]
forw_top = forw_ranked[0:min(forw_ranked.shape[0],top_value)]
rev_top = rev_ranked[0:min(rev_ranked.shape[0],top_value)]


def polyT_flag(spacer_arg, pbs_arg, rtt_arg): #examine
    flag = "OK"
    extension = pbs_arg + rtt_arg
    if ("TTTT" in spacer_arg) or ("AAAA" in spacer_arg):
        flag = "FLAG"
    if ("TTTT" in extension) or ("AAAA" in extension):
        flag = "FLAG"
    return flag


def pbs_rtt_GCrich_flag(pbs_arg, rtt_arg):
    flag = "OK"
    WINDOW = 5
    extension = pbs_arg + rtt_arg
    for i in range(len(extension)-WINDOW+1):
        GC = fun.get_GCcount(extension[i:(i+WINDOW)])
        if GC >=3:
            flag = "FLAG"
    return flag


def make_pair(forw_arg, rev_arg, metrics_all_arg,design_arg,lmin_arg,lmax_arg,seq_arg,pbs_len_arg,rtt_len_arg=None,twin_rtt_arg=None):
    default_columns = ["del_start", "del_end","del_size","FWD_spacer","RVS_spacer","FWD_PBS","FWD_PBS_length","RVS_PBS","RVS_PBS_length","FWD_RTT","FWD_RTT_length","RVS_RTT","RVS_RTT_length"]
    flags_columns_raw = ["poly_T_4","SpacerGC_25_75","PBS_GC_30_60","PBS_RTT_GCrich"]
    DEF_COLUMN_NUM = 13
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
    print(column_names)
    pair = pd.DataFrame(columns = column_names)
    if design_arg == "TwinPE":
        rev_rtt = fun.reverse_complement(twin_rtt_arg)
        for x in range(forw_arg.shape[0]):
            for y in range(rev_arg.shape[0]):
                size = int(rev_arg.loc[y, 'del_end']) -  int(forw_arg.loc[x, 'del_start']) + 1
                forw_spacer = forw_arg.loc[x,'target']
                rev_spacer = rev_arg.loc[y,'target']
                forw_pbs = fun.twinpe_forw_each(seq_arg, forw_arg, x)[1]
                rev_pbs = fun.twinpe_rev_each(seq_arg, rev_arg, y)[1]
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
                    twin_rtt_arg,
                    len(twin_rtt_arg),
                    rev_rtt,
                    len(rev_rtt)]
                    metrics_part = []
                    for metric in metrics_:
                        metrics_part.append(forw_arg.loc[x,metric])
                        metrics_part.append(rev_arg.loc[y,metric])
                    flags_part = [ polyT_flag(forw_spacer, forw_pbs, twin_rtt_arg),
                    polyT_flag(rev_spacer,rev_pbs,rev_rtt),
                    fun.spacerGC_flag(forw_spacer, spacerGC_low, spacerGC_high),
                    fun.spacerGC_flag(rev_spacer, spacerGC_low, spacerGC_high),
                    fun.pbsGC_flag(forw_pbs, pbsGC_low, pbsGC_high),
                    fun.pbsGC_flag(rev_pbs, pbsGC_low, pbsGC_high),
                    pbs_rtt_GCrich_flag(forw_pbs, twin_rtt_arg),
                    pbs_rtt_GCrich_flag(rev_pbs, rev_rtt)]
                    pair.loc[pair.shape[0],:] = default_part + metrics_part + flags_part
    return pair


#logger.info("Searching for pairs with deletion length in range and assessing them...")
#args.metrics = ["DeepPE","DeepSpCas9","CFDscore", "MITscore", "mismatch_hit"]
#metrics_output, filterby, rankby_each, rankby_pair, metrics_all = fun.validate_metrics(args.metrics, args.filterby,args.rankby_each, args.rankby_pair, design)
#import functions as fun
pair = make_pair(forw_top, rev_top, metrics_all,design,lmin,lmax,seq,pbs_len,rtt_len_arg=rtt_len,twin_rtt_arg=twin_rtt)
print(pair)
logger.info(str(pair.shape[0]) + " valid pairs found...") 

logger.info("Ranking pairs...") 
pair_ranked = fun.rank_pair(pair, rankby_pair, DEF_COLUMN_NUM)

pair_ranked.to_csv(outdir + "/" + outelement + "_pairs_base_SoftwareName.csv",sep=",",index=False)
logger.info("Output files saved in ./" + outdir)
fun.remove_tmp_outdir(tmp_outdir)