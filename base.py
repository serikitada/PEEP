import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import functions as fun
import argparse
import pandas as pd
import json


def main():
    spacer_len = 20
    context_nt = 30
    gap = 6 #the gap between the positions where actual deletion take place and the ending/starting edge of forward/reverse spacer+PAM sequences
    fork_value = 300 #the threshold value for the difference in the search ranges of the forward and reverse spacers, to decide whether to run flashfry once for both or twice respectively considering the computational cost.
    top_value = 50 #the cutoff number to be used to narrow down the number of forward/reverse spacer candidates before making pairs, when too many candidate spacers exist 

    defaults = {
        'metrics': ["DeepPE","DeepSpCas9","MITscore","CFDscore","mismatch_hit"], 
        'filter_by': [["mismatch_hit",0,2],["MITscore",50]],
        'erankby': ['DeepSpCas9', 'CFDscore'], 
        'prankby': [["DeepPE", "product"],["CFDscore","sum"]],
        'pbsl': 11, 
        'rttl': 30,
        'scaf':"gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgc",
        'twin_rtt': None
        }
    
    #loading user-defined pegRNA properties --------------------------
    parser = argparse.ArgumentParser(description='PEEP demo')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')  

    #module specific arguments
    optional.add_argument('-smin', '--start_min', dest = 'smin', type = int, help ='Minimum starting position of the deletion, index starting from 1', required=False)
    optional.add_argument('-smax', '--start_max', dest = 'smax', type = int, help = 'Maximum starting position of the deletion', required=False)
    optional.add_argument('-emin', '--end_min', dest = 'emin', type = int, help ='Minimum ending position of the deletion', required=False)
    optional.add_argument('-emax', '--end_max', dest = 'emax', type = int, help = 'Maximum ending position of the deletion', required=False)
    optional.add_argument('-lmin', '--length_min', dest = 'lmin', type = int, help ='Minimum length of the deletion', required=False)
    optional.add_argument('-lmax', '--length_max', dest = 'lmax', type = int, help = 'Maximum length of the deletion', required=False)
    #optional.add_argument('-n', '--numbers_of_pairs', dest = 'npair', type = int, help ='Number of pairs to return', required=False)

    #common arguments 
    required.add_argument('-f', '--fasta', dest = 'fasta', type=str, help ='(The path to)the FASTA file', required=True)
    required.add_argument('-ref', '--ref_genome', dest = 'ref', type = str, help ='(The path to) the reference genome already indexed by "index.py"', required=True)
    required.add_argument('-design', '--design', dest='design', choices=['PRIME-Del','twinPE'],help='Methods for designing the whole pegRNAs if required',required=True)
    optional.add_argument('-o', '--outdir', dest = 'outdir', type = str, help ='Path to output directory', required=False)
    optional.add_argument('-name', '--outname', dest = 'name', type = str, help ='Name of the output files except filename extensions', required=False)
    optional.add_argument('-m', '--metrics', dest = 'metrics',  nargs='+', type=str, default = defaults['metrics'], help = 'Metrics to show in the output (metrics that are used in filtering and ranking will also be shown in addition to the metrics specified by this argument)') 
    optional.add_argument('-filter', '--filter_by', dest = 'filterby', type = json.loads,  help ='Metric to use for filtering out bad spacers',default=defaults['filter_by'], required=False) 
    optional.add_argument('-erankby', '--rankby_each', dest='rankby_each', type=str, nargs='+', default=defaults['erankby'], help="Metric to use for sorting spacers, in the order of importance") #consider more
    optional.add_argument('-prankby', '--rankby_pair', dest='rankby_pair', type=json.loads, default=defaults['prankby'], help="Metric to use for sorting pairs in the order of importance") #consider more
    optional.add_argument('-pbsl', '--pbs_len', dest = 'pbs_len', type=int, default = defaults['pbsl'], help = 'PBS length in nt') 
    optional.add_argument('-rttl', '--rtt_len', dest = 'rtt_len', type=int, default = defaults['rttl'], help = 'RTT length in nt') 
    optional.add_argument('-scaf', '--scaffold', dest = 'scaf', type=str, default = defaults['scaf'], help = 'Scaffold of the prime editor (PE2 scaffold by default)') 
    optional.add_argument('-twin_rtt', '--twinPE_rtt', dest = 'twin_rtt', type=str, default = defaults['twin_rtt'], help = 'RTT sequence for Twin-PE') 

    args = parser.parse_args()

    #validating and reflecting arguments--------------------------------------------------------
    outelement = fun.get_outelement(args.name, '{args.fasta}'.format(args=args), args.design)
    outdir = fun.output_path(args.outdir, outelement)
    tmp_outdir = fun.tmp_output_path(outdir)
    logger = fun.get_log(outdir+"/"+outelement+"_single_PEEP.log", args) #making a log file
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
    forw_scored = fun.score_each(seq, forw_spacers_path_woex,ref,metrics_all, scaf, pbs_len, design, gap, sminr,  rtt_len,twin_rtt)
    logger.info("Scoring reverse spacers...")
    rev_scored = fun.score_each(seq, rev_spacers_path_woex,ref,metrics_all, scaf, pbs_len, design, gap, eminr, rtt_len,twin_rtt)

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
    ###outputs following columns: default columns, metrics-related columns, and flag-related columns (see documentation).
    ###default_columns: ["del_start", "del_end","del_length","FWD_spacer","RVS_spacer","FWD_PBS","FWD_PBS_length","RVS_PBS","RVS_PBS_length","FWD_RTT","FWD_RTT_length","RVS_RTT","RVS_RTT_length"]
    ###metrics-related columns: default or user-defined metrics for forward and reverse spacers, respectively
    ###flag-related columns: ["poly_T_4","SpacerGC_25_75","PBS_GC_30_60","PBS_RTT_GCrich"]
    forw_top = forw_ranked[0:min(forw_ranked.shape[0],top_value)]
    rev_top = rev_ranked[0:min(rev_ranked.shape[0],top_value)]
    forw_top.to_csv(outdir + "/opt_forw_" + outelement + "_single_PEEP.csv",sep=",",index=False)
    rev_top.to_csv(outdir + "/opt_rev_" + outelement + "_single_PEEP.csv",sep=",",index=False)

    logger.info("Searching for pairs with deletion length in range and assessing them...")
    pair = fun.make_pair(forw_top, rev_top, metrics_all,design,lmin,lmax,seq,pbs_len,scaf, rtt_len_arg=rtt_len,twin_rtt_arg=twin_rtt)
    logger.info(str(pair.shape[0]) + " valid pairs found with its deletion length in range...") 

    logger.info("Ranking pairs...") 
    pair_ranked = fun.rank_pair(pair, rankby_pair)

    pair_ranked.to_csv(outdir + "/" + outelement + "_single_PEEP.csv",sep=",",index=False)
    logger.info("Output files saved in ./" + outdir)
    fun.remove_tmp_outdir(tmp_outdir)

if __name__ == '__main__':
    main()

