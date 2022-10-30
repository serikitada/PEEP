import functions as fun


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
    for x in range(forw_arg.shape[0]):
        for y in range(rev_arg.shape[0]):
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
                flags_part = [ polyT_flag(forw_spacer, forw_pbs, forw_rtt),
                polyT_flag(rev_spacer,rev_pbs,rev_rtt),
                spacerGC_flag(forw_spacer, spacerGC_low, spacerGC_high),
                spacerGC_flag(rev_spacer, spacerGC_low, spacerGC_high),
                pbsGC_flag(forw_pbs, pbsGC_low, pbsGC_high),
                pbsGC_flag(rev_pbs, pbsGC_low, pbsGC_high),
                pbs_rtt_GCrich_flag(forw_pbs, forw_rtt),
                pbs_rtt_GCrich_flag(rev_pbs, rev_rtt)]
                    
                pair.loc[pair.shape[0],:] = default_part + metrics_part + flags_part
    return pair