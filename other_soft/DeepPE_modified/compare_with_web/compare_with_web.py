import os,sys
import pandas as pd
FILE_DIRECTORY = os.path.dirname(os.path.realpath(__file__))
sys.path.append(FILE_DIRECTORY)
parent = os.path.dirname(FILE_DIRECTORY)
sys.path.append(parent)
import DeepPE_main

web = pd.read_csv(FILE_DIRECTORY+"/web_example_DeepPE_Scores.txt",sep="\t")
score = []
for row in range(web.shape[0]): 
    if row%10==0: print(row)
    seq_47nt = web.iloc[row,2]
    pbs_len = web.iloc[row,5]
    rtt_len = web.iloc[row,6]
    rtt = str(web.iloc[row,7])[0:rtt_len]
    pbs = str(web.iloc[row,7])[rtt_len:(rtt_len+pbs_len)]
    #print(seq_47nt,pbs_len,rtt_len,rtt,pbs)
    score.append(DeepPE_main.main(seq47_arg=seq_47nt, pbs_arg=pbs, rtt_arg=rtt))

df = pd.DataFrame(score)
df.to_csv(FILE_DIRECTORY+"/command_Scores_woGtoC.txt", sep="\t",index=False)

