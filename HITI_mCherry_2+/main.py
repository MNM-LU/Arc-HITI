from cmath import nan
from heapq import merge
import pandas as pd
import os 

#############
#User configurable variables
#############
#scripts_dir is the root directory of the analysis folder which contains the scripts_hiti.py script. All functions and classes
#used in the analysis have been stored in this file.
scripts_dir="/media/data/AtteR/projects/hiti/pipeline_output_reorg/hiti-arc-analysis"
os.chdir(scripts_dir)
from scripts_hiti import *
from scripts_hiti import analyze_all
#sample directory (inside which all the subdirectories exist)
sample_dir=scripts_dir + "/HITI_mCherry_2+"
os.chdir(sample_dir)
#path to the analysis folder
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP1_3p"

#where the programs bbduk and starcode are found
program_path="/media/data/AtteR/Attes_bin"

#############

#SP1 3' 
read_fwd = True
lliteral=" literal=CGGCGGCATGGACGAG"
rliteral=" literal=GTCTTCTACCGTCTGGAGAGG"
direc="3p"
filterlitteral="CGGCGGCATGGACGAGCTGTACAAGCCGAACGTTCGGAGCCGCAGCACCGACGACCAG"
target_sequence="ctgtacaagccgaacGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCT"
target_sequence=target_sequence.upper()

#########
#read preprocessing for each sample: trim, record read counts before and after trimming, cluster the reads 
#using starcode, calculate percentage, merge into the full dataframe containing read count-and percentage for
#each sample. 

df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral,rliteral,read_fwd, direc, program_path)
df_trim_full=calculate_perc_sd(df_full)

df_trim_full.iloc[0,0]
result="unaligned/mCherry_3p_2+.fasta"
save_fasta(result, df_trim_full, target_sequence)
####################
csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
df_trim_full.to_csv(csv_file)
full_df_trim=pd.read_csv(csv_file, index_col=[0])

full_df_trim.iloc[0,0]
nt='CTGTACAAGCCGAACAGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATA'

nt='CTGTACAAGCCGAACAGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATA'

ts=(str(Seq(target_sequence[0:]).translate()))
ts
aa=(str(Seq(full_df_trim.iloc[0,0][0:]).translate()))
aa=(str(Seq('CTGTACAAGCCGAACAGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATA'[0:]).translate()))
aa
LYKPNSSEPQHRRPDGAGP
LYKPNSSEPQHRRPDGAGP
'''
NT: CTGTACAAGCCGAACAGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATA ---AA: LYKPNSSEPQHRRPDGAGP
NT: CTGTACAAGCCGAACTGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATA---AA: LYKPNCSEPQHRRPDGAGP
NT: CTGTACAAGCCGAACATGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCAT---AA: LYKPNMFGAAAPTTRWSWT
NT: CTGTACAAGCCGAACAGGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCAT---AA: LYKPNRFGAAAPTTRWSWT
'''

LYKPNVRSRSTDDQMELDHMTTGGLHAYP
LYKPN-R-R--D----------GG----P

matches=SequenceMatcher(None,ts,aa)
matches.get_matching_blocks()
range_line=0
seqs=[]
alignments_per_ref=[]
for i in range(len(matches.get_matching_blocks())):
    match=matches.get_matching_blocks()[i]
    seqs.append(len(ts[range_line:match.a])*"-"+aa[match.b:match.b+match.size])
    range_line=match.a+match.size
alignments_per_ref.append(''.join(seqs))




###################
output_path="aligned/NT/"
result=output_path + "mCherry_3p_2+_prim.fasta"
test_res=aligner(df_trim_full, target_sequence, "align_local2", result, output_path, lliteral, rliteral, 3,1)
####################
#AA
####################
corr_frame=0
result="unaligned/mCherry_3p_2+.fasta"
output_html="aligned/AA/mCherry_3p_2+.html"
out_csv="aligned/AA/mCherry_3p_2+.csv"

####################


####################
translate_NT(result, corr_frame,direc, out_csv)
####################



####################
#SP1 5' 
assay_end = '5p'
filterlitteral = 'GCCTAGGCTAAGAACTCCTCCGCGCCACCATGGTGAGCAAGGGCGAGGAGGATAACATGG'
#rev compl of prev rliteral
lliteral=' literal=CCATGTTATCCTCCTCGCCC'
rliteral = ' literal=GTGTCTCCGGTCCCCAAAAT'
read_fwd = True
direc="5p"
#rev compl
target_sequence="tgctcaccatggtggcgcggaggagttcttagcctAGGCTAAGAACTCCTCTGAGGCAGAAGCCGGAAGGGAGCAGAGCCGGCGGCTGCAGCG"
target_sequence=target_sequence.upper()
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP1_5p"

############
#read preprocessing for each sample: trim, record read counts before and after trimming, cluster the reads 
#using starcode, calculate percentage, merge into the full dataframe containing read count-and percentage for
#each sample.
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral, rliteral, read_fwd, direc, program_path)
df_trim_full2=calculate_perc_sd(df_full)
result="unaligned/mCherry_5p_2+.fasta"
save_fasta(result, df_trim_full2, target_sequence)

csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
df_trim_full2.to_csv(csv_file)

#NT
####################
output_path="aligned/NT/"
result=output_path+"mCherry_5p_2+_prim.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local2", result, output_path,lliteral, rliteral, 4,2)
####################



#AA
#446%3
#need to make a reverse complement of an imported read?
####################
# corr_frame=2
# result="unaligned/mCherry_5p_2+.fasta"
# output_html="aligned/AA/mCherry_5p_2+_AA.html"
# out_csv="aligned/AA/mCherry_5p_2+_AA.csv"
# translate_NT(result, corr_frame,direc, out_csv)
####################
