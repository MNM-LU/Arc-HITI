from cmath import nan
from heapq import merge
import pandas as pd
import os 
os.getcwd()
os.chdir("/media/data/AtteR/projects/hiti/pipeline_output_reorg/hiti-arc-analysis")
from scripts_hiti import *
from scripts_hiti import import_reads_process_mini
from scripts_hiti import calculate_perc_sd2

os.chdir("/media/data/AtteR/projects/hiti/pipeline_output_reorg/hiti-arc-analysis/SP1_mCherry")


####################
#SP1 3' 
'''
SP1
3p
ctgtacaagccgaacgttcggagccgcagcaccgacgaccagATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAG
lliteral cggcggcatggacgag rliteral GTCTTCTACCGTCTGGAGAGG

'''
trimmed_df="/media/data/AtteR/projects/hiti/dfs/full_df_trim_mcherry_p3_HITI_SP1.csv"
transgene = 'mCherry'
assay_end = '3p'
read_fwd = True
lliteral=" literal=CGGCGGCATGGACGAG"
rliteral=" literal=GTCTTCTACCGTCTGGAGAGG"

base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP1_3p"
filterlitteral = 'CTCCCTCCACACGTGCATCTCACGCTTGACCCAGCGCTCCAGGTTGGCGATGGT' #region prior to r2 primer

target_sequence="ctgtacaagccgaacgttcggagccgcagcaccgacgaccagATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAG"
target_sequence=target_sequence.upper()

#########
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral,rliteral,read_fwd)
df_trim_full2=calculate_perc_sd2(df_full)
result="unaligned/Exp2_3p_mcherry_SP1.fasta"
save_fasta(result, df_trim_full2, target_sequence)
####################
csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
df_trim_full2.to_csv(csv_file)

#NT
output_path="aligned/"
#result="aligned/Exp2_3p_mcherry_SP1_local3.fasta"
#test_res=aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, -3,-1)
result=output_path + "Exp2_3p_mcherry_SP1_local2_prim.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local2", result, output_path, lliteral, rliteral, 3,1)
####################
#AA
####################
corr_frame=0
result="unaligned/Exp2_3p_mcherry_SP1.fasta"
output_html="aligned/AA/Exp2_3p_mcherry_mcherry_SP1_AA.html"
out_csv="aligned/AA/Exp2_3p_mcherry_SP1_AA.csv"
df_aa=translate_nt_aa_csv(result,corr_frame, out_csv)
####################


####################
#SP1 5' 
'''
5p
gcagagctcaagcgagttctcccgcagccgcagtctctgggcctctctagcttcagcggcgacgagcctgccacactcgctaagctcctccggcaccgcacacctgccactgccgctgcagccgccggctctgctcccttccggcttctgcctcagaggagttcttagcctaggctaagaactcctccgcgccaccatggtgagcaa
rliteral gggcgaggaggataacatgg
'''
assay_end = '5p'
filterlitteral = 'CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAG'
rliteral = ' literal=GGGCGAGGAGGATAACATGG'
lliteral = ' literal=GTGTCTCCGGTCCCCAAAAT'


target_sequence= "CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAGGGGCGGGTCCTTGAGGGCACCCACGGGAGGGGAGCGAGTAGGCGCGGAAGGCGGGGCCTGCGGCAGGAGAGGGCGCGGGCGGGCTCTGGCGCGGAGCCTGGGCGCCGCCAATGGGAGCCAGGGCTCCACGAGCTGCCGCCCACGGGCCCCGCGCAGCATAAATAGCCGCTGGTGGCGGTTTCGGTGCAGAGCTCAAGCGAGTTCTCCCGCAGCCGCAGTCTCTGGGCCTCTCTAGCTTCAGCGGCGACGAGCCTGCCACACTCGCTAAGCTCCTCCGGCACCGCACACCTGCCACTGCCGCTGCAGCCGCCGGCTCTGCTCCCTTCCGGCTTCTGCCTCAGAGGAGTTCTTAGCCTaggctaagaactcctccgcgccaccatggtgagcaa"
target_sequence=target_sequence.upper()
read_fwd = True

base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP1_5p"
#########
#gives an error: no columns to parse from file when trying to import starcode file as a df and then calculate count and perc

df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral, rliteral, read_fwd)
df_trim_full2=calculate_perc_sd2(df_full)
result="unaligned/Exp2_5p_mcherry_SP1.fasta"
save_fasta(result, df_trim_full2, target_sequence)

csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
df_trim_full2.to_csv(csv_file)

#NT
####################
output_path="aligned/NT/"
# result=output_path + "Exp2_5p_mcherry_mcherry_SP1_local3.fasta"
# test_res=aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, -3,-1)
result=output_path+"Exp2_5p_mcherry_SP1_local2_prim.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local2", result, output_path,lliteral, rliteral, 4,2)
####################

#atm if theres a mismatch with NTs in rows colours the first NT mismatch is coloured correctly as red but the next one is not!


#AA
#446%3
####################
corr_frame=2
result="unaligned/Exp2_5p_mcherry_SP1.fasta"
output_html="aligned/AA/Exp2_5p_mcherry_SP1_AA.html"
out_csv="aligned/AA/Exp2_5p_mcherry_SP1_AA.csv"
df_aa=translate_nt_aa_csv(result,corr_frame, out_csv)
#translate_nt_aa_hiti2(result, corr_frame, output_html)
####################



