from cmath import nan
from heapq import merge
import pandas as pd
import os 
os.getcwd()
os.chdir("/media/data/AtteR/projects/hiti/pipeline_output_reorg/hiti-arc-analysis")
from scripts_hiti import *
from scripts_hiti import import_reads_process_mini
from scripts_hiti import calculate_perc_sd2

os.chdir("/media/data/AtteR/projects/hiti/pipeline_output_reorg/hiti-arc-analysis/SP4_mCherry")

#############
#SP4 3' 
trimmed_df="/media/data/AtteR/projects/hiti/dfs/full_df_trim_mcherry_p3_HITI_SP4.csv"
transgene = 'mCherry'
assay_end = '3p'
read_fwd = True
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP4_3p"
filterlitteral = 'CTCCCTCCACACGTGCATCTCACGCTTGACCCAGCGCTCCAGGTTGGCGATGGT' #region prior to r2 primer
lliteral=" literal=CGGCGGCATGGACGAG"
rliteral=" literal=CTGGGTCAAGCGTGAGATG"

target_sequence="ctgtacaagATGGAGCTGGACCATATGACCCGGTGCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCG"
target_sequence=target_sequence.upper()
#########
#########

#The counts and percs for each animal regarding cluster seq are identical

#########
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral, rliteral, read_fwd)
#starcode produces files with exactly same seqs, counts and percs
#files do have different number of reads but for some reason starcode clusters all of them as the same
df_trim_full2=calculate_perc_sd2(df_full)
result="unaligned/Exp2_3p_mcherry_SP4.fasta"
save_fasta(result, df_trim_full2, target_sequence)

csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
df_trim_full2.to_csv(csv_file)

#NT
####################
output_path="aligned/NT/"
# result=output_path+"Exp2_3p_mcherry_mcherry_SP4_local3.fasta"
# aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, -3,-1)
result=output_path + "Exp2_3p_mcherry_mcherry_SP4_local2_prim.fasta"
aligner(df_trim_full2, target_sequence, "align_local2", result, output_path,lliteral, rliteral, 3,1)
#AA
####################
corr_frame=0
result="unaligned/Exp2_3p_mcherry_SP4.fasta"
output_html="aligned/AA/Exp2_3p_mcherry_SP4_AA.html"
out_csv="aligned/AA/Exp2_3p_mcherry_SP4_AA.csv"
df_aa=translate_nt_aa_csv(result,corr_frame, out_csv)
####################
#take each of the columns, align them and save into output files

#5'
#############
#SP4 5' 
'''
SP4
5p
gcagagctcaagcgagttctcccgcagccgcagtctctgggcctctctagcttcagcggcgacgagcctgccacactcgctaagctcctccggcaccgcacacctgccactgccgctgcagccgccggctctgctcccttccggcttctgcctcagaggagttcttagcctgttcggagccgcagcaccgacgaccagATGGAGCTGGACCATATGACGTCATATGGTCCAGCTCGGgtgagcaa
lliteral 
rliteral gggcgaggaggataacatgg

'''
transgene = 'mCherry'
assay_end = '5p'
lliteral = ' literal=GTGTCTCCGGTCCCCAAAAT'
rliteral = ' literal=GGGCGAGGAGGATAACATGG'
filterlitteral = 'CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAG'
filterlitteral='AGCCTGTTAACAGGCGCGCCACCATGGTGAGCAAGGGCGAGGAGGATAACATGG'
target_sequence = "GTGTCTCCGGTCCCCAAAATCCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAGGGGCGGGTCCTTGAGGGCACCCACGGGAGGGGAGCGAGTAGGCGCGGAAGGCGGGGCCTGCGGCAGGAGAGGGCGCGGGCGGGCTCTGGCGCGGAGCCTGGGCGCCGCCAATGGGAGCCAGGGCTCCACGAGCTGCCGCCCACGGGCCCCGCGCAGCATAAATAGCCGCTGGTGGCGGTTTCGGTGCAGAGCTCAAGCGAGTTCTCCCGCAGCCGCAGTCTCTGGGCCTCTCTAGCTTCAGCGGCGACGAGCCTGCCACACTCGCTAAGCTCCTCCGGCACCGCACACCTGCCACTGCCGCTGCAGCCGCCGGCTCTGCTCCCTTCCGGCTTCTGCCTCAGAGGAGTTCTTAGCCTGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACGTCATATGGTCCAGCTCGGgtgagcaagggcgaggaggataacatgg"
target_sequence="CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAGGGGCGGGTCCTTGAGGGCACCCACGGGAGGGGAGCGAGTAGGCGCGGAAGGCGGGGCCTGCGGCAGGAGAGGGCGCGGGCGGGCTCTGGCGCGGAGCCTGGGCGCCGCCAATGGGAGCCAGGGCTCCACGAGCTGCCGCCCACGGGCCCCGCGCAGCATAAATAGCCGCTGGTGGCGGTTTCGGTGCAGAGCTCAAGCGAGTTCTCCCGCAGCCGCAGTCTCTGGGCCTCTCTAGCTTCAGCGGCGACGAGCCTGCCACACTCGCTAAGCTCCTCCGGCACCGCACACCTGCCACTGCCGCTGCAGCCGCCGGCTCTGCTCCCTTCCGGCTTCTGCCTCAGAGGAGTTCTTAGCCTGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACGTCATATGGTCCAGCTCGGgtgagcaa"
target_sequence=target_sequence.upper()
read_fwd = True
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP4_5p"
#########
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral, rliteral, read_fwd)
df_trim_full2=calculate_perc_sd2(df_full)
result="unaligned/Exp2_5p_mcherry_SP4.fasta"
save_fasta(result, df_trim_full2, target_sequence)

csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
df_trim_full2.to_csv(csv_file)

#NT
####################
output_path="aligned/NT/"
#result="aligned/Exp2_5p_mcherry_mcherry_SP4_local3.fasta"
# test_res=aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, -3,-1)
result=output_path + "Exp2_5p_mcherry_SP4_local2_prim.fasta"
#gop and gep got rid of an annoying shift in NTs which made part of the alignments off so changed to 4 and
aligner(df_trim_full2, target_sequence, "align_local2", result, output_path, lliteral, rliteral,4,2)
####################

#AA
#484%3
####################
corr_frame=1
result="unaligned/Exp2_5p_mcherry_SP4.fasta"
output_html="aligned/AA/Exp2_5p_mcherry_SP4_AA.html"
out_csv="aligned/AA/Exp2_5p_mcherry_SP4_AA.csv"
df_aa=translate_nt_aa_csv(result,corr_frame, out_csv)
####################
