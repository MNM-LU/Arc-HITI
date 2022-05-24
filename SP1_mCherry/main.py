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
sample_dir=scripts_dir + "/SP1_mCherry"
os.chdir(sample_dir)
#path to the analysis folder
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP1_3p"
#############

#SP1 3' 
read_fwd = True
lliteral=" literal=CGGCGGCATGGACGAG"
rliteral=" literal=GTCTTCTACCGTCTGGAGAGG"
direc="3p"
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP1_3p"
filterlitteral="CGGCGGCATGGACGAGCTGTACAAGCCGAACGTTCGGAGCCGCAGCACCGACGACCAG"
target_sequence="ctgtacaagccgaacGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCT"
target_sequence=target_sequence.upper()

#########
#read preprocessing for each sample: trim, record read counts before and after trimming, cluster the reads 
#using starcode, calculate percentage, merge into the full dataframe containing read count-and percentage for
#each sample.
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral,rliteral,read_fwd, direc)
df_trim_full2=calculate_perc_sd(df_full)
result="unaligned/Exp2_3p_mcherry_SP1.fasta"
save_fasta(result, df_trim_full2, target_sequence)
####################
csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
df_trim_full2.to_csv(csv_file)

#NT
output_path="aligned/NT/"
result=output_path + "Exp2_3p_mcherry_SP1_local2_prim.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local2", result, output_path, lliteral, rliteral, 3,1)
####################
#AA
####################
corr_frame=0
result="unaligned/Exp2_3p_mcherry_SP1.fasta"
output_html="aligned/AA/Exp2_3p_mcherry_mcherry_SP1_AA.html"
out_csv="aligned/AA/Exp2_3p_mcherry_SP1_AA.csv"
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
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral, rliteral, read_fwd, direc)
df_trim_full2=calculate_perc_sd(df_full)
result="unaligned/Exp2_5p_mcherry_SP1.fasta"
save_fasta(result, df_trim_full2, target_sequence)

csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
df_trim_full2.to_csv(csv_file)

#NT
####################
output_path="aligned/NT/"
result=output_path+"Exp2_5p_mcherry_SP1_local2_prim.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local2", result, output_path,lliteral, rliteral, 4,2)
####################



#AA
#446%3
#need to make a reverse complement of an imported read?
####################
corr_frame=2
result="unaligned/Exp2_5p_mcherry_SP1.fasta"
output_html="aligned/AA/Exp2_5p_mcherry_SP1_AA.html"
out_csv="aligned/AA/Exp2_5p_mcherry_SP1_AA.csv"
translate_NT(result, corr_frame,direc, out_csv)
####################



