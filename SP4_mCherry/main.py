from cmath import nan
from heapq import merge
from locale import ABDAY_2
from xml.dom.expatbuilder import theDOMImplementation
import pandas as pd
import os 
os.getcwd()
os.chdir("/media/data/AtteR/projects/hiti/pipeline_output_reorg/hiti-arc-analysis")
from scripts_hiti import *
from scripts_hiti import import_reads_process_mini
from scripts_hiti import calculate_perc_sd2


#############
#User configurable variables
#scripts_dir is the root directory of the analysis folder which contains the scripts_hiti.py script. All functions and classes
#used in the analysis have been stored in this file.
scripts_dir="/media/data/AtteR/projects/hiti/pipeline_output_reorg/hiti-arc-analysis"
os.chdir("/media/data/AtteR/projects/hiti/pipeline_output_reorg/hiti-arc-analysis")
#sample directory (inside which all the subdirectories exist)
sample_dir=scripts_dir + "/SP4_mCherry"
os.chdir(sample_dir)

#path to the analysis folder
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP4_3p"
#############

#############
#SP4 3' 
read_fwd = True
direc = "3p"
filterlitteral = 'CGGCGGCATGGACGAGCTGTACAAGATGGAGCTGGACCATATGACCCGGTGCACCGGCG'
lliteral = ' literal=CGGCGGCATGGACGAG'
rliteral = ' literal=CTGGGTCAAGCGTGAGATG'
target_sequence="ctgtacaagATGGAGCTGGACCATATGACCCGGTGCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGC"
target_sequence=target_sequence.upper()


############
#read preprocessing for each sample: trim, record read counts before and after trimming, cluster the reads 
#using starcode, calculate percentage, merge into the full dataframe containing read count-and percentage for
#each sample.
df_full=import_reads_process_mini(base_path, target_sequence,filterlitteral,lliteral,rliteral,read_fwd, direc)
df_trim_full=calculate_perc_sd(df_full)
result="unaligned/Exp2_3p_mcherry_SP4.fasta"
save_fasta(result, df_trim_full, target_sequence)

csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
df_trim_full.to_csv(csv_file)

#NT
####################
output_path="aligned/NT/"
result=output_path + "Exp2_3p_mcherry_SP4_local2_prim.fasta"

aligner(df_trim_full, target_sequence, "align_local2", result, output_path,lliteral, rliteral, 3,1)


#AA
####################
corr_frame=0
result="unaligned/Exp2_3p_mcherry_SP4.fasta"
output_html="aligned/AA/Exp2_3p_mcherry_SP4_AA.html"
out_csv="aligned/AA/Exp2_3p_mcherry_SP4_AA.csv"
aa_fasta="unaligned/Exp2_3p_mcherry_SP4_AA.fasta"
#translate and save results as csv
#df_aa=translate_nt_aa_csv(result,corr_frame, out_csv)

translate_NT(result, corr_frame,direc, out_csv)





#5'
#############
#User configurable variables

#path to the analysis folder
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP4_5p"
#############


#SP4 5'
direc = "5p"
lliteral = ' literal=CCATGTTATCCTCCTCGCCC'
rliteral = ' literal=ATTTTGGGGACCGGAGACAC'
filterlitteral="CCATGTTATCCTCCTCGCCCTTGCTCACCCGAGCTGGACCATATGACGTCATATGGT"
target_sequence="tgctcacCCGAGCTGGACCATATGACGTCATATGGTCCAGCTCCATCTGGTCGTCGGTGCTGCGGCTCCGAACAGGCTAAGAACTCCTCTGAGGCAGAAGCCGGAAGGGAGCAGAGCCGGCGGCTGCAGCGGCAGTGGCAGGTGT"
target_sequence=target_sequence.upper()
read_fwd = True

############
#read preprocessing for each sample: trim, record read counts before and after trimming, cluster the reads 
#using starcode, calculate percentage, merge into the full dataframe containing read count-and percentage for
#each sample.

df_full=import_reads_process_mini(base_path, target_sequence,filterlitteral,lliteral,rliteral,read_fwd, direc)
df_trim_full=calculate_perc_sd(df_full)
result="unaligned/Exp2_5p_mcherry_SP4.fasta"
save_fasta(result, df_trim_full, target_sequence)

csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
df_trim_full.to_csv(csv_file)

#NT
####################
output_path="aligned/NT/"
result=output_path + "Exp2_5p_mcherry_SP4_local2_prim.fasta"
aligner(df_trim_full, target_sequence, "align_local2", result, output_path, lliteral, rliteral,4,2)
####################

#AA
#484%3
####################
corr_frame=1
result="unaligned/Exp2_5p_mcherry_SP4.fasta"
output_html="aligned/AA/Exp2_5p_mcherry_SP4_AA.html"
out_csv="aligned/AA/Exp2_5p_mcherry_SP4_AA.csv"
translate_NT(result, corr_frame,direc, out_csv)

####################
