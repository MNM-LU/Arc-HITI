from cmath import nan
from heapq import merge
import pandas as pd
import os 

#############
#User configurable variables
#############
#scripts_dir is the root directory of the analysis folder which contains the scripts_hiti.py script. All functions and classes
#used in the analysis have been stored in this file.
scripts_dir="/media/data/AtteR/projects/hiti/pipeline_output_reorg/Arc-HITI/"
os.chdir(scripts_dir)
from scripts_hiti import *
from scripts_hiti import analyze_all
#sample directory (inside which all the subdirectories exist)
sample_dir=scripts_dir + "/HITI_mCherry_1+"
os.chdir(sample_dir)
#path to the analysis folder
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP1_3p"

#where the programs bbduk and starcode are found
#program_path="/media/data/AtteR/Attes_bin/"

#############

#SP1 3' 
read_fwd = True
lliteral=" literal=CGGCGGCATGGACGAG"
rliteral=" literal=GTCTTCTACCGTCTGGAGAGG"
direc="3p"
filterlitteral="CGGCGGCATGGACGAGCTGTACAAGCCGAACGTTCGGAGCCGCAGCACCGACGACCAG"
target_sequence="ctgtacaagccgaacGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGA"
target_sequence=target_sequence.upper()
target_sequence="CTGTACAAGCCGAAC---GTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGA"

#########
#read preprocessing for each sample: trim, record read counts before and after trimming, cluster the reads 
#using starcode, calculate percentage, merge into the full dataframe containing read count-and percentage for
#each sample. 

df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral,rliteral,read_fwd, direc)
df_trim_full=calculate_perc_sd(df_full,3)

#save_fasta(result, df_trim_full, target_sequence)
####################
result="unaligned/mCherry_3p_2+.fasta"
save_fasta(result, df_trim_full, target_sequence)

csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
csv_file="unaligned/mCherry_3p_2+.csv"
df_trim_full.to_csv(csv_file)

#full_df_trim=pd.read_csv(csv_file, index_col=[0])


###################
output_path="aligned/NT/"
result=output_path + "mCherry_3p_2+_prim.fasta"
test_res=aligner(df_trim_full, target_sequence, "align_local2", result, output_path, lliteral, rliteral, 4,2)
####################
#AA
####################
corr_frame=0
aa_primer_frame=1
result="unaligned/mCherry_3p_2+.fasta"
output_html="aligned/AA/mCherry_3p_2+.html"
out_csv="aligned/AA/mCherry_3p_2+.csv"

translate_NT(result, corr_frame,direc, out_csv, lliteral.split("=")[1], aa_primer_frame)
####################



####################
#SP1 5' 
filterlitteral = 'GCCTAGGCTAAGAACTCCTCCGCGCCACCATGGTGAGCAAGGGCGAGGAGGATAACATGG'
filterlitteral = 'CCATGTTATCCTCCTCGCCCTTGCTCACCATGGTGGCGCGGAGGAGTTCTTAGCCTAGGC'
#rev compl of prev rliteral
lliteral=' literal=CCATGTTATCCTCCTCGCCC'
rliteral = ' literal=GTGTCTCCGGTCCCCAAAAT'
read_fwd = False
direc="5p"
#rev compl
target_sequence="ttgctcaccatggtggcgcggaggagttcttagcctAGGCTAAGAACTCCTCTGAGGCAGAAGCCGGAAGGGAGCAGAGCCGGCGGCTGCAGCG"
target_sequence=target_sequence.upper()
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP1_5p"

############
#read preprocessing for each sample: trim, record read counts before and after trimming, cluster the reads 
#using starcode, calculate percentage, merge into the full dataframe containing read count-and percentage for
#each sample.
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral, rliteral, read_fwd, direc)
df_full
df_trim_full2=calculate_perc_sd(df_full,3)
df_trim_full2['percent_mean'].sum()

result="unaligned/mCherry_5p_1+.fasta"
save_fasta(result, df_trim_full2, target_sequence)

csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
df_trim_full2.to_csv(csv_file)

df_trim_full2=pd.read_csv(csv_file,index_col=[0])
df_trim_full2
#NT
####################
output_path="aligned/NT/"
result=output_path+"mCherry_5p_1+_prim.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local2", result, output_path,lliteral, rliteral, 4,2, "5p", "True")
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
