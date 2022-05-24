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
#sample directory (inside which all the subdirectories exist)
sample_dir=scripts_dir + "/LP1_mCherry/"
os.chdir(sample_dir)
#path to the analysis folder
base_path = '/media/data/AtteR/projects/hiti/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
export_path=sample_dir = "trimmed_data/"
#############


#mCherry 3p
#############
transgene="mCherry"
read_fwd = True
direc="3p"
animal_list = [7, 8, 9, 10, 11, 12]
filterlitteral = 'CTCCCTCCACACGTGCATCTCACGCTTGACCCAGCGCTCCAGGTTGGCGATGGT'
lliteral = ' literal=CGGCGGCATGGACGAG'
rliteral = ' literal=CATATGACCACCGG' 
#original, full length
target_sequence = "CTGTACAAGGtcggtgctgcggctccgCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGC"
target_sequence=target_sequence.upper()

#read preprocessing for each sample: trim, record read counts before and after trimming, cluster the reads 
#using starcode, calculate percentage, merge into the full dataframe containing read count-and percentage for
#each sample.
full_df=analyze_all(base_path, transgene, filterlitteral,lliteral,rliteral,export_path,read_fwd, animal_list, target_sequence, direc)
result="unaligned/Exp1_3p_mcherry_LP1.fasta"
full_df_trim=calculate_perc_sd(full_df)
save_fasta(result, full_df_trim, target_sequence)
csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
full_df_trim.to_csv(csv_file)
full_df_trim=pd.read_csv(csv_file, index_col=[0])

#NT
####################
output_path="aligned/NT/"
result=output_path + "Exp1_3p_mcherry_LP1_local2_prim.fasta"
test_res=aligner(full_df_trim, target_sequence, "align_local2", result, output_path, lliteral, rliteral, 3,1)
####################

#AA
#405%3
####################
corr_frame=1
result="unaligned/Exp1_3p_mcherry_LP1.fasta"
out_csv="aligned/AA/Exp2_3p_mcherry_LP1_AA.csv"
output_html="aligned/AA/Exp2_3p_mcherry_LP1_AA.html"

############
translate_NT(result, corr_frame,direc, out_csv)


#mCherry 5p
############
transgene='mCherry'
assay_end='5p'
read_fwd=True
direc='5p'
#filterlitteral = 'CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAG'
filterlitteral='CCATGTTATCCTCCTCGCCCTTGCTCACCATGGTGGCGCGcctgttAACAGGCTAAGAAC'
lliteral = ' literal=CCATGTTATCCTCCTCGCCC'

#rliteral = ' literal=CCTCTGAGGCAGAA' old from the original script
rliteral = ' literal=GTGTCTCCGGTCCCCAAAAT'
filterlitteral = 'CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAG'

#rev end
filterlitteral = 'CTCCTCGCCCTTGCTCACCATGGTGGCGCGcctgttAACAGGCTAAGAACTCCTCTGA'

#lliteral = ' literal=TTATCCTCCTCGCCC'
rliteral = ' literal=CCTCTGAGGCAGAA'
base_path = '/media/data/AtteR/projects/hiti/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
export_path = '/media/data/AtteR/projects/hiti/output/'
target_sequence = "TTGCTCACCATGGTGGCGCGCCTGTTAACAGGCTAAGAACTCCTCTGAGGCAGAAGCCGGAAGGGAGCAGAGCCGGCGGCTGCAGCGGCAGTGGCAGGTGTGCGGTGCCGGAGGAGCTTAGCGAGTGTGGCAGGCTCGTCGCCGCTGAAGCTAGAGAGGCCCAGAGACTGCGGCTGCGGGAGAACTCGCTTGAGCTCTGCACCGAAACCGCCACCAGCGGCTATTTATGCTGCGCGGGGCCCGTGGGCGGCAGCTCGTGGAGCCCTGGCTCCCATTGGCGGCGCCCAGGCTCCGCGCCAGAGCCCGCCCGCGCCCTCTCCTGCCGCAGGCCCCGCCTTCCGCGCCTACTCGCTCCCCTCCCGTGGGTGCCCTCAAGGACCCGCCCCTTCACAGCCCCGAGTGACTAATGTGCTCTGCTGCGCGCCTCCCACCGGGAGGGATTTTGGGGACCGGAGACAC"
animal_list = [1, 2, 3, 4, 5, 6]
target_sequence = target_sequence.upper()

############
#read preprocessing for each sample: trim, record read counts before and after trimming, cluster the reads 
#using starcode, calculate percentage, merge into the full dataframe containing read count-and percentage for
#each sample.

full_df=analyze_all(base_path, transgene, filterlitteral,lliteral,rliteral,export_path,read_fwd, animal_list, target_sequence, direc)
full_df_trim=calculate_perc_sd(full_df)
result="unaligned/Exp1_5p_mcherry_LP1.fasta"
save_fasta(result, full_df_trim, target_sequence)
#########
csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
full_df_trim.to_csv(csv_file)

#NT
####################
output_path="aligned/NT/"
result=output_path+"Exp1_5p_mcherry_LP1_local2_prim.fasta"
aligner(full_df_trim, target_sequence, "align_local2", result, output_path, lliteral, rliteral,3,1)
####################


#AA
#439%3
corr_frame=1
result="unaligned/Exp1_5p_mcherry_LP1.fasta"
out_csv="aligned/AA/Exp1_5p_mcherry_LP1_AA.csv"
output_html="aligned/AA/Exp1_5p_mcherry_LP1_AA.html"

translate_NT(result, corr_frame,direc, out_csv)

