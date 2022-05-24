from cmath import nan
from heapq import merge
import pandas as pd
import os 
os.getcwd()
os.chdir("/media/data/AtteR/projects/hiti/pipeline_output_reorg/hiti-arc-analysis")
from scripts_hiti import *
from scripts_hiti import analyze_all

os.chdir("/media/data/AtteR/projects/hiti/pipeline_output_reorg/hiti-arc-analysis/SATI_mCherry/")

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
filterlitteral = 'CTCCTCGCCCTTGCTCACCATGGTGGCGCGcctgttAACAGGCTAAGAACTCCTCTGA'
#lliteral = ' literal=TTATCCTCCTCGCCC'
export_path = '/media/data/AtteR/projects/hiti/output/'
target_sequence = "TTGCTCACCATGGTGGCGCGCCTGTTAACAGGCTAAGAACTCCTCTGAGGCAGAAGCCGGAAGGGAGCAGAGCCGGCGGCTGCAGCGGCAGTGGCAGGTGTGCGGTGCCGGAGGAGCTTAGCGAGTGTGGCAGGCTCGTCGCCGCTGAAGCTAGAGAGGCCCAGAGACTGCGGCTGCGGGAGAACTCGCTTGAGCTCTGCACCGAAACCGCCACCAGCGGCTATTTATGCTGCGCGGGGCCCGTGGGCGGCAGCTCGTGGAGCCCTGGCTCCCATTGGCGGCGCCCAGGCTCCGCGCCAGAGCCCGCCCGCGCCCTCTCCTGCCGCAGGCCCCGCCTTCCGCGCCTACTCGCTCCCCTCCCGTGGGTGCCCTCAAGGACCCGCCCCTTCACAGCCCCGAGTGACTAATGTGCTCTGCTGCGCGCCTCCCACCGGGAGGGATTTTGGGGACCGGAGACAC"

base_path = '/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SATI_5p'
export_path = '/media/data/AtteR/projects/hiti/pipeline_output_reorg/'
#target_sequence = "CCATGTTATCCTCCTCGCCCTTGCTCACCATGGTGGCGCGCCTGTTAACAGGCTAAGAACTCCTCTGAGGCAGAAGCCGGAAGGGAGCAGAGCCGGCGGCTGCAGCGGCAGTGGCAGGTGTGCGGTGCCGGAGGAGCTTAGCGAGTGTGGCAGGCTCGTCGCCGCTGAAGCTAGAGAGGCCCAGAGACTGCGGCTGCGGGAGAACTCGCTTGAGCTCTGCACCGAAACCGCCACCAGCGGCTATTTATGCTGCGCGGGGCCCGTGGGCGGCAGCTCGTGGAGCCCTGGCTCCCATTGGCGGCGCCCAGGCTCCGCGCCAGAGCCCGCCCGCGCCCTCTCCTGCCGCAGGCCCCGCCTTCCGCGCCTACTCGCTCCCCTCCCGTGGGTGCCCTCAAGGACCCGCCCCTTCACAGCCCCGAGTGACTAATGTGCTCTGCTGCGCGCCTCCCACCGGGAGGGATTTTGGGGACCGGAGACAC"
target_sequence = target_sequence.upper()
############
#read preprocessing for each sample: trim, record read counts before and after trimming, cluster the reads 
#using starcode, calculate percentage, merge into the full dataframe containing read count-and percentage for
#each sample.

df_full=import_reads_process_mini(base_path, target_sequence,filterlitteral,lliteral,rliteral,read_fwd, direc)
df_trim_full=calculate_perc_sd(df_full)
result="unaligned/Exp2_5p_mcherry_SATI.fasta"
save_fasta(result, df_trim_full, target_sequence)

csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
df_trim_full.to_csv(csv_file)

#NT
####################
output_path="aligned/NT/"
result=output_path + "Exp2_5p_mcherry_SATI_local2_prim.fasta"
aligner(df_trim_full, target_sequence, "align_local2", result, output_path, lliteral, rliteral,4,2)
####################

#AA
#484%3
####################
corr_frame=1
result="unaligned/Exp2_5p_mcherry_SATI.fasta"
output_html="aligned/AA/Exp2_5p_mcherry_SATI_AA.html"
out_csv="aligned/AA/Exp2_5p_mcherry_SATI_AA.csv"
translate_NT(result, corr_frame,direc, out_csv)

####################
