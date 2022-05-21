from cmath import nan
from heapq import merge
import pandas as pd
import os 
os.getcwd()
os.chdir("/media/data/AtteR/projects/hiti/pipeline_output_reorg/hiti-arc-analysis")
from scripts_hiti import *
from scripts_hiti import analyze_all
from scripts_hiti import calculate_perc_sd2
os.chdir("/media/data/AtteR/projects/hiti/pipeline_output_reorg/hiti-arc-analysis/LP1_GFP")

#==EXP1==#
#GFP 3p

#############
transgene = 'GFP'
assay_end = '3p'
animal_list = [20, 21, 22, 23, 24] 

read_fwd = True
# filterlitteral = 'GGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGG'
# lliteral = ' literal=GGACGACGGCAACTACAAGA' 
# rliteral = ' literal=GGAGCTGGACCATATGACCAC'
filterlitteral = 'GGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGG'
lliteral = ' literal=GTGGTCATATGGTCCAGCTCC'
rliteral = ' literal=TCGTCCATGCCGAG'
target_sequence = "TTAGCTTCTGCCTCAGAGGAGTTCTTAGCCTGTTAACAGGCGCGCCACCATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGC"

base_path = '/media/data/AtteR/projects/hiti/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
export_path = '/media/data/AtteR/projects/hiti/pipeline_output_reorg/'
target_sequence=target_sequence.upper()
direc="3p"
full_df=analyze_all(base_path, transgene, filterlitteral,lliteral,rliteral,export_path,read_fwd, animal_list, target_sequence, direc)

full_df_trim=calculate_perc_sd2(full_df)
result="unaligned/Exp1_3p_GFP_LP1.fasta"
csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
full_df_trim.to_csv(csv_file)
save_fasta(result, full_df_trim, target_sequence)
full_df_trim=pd.read_csv(csv_file,  index_col=[0])
#########
#NT
####################
output_path="aligned/NT/"
result=output_path + "Exp1_3p_GFP_LP1_local2.fasta"
aligner(full_df_trim, target_sequence, "align_local2", result, output_path, lliteral, rliteral, 3,1)
####################


#AA
####################
corr_frame=2
result="unaligned/Exp1_3p_GFP_LP1.fasta"
out_csv="aligned/AA/Exp1_3p_GFP_LP1_AA.csv"
output_html="aligned/AA/Exp1_3p_GFP_LP1_AA.html"
translate_NT(result, corr_frame,direc, out_csv)

####################




#GFP 5p
#############
transgene = 'GFP'
read_fwd = True

#(2561-2210)%3
animal_list = [13, 14, 15, 16, 17, 18] 
filterlitteral = 'ctccgCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCAC'
lliteral = ' literal=GTGGTCATATGGTCCAGCTCC'
rliteral = ' literal=GGACGACGGCAACTACAAGA'
base_path = '/media/data/AtteR/projects/hiti/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
export_path = '/media/data/AtteR/projects/hiti/pipeline_output_reorg/'
target_sequence = "GTGGTCATATGGTCCAGCTCCATCTGGTCGTCGGTGCTGCGGCTCCGcggagccgcagcaccgaCCTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGAACTCCAGCAGGACCATGTGATCGCGCT"

target_sequence=target_sequence.upper()
#########
#full_df=analyze_all()

direc="5p"
full_df=analyze_all(base_path, transgene, filterlitteral,lliteral,rliteral,export_path,read_fwd, animal_list, target_sequence, direc)

full_df_trim=calculate_perc_sd2(full_df)
result="unaligned/Exp1_5p_GFP_LP1.fasta"
csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
full_df_trim.to_csv(csv_file)
save_fasta(result, full_df_trim, target_sequence)
full_df_trim=pd.read_csv(csv_file,  index_col=[0])
#########

#NT
####################
output_path="aligned/NT/"
result=output_path+"Exp1_5p_GFP_LP1_local2.fasta"
test_res=aligner(full_df_trim, target_sequence, "align_local2", result, output_path, lliteral, rliteral,3,1)
####################

#AA
####################
corr_frame=1
result="unaligned/Exp1_5p_GFP_LP1.fasta"
out_csv="aligned/AA/Exp1_5p_GFP_LP1_AA.csv"
output_html="aligned/AA/Exp1_5p_GFP_LP1_AA.html"

translate_NT(result, corr_frame,direc, out_csv)
#############

