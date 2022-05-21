from cmath import nan
from heapq import merge
import pandas as pd
import os 
os.getcwd()
os.chdir("/media/data/AtteR/projects/hiti/pipeline_output_reorg/hiti-arc-analysis")
from scripts_hiti import *
from scripts_hiti import analyze_all
from scripts_hiti import calculate_perc_sd2
from scripts_hiti import translate_nt_aa_csv

os.chdir("/media/data/AtteR/projects/hiti/pipeline_output_reorg/hiti-arc-analysis/LP1_mCherry")


#mCherry 3p
#############
transgene = 'mCherry'
assay_end = '3p'

read_fwd = True
filterlitteral='CTGTACAAGGtcggtgctgcggctccgCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCAC' 
direc="3p"
animal_list = [7, 8, 9, 10, 11, 12]
filterlitteral = 'CTCCCTCCACACGTGCATCTCACGCTTGACCCAGCGCTCCAGGTTGGCGATGGT'
lliteral = ' literal=GGACGACGGCAACTACAAGA'
rliteral = ' literal=GGAGCTGGACCATATGACCAC' 
#filterlitteral = 'GGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGG'
#lliteral = ' literal=GTGGTCATATGGTCCAGCTCC'
#rliteral = ' literal=TCGTCCATGCCGAG'

#base_path = '/home/lcadmin/mnm-lts/SequencingFiles/arc_hiti_asyn/HITI1-8977973/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
#export_path = '/home/lcadmin/mnm-lts/HITI-analysis/'
base_path='/media/data/AtteR/projects/hiti/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
export_path='/media/data/AtteR/projects/hiti/pipeline_output_reorg/'
#original, full length
#target_sequence="CTGTACAAGGtcggtgctgcggctccgcggagccgcagcaccgacgaccagATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAG"
target_sequence = "CGGCGGCATGGACGAGCTGTACAAGGTCGGTGCTGCGGCTCCGCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAGGTCTTCTACCGTCTGGAGAGG"

#target_sequence="CTGTACAAGGTCGGTGCTGCGGCTCCGCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCACCGG"
target_sequence=target_sequence.upper()
full_df=analyze_all(base_path, transgene, filterlitteral,lliteral,rliteral,export_path,read_fwd, animal_list, target_sequence, direc)

# data_dict=create_datadict(base_path, transgene, animal_list)
# full_df=import_reads_process(data_dict, transgene,assay_end,filterlitteral,lliteral,rliteral,export_path,read_fwd)
full_df_trim=calculate_perc_sd(full_df)
result="unaligned/Exp1_3p_mcherry_LP1.fasta"
save_fasta(result, full_df_trim, target_sequence)
csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
full_df_trim.to_csv(csv_file)

#########

#NT
####################
output_path="aligned/NT/"
#result="aligned/Exp1_3p_mcherry_LP1_local3.fasta"
#test_res=aligner(full_df_trim_orig, target_sequence, "align_local3", result, output_path, -3,-1)
result=output_path + "Exp1_3p_mcherry_LP1_local2_prim.fasta"
test_res=aligner(full_df_trim, target_sequence, "align_local2", result, output_path, lliteral, rliteral, 3,1)
####################

#AA
#405%3
####################
corr_frame=0
result="unaligned/Exp1_3p_mcherry_LP1.fasta"
out_csv="aligned/AA/Exp2_3p_mcherry_LP1_AA.csv"
output_html="aligned/AA/Exp1_5p_GFP_LP1_AA.html"

translate_NT(result, corr_frame,direc, out_csv)

#df_aa=translate_nt_aa(result, 1)
#visualise_aa_hybrid_alignments(df_aa, output_html)

# df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral,read_fwd)


#mCherry 5p
############
transgene='mCherry'
assay_end='5p'
read_fwd=True
direc='5p'
#filterlitteral='CTCCCTCCACACGTGCATCTCACGCTTGACCCAGCGCTCCAGGTTGGCGATGGT' #region prior to r2 primer
# lliteral = ' literal=GTGTCTCCGGTCCCCAAAAT' 
# rliteral = ' literal=GGGCGAGGAGGATAACATGG'
filterlitteral = 'CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAG'
lliteral = ' literal=TTATCCTCCTCGCCC'
rliteral = ' literal=CCTCTGAGGCAGAA'

#base_path = '/home/lcadmin/mnm-lts/SequencingFiles/arc_hiti_asyn/HITI1-8977973/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
#export_path = '/home/lcadmin/mnm-lts/HITI-analysis/'
base_path = '/media/data/AtteR/projects/hiti/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
export_path = '/media/data/AtteR/projects/hiti/pipeline_output_reorg/'
animal_list = [1, 2, 3, 4, 5, 6]
#target_sequence = "CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAGGGGCGGGTCCTTGAGGGCACCCACGGGAGGGGAGCGAGTAGGCGCGGAAGGCGGGGCCTGCGGCAGGAGAGGGCGCGGGCGGGCTCTGGCGCGGAGCCTGGGCGCCGCCAATGGGAGCCAGGGCTCCACGAGCTGCCGCCCACGGGCCCCGCGCAGCATAAATAGCCGCTGGTGGCGGTTTCGGTGCAGAGCTCAAGCGAGTTCTCCCGCAGCCGCAGTCTCTGGGCCTCTCTAGCTTCAGCGGCGACGAGCCTGCCACACTCGCTAAGCTCCTCCGGCACCGCACACCTGCCACTGCCGCTGCAGCCGCCGGCTCTGCTCCCTTCCGGCTTCTGCCTCAGAGGAGTTCTTAGCCTGTTaacaggCGCGCCACCATGGTGAGCAAGGGCGAGGAGGATAACATGG"
# target_sequence = "CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAGGGGCGGGTCCTTGAGGGCACCCACGGGAGGGGAGCGAGTAGGCGCGGAAGGCGGGGCCTGCGGCAGGAGAGGGCGCGGGCGGGCTCTGGCGCGGAGCCTGGGCGCCGCCAATGGGAGCCAGGGCTCCACGAGCTGCCGCCCACGGGCCCCGCGCAGCATAAATAGCCGCTGGTGGCGGTTTCGGTGCAGAGCTCAAGCGAGTTCTCCCGCAGCCGCAGTCTCTGGGCCTCTCTAGCTTCAGCGGCGACGAGCCTGCCACACTCGCTAAGCTCCTCCGGCACCGCACACCTGCCACTGCCGCTGCAGCCGCCGGCTCTGCTCCCTTCCGGCTTCTGCCTCAGAGGAGTTCTTAGCCTGTTaacaggCGCGCCACCATGGTGAGCAA"
# target_sequence = "CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAGGGGCGGGTCCTTGAGGGCACCCACGGGAGGGGAGCGAGTAGGCGCGGAAGGCGGGGCCTGCGGCAGGAGAGGGCGCGGGCGGGCTCTGGCGCGGAGCCTGGGCGCCGCCAATGGGAGCCAGGGCTCCACGAGCTGCCGCCCACGGGCCCCGCGCAGCATAAATAGCCGCTGGTGGCGGTTTCGGTGCAGAGCTCAAGCGAGTTCTCCCGCAGCCGCAGTCTCTGGGCCTCTCTAGCTTCAGCGGCGACGAGCCTGCCACACTCGCTAAGCTCCTCCGGCACCGCACACCTGCCACTGCCGCTGCAGCCGCCGGCTCTGCTCCCTTCCGGCTTCTGCCTCAGAGGAGTTCTTAGCCTGTTaacaggCGCGCCACCATGGTGAGCAA"
target_sequence = "CCATGTTATCCTCCTCGCCCTTGCTCACCATGGTGGCGCGCCTGTTAACAGGCTAAGAACTCCTCTGAGGCAGAAGCCGGAAGGGAGCAGAGCCGGCGGCTGCAGCGGCAGTGGCAGGTGTGCGGTGCCGGAGGAGCTTAGCGAGTGTGGCAGGCTCGTCGCCGCTGAAGCTAGAGAGGCCCAGAGACTGCGGCTGCGGGAGAACTCGCTTGAGCTCTGCACCGAAACCGCCACCAGCGGCTATTTATGCTGCGCGGGGCCCGTGGGCGGCAGCTCGTGGAGCCCTGGCTCCCATTGGCGGCGCCCAGGCTCCGCGCCAGAGCCCGCCCGCGCCCTCTCCTGCCGCAGGCCCCGCCTTCCGCGCCTACTCGCTCCCCTCCCGTGGGTGCCCTCAAGGACCCGCCCCTTCACAGCCCCGAGTGACTAATGTGCTCTGCTGCGCGCCTCCCACCGGGAGGGATTTTGGGGACCGGAGACAC"
target_sequence = target_sequence.upper()
full_df=analyze_all(base_path, transgene, filterlitteral,lliteral,rliteral,export_path,read_fwd, animal_list, target_sequence, direc)

full_df_trim_orig=calculate_perc_sd(full_df)
result="unaligned/Exp1_5p_mcherry_LP1.fasta"
save_fasta(result, full_df_trim_orig, target_sequence)
#########
csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
full_df_trim_orig.to_csv(csv_file)

#NT
####################
output_path="aligned/NT/"
# result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_5p_mcherry_LP1_local3.fasta"
# aligner(full_df_trim_orig, target_sequence, "align_local3", result, output_path, -3,-1)
result=output_path+"Exp1_5p_mcherry_LP1_local2_prim.fasta"
aligner(full_df_trim_orig, target_sequence, "align_local2", result, output_path, lliteral, rliteral,3,1)
####################



#AA
#439%3
corr_frame=1
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_5p_mcherry_LP1.fasta"
result="unaligned/Exp1_5p_mcherry_LP1.fasta"
out_csv="aligned/AA/Exp1_5p_mcherry_LP1_AA.csv"
output_html="aligned/AA/Exp1_5p_mcherry_LP1_AA.html"

translate_NT(result, corr_frame,direc, out_csv)

# df_aa=translate_nt_aa(result, 1)
# visualise_aa_hybrid_alignments(df_aa, output_html)
# df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral,read_fwd)
#############
