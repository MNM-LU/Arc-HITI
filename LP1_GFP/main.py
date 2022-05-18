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
filterlitteral = 'GGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGG'
lliteral = ' literal=GGACGACGGCAACTACAAGA'
rliteral = ' literal=GGAGCTGGACCATATGACCAC'
base_path = '/media/data/AtteR/projects/hiti/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
export_path = '/media/data/AtteR/projects/hiti/pipeline_output_reorg/'
target_sequence = "CCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGGtcggtgctgcggctccgcggagccgcagcaccgacgaccagAT"
target_sequence=target_sequence.upper()
#full_df=analyze_all()
full_df=analyze_all(base_path,transgene,assay_end,filterlitteral,lliteral,rliteral,export_path,read_fwd, animal_list, target_sequence)

# data_dict=create_datadict(base_path,transgene, animal_list)
# full_df=import_reads_process(data_dict, transgene,assay_end,filterlitteral,lliteral,rliteral,export_path,read_fwd)
full_df_trim=calculate_perc_sd2(full_df)
result="unaligned/Exp1_3p_GFP_LP1_TB.fasta"
csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
full_df_trim.to_csv(csv_file)
save_fasta(result, full_df_trim, target_sequence)
full_df_trim=pd.read_csv(csv_file,  index_col=[0])
#########
#NT
####################
output_path="aligned/"
#result="aligned/Exp1_3p_GFP_LP1_local3.fasta"
#test_res=aligner(full_df_trim, target_sequence, "align_local3", result, output_path, -3,-1)
result=output_path + "Exp1_3p_GFP_LP1_local2_TB.fasta"
aligner(full_df_trim, target_sequence, "align_local2", result, output_path, lliteral, rliteral, 3,1)
####################


#AA
####################
# corr_frame=1
# result="aligned/Exp1_3p_GFP_LP1.fasta"
# df_aa=translate_nt_aa(result, 1)
# output_html="aligned/Exp1_3p_GFP_LP1_AA.html"
# visualise_aa_hybrid_alignments(df_aa, output_html)

# df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral,read_fwd)
####################




#GFP 5p
#############
transgene = 'GFP'
assay_end = '5p'
read_fwd = True
animal_list = [13, 14, 15, 16, 17, 18] 
filterlitteral = 'GCGCGGGTCTTGTAGTTGCCGTCGTCCTTGAAGAAGATGGTGCGCTCCTGGACG'
lliteral = ' literal=GCTTCTGCCTCAGAGGAGTTCT'
rliteral = ' literal=CGAGGTGAAGTTCGAGGGC'
base_path = '/media/data/AtteR/projects/hiti/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
export_path = '/media/data/AtteR/projects/hiti/pipeline_output_reorg/'
target_sequence = "tagcctgttaacaggCGCGCCACCATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGC"
target_sequence=target_sequence.upper()
#########
full_df=analyze_all()

full_df=analyze_all(base_path,transgene,assay_end,filterlitteral,lliteral,rliteral,export_path,read_fwd, animal_list, target_sequence)

# data_dict=create_datadict(base_path,transgene, animal_list)
# full_df=import_reads_process(data_dict, transgene,assay_end,filterlitteral,lliteral,rliteral,export_path,read_fwd)
full_df_trim=calculate_perc_sd2(full_df)
result="unaligned/Exp1_5p_GFP_LP1_TB.fasta"
csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
full_df_trim.to_csv(csv_file)
save_fasta(result, full_df_trim, target_sequence)
full_df_trim=pd.read_csv(csv_file,  index_col=[0])
#########

#NT
####################
output_path="aligned/"
#result="aligned/Exp1_5p_GFP_LP1_local3.fasta"
#test_res=aligner(full_df_trim, target_sequence, "align_local3", result, output_path, -3,-1)
result=output_path+"Exp1_5p_GFP_LP1_local2_TB.fasta"
test_res=aligner(full_df_trim, target_sequence, "align_local2", result, output_path, lliteral, rliteral,3,1)
####################

#AA
####################
# corr_frame=1
# result="aligned/Exp1_5p_GFP_LP1.fasta"
# df_aa=translate_nt_aa(result, 1)
# output_html="aligned/Exp1_5p_GFP_LP1_AA.html"
# visualise_aa_hybrid_alignments(df_aa, output_html)
# df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral,read_fwd)
#############

