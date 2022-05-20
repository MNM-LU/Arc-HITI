from cmath import nan
from heapq import merge
from xml.dom.expatbuilder import theDOMImplementation
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
direc = "3p"

base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP4_3p"
filterlitteral = 'CTCCCTCCACACGTGCATCTCACGCTTGACCCAGCGCTCCAGGTTGGCGATGGT' #region prior to r2 primer
filterlitteral='CGGCGGCATGGACGAGCTGTACAAGATGGAGCTGGACCATATGACCCGGTGCACCGGCGGC'
#filterlitteral= 'cggcggcatggacgagctgtacaagATGGAGCTGGACCATATGACCCGGTGCACCGGCGGCC'
lliteral=" literal=CGGCGGCATGGACGAG"
rliteral=" literal=CTGGGTCAAGCGTGAGATG"

#target_sequence="ctgtacaagATGGAGCTGGACCATATGACCCGGTGCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCG"
#target_sequence="ctgtacaagATGGAGCTGGACCATATGACCCGGTGCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCC"
target_sequence="ctgtacaagATGGAGCTGGACCATATGACCCGGTGCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAAC"
target_sequence=target_sequence.upper()
#prim incl.
# target_sequence="cggcggcatggacgagctgtacaagATGGAGCTGGACCATATGACCCGGTGCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCG"
# target_sequence=target_sequence.upper()
#########
#########

#The counts and percs for each animal regarding cluster seq are identical

#########
#df_full=import_reads_process_noprimtrim(base_path,target_sequence,filterlitteral,read_fwd)

df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral, rliteral, read_fwd, direc)
#starcode produces files with exactly same seqs, counts and percs
#files do have different number of reads but for some reason starcode clusters all of them as the same
df_trim_full2=calculate_perc_sd2(df_full)
#A-the original full ref between primers, B-90ish bps from leftliteral, C-no trim, target samme as in B

result="unaligned/Exp2_3p_mcherry_SP4.fasta"
save_fasta(result, df_trim_full2, target_sequence)

csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
df_trim_full2.to_csv(csv_file)

#NT
####################
output_path="aligned/NT/"
# result=output_path+"Exp2_3p_mcherry_mcherry_SP4_local3.fasta"
# aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, -3,-1)
result=output_path + "Exp2_3p_mcherry_SP4_local2_prim.fasta"
aligner(df_trim_full2, target_sequence, "align_local2", result, output_path,lliteral, rliteral, 3,1)
#AA
####################
corr_frame=0
result="unaligned/Exp2_3p_mcherry_SP4.fasta"
output_html="aligned/AA/Exp2_3p_mcherry_SP4_AA.html"
out_csv="aligned/AA/Exp2_3p_mcherry_SP4_AA.csv"
df_aa=translate_nt_aa_csv(result,corr_frame, out_csv)


translation_new(df_aa, output_html)

#AAs against the ref as values

#######################
#alignments work this way:
r="LYKMELDHMTRCTGGLHAYPAPRGGPAAKPNVILQIGKCRAEMLEHVRRTHRHLLTEVSKQVERELK"
a="LYKMELDHMTRC-G-----------------------------------------------------"
a="-Y-------TR--GGLHA-------------------------------------------------"
s="LYKMELDHMTRCG"
s="YTRGGLHA"
matches=SequenceMatcher(None, r,s)
seqs=[]
#you use range_line so that when you fill the remnants from left side of the match, you wont keep adding from
#beginning since in the end, we merge the seq list elements into the whole alignment of the amplicon against the ref
range_line=0
for i in range(len(matches.get_matching_blocks())):
    match=matches.get_matching_blocks()[i]
    aa_seq=str(s)[match.b:match.b+match.size]
    aa_seq=str(aa_seq)
    odds=["['", "[", "]","', '", ",", "',", "'", "',"]
    if aa_seq in odds:
        continue
    else:
        if "['" in seqs:
            seqs.index("[")
            seqs.replace("[","")
        seqs.append(len(r[range_line:match.a])*"-"+aa_seq)
        range_line=match.a+match.size
    #if there are empty elements, remove them as they mess up the next stage when joining the elements into a string
mer_seq=''.join(seqs)
alignments_per_ref.append(mer_seq)
ref_x_alignment[frame_ref]=alignments_per_ref

df=pd.DataFrame.from_dict(ref_x_alignment)
#######################


#########################

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
lliteral = ' literal=CCATGTTATCCTCCTCGCCC'
rliteral = ' literal=ATTTTGGGGACCGGAGACAC'
direc = "5p"
#make into reverse lit 
filterlitteral = 'CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAG'
filterlitteral='CATATGACGTCATATGGTCCAGCTCGGGTGAGCAAGGGCGAGGAGGATAACATGG'

#revcompl
filterlitteral="CCATGTTATCCTCCTCGCCCTTGCTCACCCGAGCTGGACCATATGACGTCATATGGT"
#target_sequence = "GTGTCTCCGGTCCCCAAAATCCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAGGGGCGGGTCCTTGAGGGCACCCACGGGAGGGGAGCGAGTAGGCGCGGAAGGCGGGGCCTGCGGCAGGAGAGGGCGCGGGCGGGCTCTGGCGCGGAGCCTGGGCGCCGCCAATGGGAGCCAGGGCTCCACGAGCTGCCGCCCACGGGCCCCGCGCAGCATAAATAGCCGCTGGTGGCGGTTTCGGTGCAGAGCTCAAGCGAGTTCTCCCGCAGCCGCAGTCTCTGGGCCTCTCTAGCTTCAGCGGCGACGAGCCTGCCACACTCGCTAAGCTCCTCCGGCACCGCACACCTGCCACTGCCGCTGCAGCCGCCGGCTCTGCTCCCTTCCGGCTTCTGCCTCAGAGGAGTTCTTAGCCTGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACGTCATATGGTCCAGCTCGGgtgagcaagggcgaggaggataacatgg"
target_sequence="TCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACGTCATATGGTCCAGCTCGGgtgagca"
#rev compl
target_sequence="tgctcacCCGAGCTGGACCATATGACGTCATATGGTCCAGCTCCATCTGGTCGTCGGTGCTGCGGCTCCG"
target_sequence=target_sequence.upper()
read_fwd = True
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP4_5p"
#########

df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral, rliteral, read_fwd, direc)
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

df_aa=translate_nt_aa_csv(result,corr_frame, out_csv, "5p")
translation_new(df_aa, output_html)

####################
