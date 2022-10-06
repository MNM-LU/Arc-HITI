from cmath import nan
from heapq import merge
from locale import ABDAY_2
from xml.dom.expatbuilder import theDOMImplementation
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
#sample directory (inside which all the subdirectories exist)
sample_dir=scripts_dir + "/HITI_mCherry_9+"
os.chdir(sample_dir)
#path to the analysis folder
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP4_3p"

#where the programs bbduk and starcode are found
program_path="/media/data/AtteR/Attes_bin"

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
df_trim_full=calculate_perc_sd(df_full,3)
result="unaligned/HITI_mCherry_3p_9+.fasta"
save_fasta(result, df_trim_full, target_sequence)

csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
df_trim_full.to_csv(csv_file)
df_trim_full=pd.read_csv(csv_file, index_col=[0])
#NT
####################
output_path="aligned/NT/"
result=output_path + "HITI_mCherry_3p_9+_prim.fasta"
aligner(df_trim_full, target_sequence, "align_local2", result, output_path,lliteral, rliteral, 3,1)


#AA
####################
corr_frame=0
aa_primer_frame=1

result="unaligned/HITI_mCherry_3p_9+.fasta"
output_html="aligned/AA/HITI_mCherry_3p_9+_AA.html"
out_csv="aligned/AA/HITI_mCherry_3p_9+_AA.csv"
aa_fasta="unaligned/HITI_mCherry_3p_9+_AA.fasta"
#translate and save results as csv
#df_aa=translate_nt_aa_csv(result,corr_frame, out_csv)

translate_NT(result, corr_frame,direc, out_csv, lliteral.split("=")[1], aa_primer_frame)











#5'
#############
#User configurable variables

#path to the analysis folder
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP4_5p"
#############


#SP4 5'
direc = "5p"
lliteral = ' literal=CCATGTTATCCTCCTCGCC'  #trimmed away one extra T from the amplicons
rliteral = ' literal=ATTTTGGGGACCGGAGACAC'

filterlitteral="CCATGTTATCCTCCTCGCCCTTGCTCACCCGAGCTGGACCATATGACGTCATATGGT"
#target_sequence="tgctcacCCGAGCTGGACCATATGACGTCATATGGTCCAGCTCCAT"
#target_sequence="tgctcacCCGAGCTGGACCATATGACGTCATATGGTCCAGCTCCATCTGGTCGTCGGTG"
target_sequence="tgctcacCCGAGCTGGACCATATGACGTCATATGGTCCAGCTCCATCTGGTCGTCGGTGCTGCGGCTCCGAACAGGCTAAGAACTCCTCTGAGGCAGAAGCCGGAAGGGAGCAGAGCCGGCGGCTGCAGCGGCAGTGGCAGGTGTGCGGTGCCGGAGGAGCTTAGCGAGTGTGGCAGGCTCGT"
target_sequence = "cttgctcacCCGAGCTGGACCATATGACGTCATATGGTCCAGCTCCATCTGGTCGTCGGTGCTGCGGCTCCGAACAGGCTAAGAACTCCTCTGAGGCAGAAGCCGGAAGGGAGCAGAGCCGGCGGCTGCAGCGGCAGTGGCAGGTGTGCGGTGCCGGAGGAGCTTAGCGAGTGTGGCAGGCTCGT"
target_sequence = "cttgctcacCCGAGCTGGACCATATGACGTCATATGGTCCAGCTCCATCTG"
target_sequence=target_sequence.upper()
read_fwd = True
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP4_5p"


############
#read preprocessing for each sample: trim, record read counts before and after trimming, cluster the reads 
#using starcode, calculate percentage, merge into the full dataframe containing read count-and percentage for
#each sample.

df_full=import_reads_process_mini(base_path, target_sequence,filterlitteral,lliteral,rliteral,read_fwd, direc)
df_trim_full=calculate_perc_sd(df_full,3)
result="unaligned/HITI_mCherry_5p_9+.fasta"
save_fasta(result, df_trim_full, target_sequence)

csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
df_trim_full.to_csv(csv_file)
df_trim_full=pd.read_csv(csv_file,  index_col=[0])

#NT
####################
output_path="aligned/NT/"
result=output_path + "HITI_mCherry_5p_9+_prim.fasta"

aligner(df_trim_full, target_sequence, "align_local2", result, output_path, lliteral, rliteral,4,2,'5p', "True")
####################

#AA
#484%3
####################
corr_frame=0
aa_primer_frame=0

result="unaligned/HITI_mCherry_5p_9+.fasta"
output_html="aligned/AA/HITI_mCherry_5p_9+_AA.html"
out_csv="aligned/AA/HITI_mCherry_5p_9+_AA.csv"

direc="5p"
translate_NT(result, corr_frame,direc, out_csv, lliteral.split("=")[1],aa_primer_frame)



##################################################################################################
# testing stuff

s = Seq("CTTGCTCACCCGAGCTGGACCATATGACGTCATATGGTCCAGCTCCATCTGGTCGT").reverse_complement()
s[2:].translate()

full_df_trim_AA=pd.read_csv(out_csv, index_col=[0])

strand_dir= {"3p": translate_3p,
        "5p": translate_5p}  

trans_init = strand_dir.get(str(dir), None)  # Get the chosen class, or None if input is bad
df_aa = trans_init(result, corr_frame).translate_nt()
for i, seq in enumerate(df_aa.iloc[:,1]):
    if seq=='':
        try:
            df_aa=df_aa.drop(df_aa.index[i]) 
        except IndexError: #indexing error that occured with gfp 3p
            continue
align_method="align_local3"
align_class = {"align_local": align_local,
        "align_local2": align_local2,
        "align_local3":align_local3,
        "align_global":align_global,
        "align_global2":align_global2}  
id_f=1
aligner_init = align_class.get(str(align_method), None)  # Get the chosen class, or None if input is bad
aa_data_align=pd.DataFrame(columns = list(df_aa.columns))
#aa_data_align=aa_data_align.append(df_aa.iloc[:,0])
aa_data_align['Seq_stats']=df_aa.iloc[:,0]
seq_info_dic=df_aa.iloc[:,0]
seq_info_dic={df_aa.columns[0]: df_aa.iloc[:,0]}
dic_aa_align=dict()
print(df_aa)
for ref_fr in df_aa.columns[1:]:
    #print(ref_fr)
    seqs_in_frame=[]
    for i, aa_seq in enumerate(df_aa[ref_fr]):
        #yield iteratively the header of the certain seq and the corresponding seq
        #header=">"+ str(id_f)+"_CluSeq:" + str((round(full_df.iloc[seq_i,-4],5))) + "_var:"+str((round(full_df.iloc[seq_i,-2],5))) +"_sd:" + str((round(full_df.iloc[seq_i,-1],5)))
        
        #seq_obj_align = aligner_init(aa_seq, ref_fr.split("|")[1], gop, gep).align()
        #seq_obj_align = re.sub(r'[(\d|\s]', '', seq_obj_align) #remove digits from the string caused by the alignment and empty spaces from the start
        matches=SequenceMatcher(None,ref_fr.split("|")[1],aa_seq)
        matches.get_matching_blocks()
        range_line=0
        seqs=[]
        alignments_per_ref=[]
        for i in range(len(matches.get_matching_blocks())):
            match=matches.get_matching_blocks()[i]
            seqs.append(len(ref_fr.split("|")[1][range_line:match.a])*"-"+aa_seq[match.b:match.b+match.size])
            range_line=match.a+match.size
        alignments_per_ref.append(''.join(seqs))
        alignments_per_ref= str(alignments_per_ref).replace("['", "").replace("']", "")
        seqs_in_frame.append(str(alignments_per_ref))
    dic_aa_align[ref_fr] = seqs_in_frame
dic_aa_align.update(seq_info_dic)
dic_aa_align.keys()
df_aa_align=pd.DataFrame(dic_aa_align)
first_column = df_aa_align.pop('Seq_stats')
# insert column using insert(position,column_name,
# first_column) function
df_aa_align.insert(0, 'Seq_stats', first_column)
df_aa_align.to_csv(out_csv)
#write into file, add translated primer

p = "/media/data/AtteR/projects/hiti/pipeline_output_reorg/hiti-arc-analysis/HITI_mCherry_9+/aligned/AA/HITI_mCherry_5p_9+_AA.csv"
df_aa_align = pd.read_csv(p, index_col=[0])
df_aa_align.columns
df_aa_align.iloc[0,2]
primer = lliteral.split("=")[1]
from subprocess import call
from Bio.Seq import Seq

aa_file = "/media/data/AtteR/projects/hiti/pipeline_output_reorg/hiti-arc-analysis/HITI_mCherry_9+/aligned/AA/fasta/HITI_mCherry_5p_9+_Frame_corr_2.fasta"
aa_file = "/media/data/AtteR/projects/hiti/pipeline_output_reorg/hiti-arc-analysis/HITI_mCherry_9+/unaligned/HITI_mCherry_5p_9+.fasta"
mview_file= "aligned/AA/html/" +aa_file.split("/")[-1].split(".")[-2] + ".html"
mview_file
mview_command='/media/data/AtteR/Attes_bin/mview -in fasta -html head -css on -reference 1 -coloring identity ' + aa_file + '>' + mview_file

# mview_command='/media/data/AtteR/Attes_bin/mview -in fasta -html head -css on -reference 1 -coloring identity ' + aa_file + '>' + mview_file

call([mview_command], shell=True)
for a, aa_seq in enumerate(list(df_aa_align.iloc[:,i])):
    if len(aa_seq)%3!=0:
        len(aa_seq)%3 * - + aa_seq

for i, ref_fr in enumerate(df_aa_align.columns[1:], start=1):
    aa_file="aligned/AA/fasta/" +result.split("/")[-1].split(".")[-2] + '_' +ref_fr.split("|")[0] + ".fasta"
    with open(aa_file, "w") as f:
        f.write(">0_Ref_frame_" + str(corr_frame) + "\n")
        if len(ref_fr.split("|")[1])%3!=0:
            f.write(len(ref_fr.split("|")[1])%3 * '-' + ref_fr.split("|")[1] + "-" + str(Seq(primer).reverse_complement()[aa_primer_frame:].translate()) + "\n")
        else:
            f.write(ref_fr.split("|")[1] + "-" + str(Seq(primer).reverse_complement()[aa_primer_frame:].translate()) + "\n")
        for a, aa_seq in enumerate(list(df_aa_align.iloc[:,i])):
            f.write(">"+ df_aa_align.iloc[a,0] + "\n")
            if len(aa_seq)%3 != 0:
                print(f'unequal by {len(aa_seq)%3}')
                f.write(len(aa_seq)%3 * "-" + aa_seq + "-" + str(Seq(primer).reverse_complement()[aa_primer_frame:].translate()) + "\n")
            else:
                f.write(aa_seq + "-" + str(Seq(primer).reverse_complement()[aa_primer_frame:].translate()) + "\n")

    mview_file= "aligned/AA/html/" +aa_file.split("/")[-1].split(".")[-2] + ".html"
    mview_command='/media/data/AtteR/Attes_bin/mview -in fasta -html head -css on -reference 1 -coloring identity ' + aa_file + '>' + mview_file
    call([mview_command], shell=True)
    print("Alignments created in html format! Files found inside aligned/AA directory")


len("EPATLAKLLRHRTPATAAAAAGSAPFRLLPQRSS*PVRSRSTDDQMELDHMTSYGPARVSK-GEEDNM")
len("--------------------------------------------DQMELDHMTSYGPARVSK-GEEDNM")
if "3" in dir:
    for i, ref_fr in enumerate(df_aa_align.columns[1:], start=1):
        aa_file="aligned/AA/fasta/" +result.split("/")[-1].split(".")[-2] + '_' +ref_fr.split("|")[0] + ".fasta"
        with open(aa_file, "w") as f:
            f.write(">0_Ref_" + str(corr_frame) + "\n")
            f.write(str(Seq(primer[aa_primer_frame:]).translate()) + "-" + ref_fr.split("|")[1] + "\n")
            for a, aa_seq in enumerate(list(df_aa_align.iloc[:,i])):
                f.write(">"+ df_aa_align.iloc[a,0] + "\n")
                f.write(str(Seq(primer[aa_primer_frame:]).translate()) + "-" + aa_seq + "\n")
        mview_file= "aligned/AA/html/" +aa_file.split("/")[-1].split(".")[-2] + ".html"
        mview_command='/media/data/AtteR/Attes_bin/mview -in fasta -html head -css on -reference 1 -coloring identity ' + aa_file + '>' + mview_file
        call([mview_command], shell=True)
        print("Alignments created in html format! Files found inside aligned/AA directory")

else:
    for i, ref_fr in enumerate(df_aa_align.columns[1:], start=1):
        aa_file="aligned/AA/fasta/" +result.split("/")[-1].split(".")[-2] + '_' +ref_fr.split("|")[0] + ".fasta"
        with open(aa_file, "w") as f:
            f.write(">0_Ref_" + str(corr_frame) + "\n")
            f.write(ref_fr.split("|")[1] + "-" + str(Seq(primer).reverse_complement()[aa_primer_frame:].translate()) + "\n")
            for a, aa_seq in enumerate(list(df_aa_align.iloc[:,i])):
                f.write(">"+ df_aa_align.iloc[a,0] + "\n")
                f.write(aa_seq + "-" + str(Seq(primer).reverse_complement()[aa_primer_frame:].translate()) + "\n")
        mview_file= "aligned/AA/html/" +aa_file.split("/")[-1].split(".")[-2] + ".html"
        mview_command='/media/data/AtteR/Attes_bin/mview -in fasta -html head -css on -reference 1 -coloring identity ' + aa_file + '>' + mview_file
        call([mview_command], shell=True)
        print("Alignments created in html format! Files found inside aligned/AA directory")


####################
