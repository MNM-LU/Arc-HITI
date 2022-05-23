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

os.chdir("/media/data/AtteR/projects/hiti/pipeline_output_reorg/hiti-arc-analysis/SP4_mCherry")

#############
#SP4 3' 
trimmed_df="/media/data/AtteR/projects/hiti/dfs/full_df_trim_mcherry_p3_HITI_SP4.csv"
transgene = 'mCherry'
read_fwd = True
direc = "3p"

base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP4_3p"

#taken from the starcode notebook
filterlitteral = 'GCGCGGGTCTTGTAGTTGCCGTCGTCCTTGAAGAAGATGGTGCGCTCCTGGACG'
lliteral = ' literal=CCTCAGAGGAGTTCT'
rliteral = ' literal=GGCGAGGAGCTGTT'


target_sequence="ctgtacaagATGGAGCTGGACCATATGACCCGGTGCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAAC"
target_sequence=target_sequence.upper()

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
#SP4 5' 
transgene = 'mCherry'
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
#gop and gep got rid of an annoying shift in NTs which made part of the alignments off so changed to 4 and
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


class translate_5p:
    def __init__(self, result, corr_frame):
        self.result=result
        self.corr_frame=corr_frame

    def translate_nt(self):
        ref_aa_cor=dict()
        aa_ampls=[]
        seq_info=[] #this will be merged later with the final dict that has the aligned AAs against the ref at certain frame
        alt_frames=[0,1,2]
        alt_frames.remove(self.corr_frame)    
        for record in SeqIO.parse(self.result, "fasta"):
            if "Ref" in record.description:
                #refs_aa_frames["Frame:" + str(alt_frame)]=str(Seq(record.seq[alt_frame:]).translate())
                rev_compl_seq=Seq(record.seq).reverse_complement()
                ref_key="Frame_corr:" + str(self.corr_frame) +"|" +str(rev_compl_seq[self.corr_frame:].translate())
            else:
                seq_info.append(record.description)
                rev_compl_seq=Seq(record.seq).reverse_complement()
                print(str(rev_compl_seq[self.corr_frame:].translate()))
                aa_ampls.append(str(rev_compl_seq[self.corr_frame:].translate()))
        ref_aa_cor[ref_key]=aa_ampls
        ref_aa=dict()
        for alt_frame in alt_frames:
            aa_ampls=[]
            for record in SeqIO.parse(self.result, "fasta"):
                if "Ref" in record.description:
                    rev_compl_seq=Seq(record.seq).reverse_complement()
                    #refs_aa_frames["Frame:" + str(alt_frame)]=str(Seq(record.seq[alt_frame:]).translate())
                    ref_key="Frame:" + str(alt_frame) +"|" +str(rev_compl_seq[alt_frame:].translate())
                else:
                    rev_compl_seq=Seq(record.seq).reverse_complement()
                    aa_ampls.append(str(rev_compl_seq[alt_frame:].translate()))
            ref_aa[ref_key]=aa_ampls
        seq_info_dic = {'Seq_stats': seq_info}
        ref_aa_cor.update(ref_aa)
        ref_aa_cor.update(seq_info_dic)
        aa_df=pd.DataFrame(ref_aa_cor)
        first_column = aa_df.pop('Seq_stats')
        # insert column using insert(position,column_name,
        # first_column) function
        aa_df.insert(0, 'Seq_stats', first_column)
        return(aa_df)

gop, gep=-3,-1
strand_dir= {"3p": translate_3p,
        "5p": translate_5p}  
from itertools import islice

trans_init = strand_dir.get(str(direc), None)  # Get the chosen class, or None if input is bad
df_aa = trans_init(result, corr_frame).translate_nt()
#df_aa.iloc[50:62,:]
len(df_aa)
for i, seq in enumerate(df_aa.iloc[:,1]):
    #print(seq)
    if seq=='':
        df_aa=df_aa.drop(df_aa.index[i]) 
len(df_aa)
df_aa.iloc[59,1]

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
for ref_fr in df_aa.columns[1:]:
    #print(ref_fr)
    seqs_in_frame=[]
    for i, aa_seq in enumerate(df_aa[ref_fr]):
        print(str(i) + "_aa_seq: " + aa_seq)
        #yield iteratively the header of the certain seq and the corresponding seq
        #header=">"+ str(id_f)+"_CluSeq:" + str((round(full_df.iloc[seq_i,-4],5))) + "_var:"+str((round(full_df.iloc[seq_i,-2],5))) +"_sd:" + str((round(full_df.iloc[seq_i,-1],5)))
        #seq_obj_align = aligner_init(full_df.iloc[seq_i,0], target_sequence, gop, gep).align()
        seq_obj_align = aligner_init(aa_seq, ref_fr.split("|")[1], gop, gep).align()
        seq_obj_align = re.sub(r'[(\d|\s]', '', seq_obj_align) #remove digits from the string caused by the alignment and empty spaces from the start
        matches=SequenceMatcher(None,ref_fr.split("|")[1],seq_obj_align)
        range_line=0
        seqs=[]
        alignments_per_ref=[]
        for i in range(len(matches.get_matching_blocks())):
            match=matches.get_matching_blocks()[i]
            seqs.append(len(ref_fr.split("|")[1][range_line:match.a])*"-"+seq_obj_align[match.b:match.b+match.size])
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

for i, ref_fr in enumerate(df_aa_align.columns[1:], start=1):
    aa_file="aligned/AA/" +result.split("/")[-1].split(".")[-2] + '_' +ref_fr.split("|")[0] + ".fasta"
    with open(aa_file, "w") as f:
        f.write(">0_Ref_" + ref_fr.split("|")[0] + "\n")
        f.write(ref_fr.split("|")[1] + "\n")
        for a, aa_seq in enumerate(list(df_aa_align.iloc[:,i])):
            f.write(">"+ df_aa_align.iloc[a,0] + "\n")
            f.write(aa_seq + "\n")
    mview_file= "aligned/AA/" +aa_file.split("/")[-1].split(".")[-2] + ".html"
    mview_command='/media/data/AtteR/Attes_bin/mview -in fasta -html head -css on -reference 1 -coloring identity ' + aa_file + '>' + mview_file
    call([mview_command], shell=True)
    print("Alignments created in html format! Files found inside aligned/AA directory")

####################
