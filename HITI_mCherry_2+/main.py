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
from scripts_hiti import analyze_all
#sample directory (inside which all the subdirectories exist)
sample_dir=scripts_dir + "/HITI_mCherry_2+"
os.chdir(sample_dir)
#path to the analysis folder
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP1_3p"
#############

#SP1 3' 
read_fwd = True
lliteral=" literal=CGGCGGCATGGACGAG"
rliteral=" literal=GTCTTCTACCGTCTGGAGAGG"
direc="3p"
filterlitteral="CGGCGGCATGGACGAGCTGTACAAGCCGAACGTTCGGAGCCGCAGCACCGACGACCAG"
target_sequence="ctgtacaagccgaacGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCT"
target_sequence=target_sequence.upper()

#########
#read preprocessing for each sample: trim, record read counts before and after trimming, cluster the reads 
#using starcode, calculate percentage, merge into the full dataframe containing read count-and percentage for
#each sample. 

df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral,rliteral,read_fwd, direc)
df_trim_full=calculate_perc_sd(df_full)

df_trim_full.iloc[0,0]
result="unaligned/mCherry_3p_2+.fasta"
save_fasta(result, df_trim_full, target_sequence)
####################
csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
df_trim_full.to_csv(csv_file)
full_df_trim=pd.read_csv(csv_file, index_col=[0])

full_df_trim.iloc[0,0]
nt='CTGTACAAGCCGAACAGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATA'

nt='CTGTACAAGCCGAACAGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATA'

ts=(str(Seq(target_sequence[0:]).translate()))
ts
aa=(str(Seq(full_df_trim.iloc[0,0][0:]).translate()))
aa=(str(Seq('CTGTACAAGCCGAACAGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATA'[0:]).translate()))
aa
LYKPNSSEPQHRRPDGAGP
LYKPNSSEPQHRRPDGAGP
'''
NT: CTGTACAAGCCGAACAGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATA ---AA: LYKPNSSEPQHRRPDGAGP
NT: CTGTACAAGCCGAACTGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATA---AA: LYKPNCSEPQHRRPDGAGP
NT: CTGTACAAGCCGAACATGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCAT---AA: LYKPNMFGAAAPTTRWSWT
NT: CTGTACAAGCCGAACAGGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCAT---AA: LYKPNRFGAAAPTTRWSWT
'''

LYKPNVRSRSTDDQMELDHMTTGGLHAYP
LYKPN-R-R--D----------GG----P

matches=SequenceMatcher(None,ts,aa)
matches.get_matching_blocks()
range_line=0
seqs=[]
alignments_per_ref=[]
for i in range(len(matches.get_matching_blocks())):
    match=matches.get_matching_blocks()[i]
    seqs.append(len(ts[range_line:match.a])*"-"+aa[match.b:match.b+match.size])
    range_line=match.a+match.size
alignments_per_ref.append(''.join(seqs))




###################
output_path="aligned/NT/"
result=output_path + "mCherry_3p_2+_prim.fasta"
test_res=aligner(df_trim_full, target_sequence, "align_local2", result, output_path, lliteral, rliteral, 3,1)
####################
#AA
####################
corr_frame=0
result="unaligned/mCherry_3p_2+.fasta"
output_html="aligned/AA/mCherry_3p_2+.html"
out_csv="aligned/AA/mCherry_3p_2+.csv"
ref_aa_cor=dict()
aa_ampls=[]
seq_info=[] #this will be merged later with the final dict that has the aligned AAs against the ref at certain frame

alt_frames=[0,1,2]
alt_frames.remove(corr_frame)    
for record in SeqIO.parse(result, "fasta"):
    if "Ref" in record.description:
        #refs_aa_frames["Frame:" + str(alt_frame)]=str(Seq(record.seq[alt_frame:]).translate())
        ref_key="Frame_corr:" + str(corr_frame) +"|" +str(Seq(record.seq[corr_frame:]).translate())
    else:
        seq_info.append(record.description)
        aa_ampls.append(str(Seq(record.seq[corr_frame:]).translate()))
        print("NT: " + record.seq + "---AA: "+ str(Seq(record.seq[corr_frame:]).translate()))


aa_ampls[0]
ref_aa_cor[ref_key]=aa_ampls
ref_aa_cor.keys()
ref_aa_cor['Frame_corr:0|LYKPNVRSRSTDDQMELDHMTTGGLHAYP'][0]
ref_aa=dict()
for alt_frame in alt_frames:
    aa_ampls=[]
    for record in SeqIO.parse(result, "fasta"):
        if "Ref" in record.description:
            #refs_aa_frames["Frame:" + str(alt_frame)]=str(Seq(record.seq[alt_frame:]).translate())
            ref_key="Frame:" + str(alt_frame) +"|" +str(Seq(record.seq[alt_frame:]).translate())
        else:
            aa_ampls.append(str(Seq(record.seq[alt_frame:]).translate()))
    ref_aa[ref_key]=aa_ampls
seq_info_dic = {'Seq_stats': seq_info}
ref_aa_cor.update(ref_aa)
ref_aa_cor.update(seq_info_dic)
aa_df=pd.DataFrame(ref_aa_cor)
first_column = aa_df.pop('Seq_stats')
# insert column using insert(position,column_name,
# first_column) function
aa_df.insert(0, 'Seq_stats', first_column)

aa_df.iloc[0,1]
##############################################################
##############################################################
class translate_3p:
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
                ref_key="Frame_corr:" + str(self.corr_frame) +"|" +str(Seq(record.seq[self.corr_frame:]).translate())
            else:
                seq_info.append(record.description)
                aa_ampls.append(str(Seq(record.seq[self.corr_frame:]).translate()))
        ref_aa_cor[ref_key]=aa_ampls
        ref_aa=dict()
        for alt_frame in alt_frames:
            aa_ampls=[]
            for record in SeqIO.parse(self.result, "fasta"):
                if "Ref" in record.description:
                    #refs_aa_frames["Frame:" + str(alt_frame)]=str(Seq(record.seq[alt_frame:]).translate())
                    ref_key="Frame:" + str(alt_frame) +"|" +str(Seq(record.seq[alt_frame:]).translate())
                else:
                    aa_ampls.append(str(Seq(record.seq[alt_frame:]).translate()))
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
class align_local2():
    aligned_data=dict()
    def __init__(self, amplicon,ref,gop=-3, gep=-1):
        self.amplicon=amplicon
        self.ref=ref
        self.gop=gop
        self.gep=gep

    def align(self):
        alignment, score, start_end_positions = local_pairwise_align_ssw(DNA(self.ref),DNA(self.amplicon),gap_open_penalty=self.gop,gap_extend_penalty = self.gep)
        out_align = ('-'*start_end_positions[0][0])+str(alignment[1])+('-'*(len(self.ref)-start_end_positions[0][1]-1))
        return(out_align)

#NT
nt="CTGTACAAGCCGAACAGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATA"
aa=Seq(nt).translate()
aa 
class align_local3():
    aligned_data=dict()
    def __init__(self, amplicon, ref,gop=-3, gep=-1):
        self.amplicon=amplicon
        self.ref=ref
        self.gop=gop
        self.gep=gep
    def align(self):
        alignments = pairwise2.align.localms(str(self.ref), str(self.amplicon), 2, -1, self.gop, self.gep)
        try:
            alignm=format_alignment(*alignments[0])
            seq_align = alignm.split("\n")[2]
            return(seq_align)

        except IndexError:
            return(None)

df_aa=full_df_trim
gop, gep=-3,-1
strand_dir= {"3p": translate_3p,
        "5p": translate_5p}  

trans_init = strand_dir.get(str(direc), None)  # Get the chosen class, or None if input is bad
df_aa = trans_init(result, corr_frame).translate_nt()
df_aa.iloc[0,1]
df_aa_c=df_aa.copy()
for i, seq in enumerate(df_aa_c.iloc[:,1]):
    if seq=='':
        df_aa=df_aa_c.drop(df_aa_c.index[i]) 

align_method="align_local3"
align_class = {"align_local": align_local,
        "align_local2": align_local2,
        "align_local3":align_local3,
        "align_global":align_global,
        "align_global2":align_global2}  
id_f=1
align_method="align_local3"
aligner_init = align_class.get(str(align_method), None)  # Get the chosen class, or None if input is bad
aa_data_align=pd.DataFrame(columns = list(df_aa.columns))
#aa_data_align=aa_data_align.append(df_aa.iloc[:,0])
aa_data_align['Seq_stats']=df_aa.iloc[:,0]
seq_info_dic={df_aa.columns[0]: df_aa.iloc[:,0]}
dic_aa_align=dict()
df_aa.columns
ref_fr='Frame_corr:0|LYKPNVRSRSTDDQMELDHMTTGGLHAYP'
for i, aa_seq in enumerate(df_aa.iloc[:,1]):
    print("i: " + str(i) + "_" + aa_seq)
    #yield iteratively the header of the certain seq and the corresponding seq
    #header=">"+ str(id_f)+"_CluSeq:" + str((round(full_df.iloc[seq_i,-4],5))) + "_var:"+str((round(full_df.iloc[seq_i,-2],5))) +"_sd:" + str((round(full_df.iloc[seq_i,-1],5)))
    #seq_obj_align = aligner_init(full_df.iloc[seq_i,0], target_sequence, gop, gep).align()
    #seq_obj_align = aligner_init(aa_seq, ts, gop, gep).align()
    seq_obj_align = re.sub(r'[(\d|\s]', '', seq_obj_align) #remove digits from the string caused by the alignment and empty spaces from the start
    matches=SequenceMatcher(None,ref_fr.split("|")[1],aa_seq)
    matches.get_matching_blocks()
    range_line=0
    seqs=[]
    alignments_per_ref=[]
    for i in range(len(matches.get_matching_blocks())):
        match=matches.get_matching_blocks()[i]
        seqs.append(len(ref_fr.split("|")[1][range_line:match.a])*"-"+seq_obj_align[match.b:match.b+match.size])
        range_line=match.a+match.size
    alignments_per_ref.append(''.join(seqs))
    alignments_per_ref= str(alignments_per_ref).replace("['", "").replace("']", "")
    #seqs_in_frame.append(str(alignments_per_ref))

    if i==1:
        break

dic_aa_align[ref_fr] = seqs_in_frame
# i+=1
# if i==1:
    #     break
dic_aa_align['Frame_corr:0|LYKPNVRSRSTDDQMELDHMTTGGLHAYP']
dic_aa_align.update(seq_info_dic)
dic_aa_align.keys()
df_aa_align=pd.DataFrame(dic_aa_align)
first_column = df_aa_align.pop('Seq_stats')
# insert column using insert(position,column_name,
# first_column) function
df_aa_align.insert(0, 'Seq_stats', first_column)
df_aa_align.to_csv(out_csv)

for i, ref_fr in enumerate(df_aa_align.columns[1:], start=1):
    aa_file="aligned/AA/fasta/" +result.split("/")[-1].split(".")[-2] + '_' +ref_fr.split("|")[0] + ".fasta"
    with open(aa_file, "w") as f:
        f.write(">0_Ref_" + ref_fr.split("|")[0] + "\n")
        f.write(ref_fr.split("|")[1] + "\n")
        for a, aa_seq in enumerate(list(df_aa_align.iloc[:,i])):
            f.write(">"+ df_aa_align.iloc[a,0] + "\n")
            f.write(aa_seq + "\n")
    mview_file= "aligned/AA/html/" +aa_file.split("/")[-1].split(".")[-2] + ".html"
    mview_command='/media/data/AtteR/Attes_bin/mview -in fasta -html head -css on -reference 1 -coloring identity ' + aa_file + '>' + mview_file
    call([mview_command], shell=True)
    print("Alignments created in html format! Files found inside aligned/AA directory")


translate_NT(result, corr_frame,direc, out_csv)
####################



####################
#SP1 5' 
assay_end = '5p'
filterlitteral = 'GCCTAGGCTAAGAACTCCTCCGCGCCACCATGGTGAGCAAGGGCGAGGAGGATAACATGG'
#rev compl of prev rliteral
lliteral=' literal=CCATGTTATCCTCCTCGCCC'
rliteral = ' literal=GTGTCTCCGGTCCCCAAAAT'
read_fwd = True
direc="5p"
#rev compl
target_sequence="tgctcaccatggtggcgcggaggagttcttagcctAGGCTAAGAACTCCTCTGAGGCAGAAGCCGGAAGGGAGCAGAGCCGGCGGCTGCAGCG"
target_sequence=target_sequence.upper()
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP1_5p"

############
#read preprocessing for each sample: trim, record read counts before and after trimming, cluster the reads 
#using starcode, calculate percentage, merge into the full dataframe containing read count-and percentage for
#each sample.
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral, rliteral, read_fwd, direc)
df_trim_full2=calculate_perc_sd(df_full)
result="unaligned/mCherry_5p_2+.fasta"
save_fasta(result, df_trim_full2, target_sequence)

csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
df_trim_full2.to_csv(csv_file)

#NT
####################
output_path="aligned/NT/"
result=output_path+"mCherry_5p_2+_prim.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local2", result, output_path,lliteral, rliteral, 4,2)
####################



#AA
#446%3
#need to make a reverse complement of an imported read?
####################
corr_frame=2
result="unaligned/mCherry_5p_2+.fasta"
output_html="aligned/AA/mCherry_5p_2+_AA.html"
out_csv="aligned/AA/mCherry_5p_2+_AA.csv"
translate_NT(result, corr_frame,direc, out_csv)
####################



