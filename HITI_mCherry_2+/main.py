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

#where the programs bbduk and starcode are found
program_path="/media/data/AtteR/Attes_bin/"

#############

#SP1 3' 
read_fwd = True
lliteral=" literal=CGGCGGCATGGACGAG"
rliteral=" literal=GTCTTCTACCGTCTGGAGAGG"
direc="3p"
filterlitteral="CGGCGGCATGGACGAGCTGTACAAGCCGAACGTTCGGAGCCGCAGCACCGACGACCAG"
target_sequence="ctgtacaagccgaacGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGA"
target_sequence=target_sequence.upper()

#########
#read preprocessing for each sample: trim, record read counts before and after trimming, cluster the reads 
#using starcode, calculate percentage, merge into the full dataframe containing read count-and percentage for
#each sample. 

df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral,rliteral,read_fwd, direc, program_path)
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
test_res=aligner(df_trim_full, target_sequence, "align_local2", result, output_path, lliteral, rliteral, 3,1)
####################
#AA
####################
corr_frame=0
result="unaligned/mCherry_3p_2+.fasta"
output_html="aligned/AA/mCherry_3p_2+.html"
out_csv="aligned/AA/mCherry_3p_2+.csv"

table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', '-':'-'
}


#add the extra "-" into the NT seq so that the alignment improves
for i, seq in enumerate(full_df_trim.iloc[:,0]):
    insert_site=seq.find('AGCCGAAC') + len('AGCCGAAC')
    seq_obj= seq[:insert_site] + "-" + seq[insert_site:]
    full_df_trim[i,0]=seq_obj
full_df_trim.iloc[2,0]
full_df_trim.columns
full_df_trim.iloc[0,0]


seq="CTGTACAAGCCGAACAGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATA"
insert_site=seq.find('AGCCGAAC') + len('AGCCGAAC')
seq=seq[:insert_site] + "-" + seq[insert_site:]


protein=""

seq=seq[0:]

def tr(seq):
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if "-" in codon:
            #change reading frame
            codon = seq[i+codon.count('-'):i+codon.count('-') + 3]
            protein+=table[codon]
        else:
            protein+=table[codon]
    return(protein)

def translate_own(seq):
    if len(seq)%3 != 0:
        diff = 3 - (len(seq) - (len(seq) / 3) * 3)
        seq+=seq+int(diff)*"-"
        protein=tr(seq)
    else:
        protein=tr(seq)
    return(protein)

class translate_3p:
    def __init__(self, result, corr_frame):
        self.result=result
        self.corr_frame=corr_frame
    def translate_nt(self):
        ref_aa_cor=dict()
        aa_ampls=[]
        seq_info=[] #this will be merged later with the final dict that has the aligned AAs against the ref at certain frame
        for record in SeqIO.parse(self.result, "fasta"):
            if "Ref" in record.description:
                #refs_aa_frames["Frame:" + str(alt_frame)]=str(Seq(record.seq[alt_frame:]).translate())
                ref_key="Frame_corr:" + str(corr_frame) +"|" +translate_own(str(record.seq[corr_frame:]))
            else:
                seq_info.append(record.description)
                aa_ampls.append(translate_own(str(record.seq[corr_frame:])))

        ref_aa_cor[ref_key]=aa_ampls
        ref_aa=dict()
        seq_info_dic = {'Seq_stats': seq_info}

        seq_info_dic = {'Seq_stats': seq_info}
        ref_aa_cor.update(ref_aa)
        ref_aa_cor.update(seq_info_dic)
        aa_df=pd.DataFrame(ref_aa_cor)
        first_column = aa_df.pop('Seq_stats')
        # insert column using insert(position,column_name,
        # first_column) function
        aa_df.insert(0, 'Seq_stats', first_column)
        return(aa_df)
def translate_NT(result, corr_frame, direc, out_csv):
    strand_dir= {"3p": translate_3p,
            "5p": translate_5p}  

    trans_init = strand_dir.get(str(direc), None)  # Get the chosen class, or None if input is bad
    df_aa = trans_init(result, corr_frame).translate_nt()
    for i, seq in enumerate(df_aa.iloc[:,1]):
        if seq=='':
            df_aa=df_aa.drop(df_aa.index[i]) 

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
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral, rliteral, read_fwd, direc, program_path)
df_trim_full2=calculate_perc_sd(df_full,3)
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
# corr_frame=2
# result="unaligned/mCherry_5p_2+.fasta"
# output_html="aligned/AA/mCherry_5p_2+_AA.html"
# out_csv="aligned/AA/mCherry_5p_2+_AA.csv"
# translate_NT(result, corr_frame,direc, out_csv)
####################
