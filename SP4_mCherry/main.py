from cmath import nan
from heapq import merge
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
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP4_3p"
filterlitteral = 'CTCCCTCCACACGTGCATCTCACGCTTGACCCAGCGCTCCAGGTTGGCGATGGT' #region prior to r2 primer
lliteral=" literal=CGGCGGCATGGACGAG"
rliteral=" literal=CTGGGTCAAGCGTGAGATG"

target_sequence="ctgtacaagATGGAGCTGGACCATATGACCCGGTGCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCG"
target_sequence=target_sequence.upper()
#########
#########

#The counts and percs for each animal regarding cluster seq are identical

#########
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral, rliteral, read_fwd)
#starcode produces files with exactly same seqs, counts and percs
#files do have different number of reads but for some reason starcode clusters all of them as the same
df_trim_full2=calculate_perc_sd2(df_full)
result="unaligned/Exp2_3p_mcherry_SP4.fasta"
save_fasta(result, df_trim_full2, target_sequence)

csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
df_trim_full2.to_csv(csv_file)

#NT
####################
output_path="aligned/NT/"
# result=output_path+"Exp2_3p_mcherry_mcherry_SP4_local3.fasta"
# aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, -3,-1)
result=output_path + "Exp2_3p_mcherry_mcherry_SP4_local2_prim.fasta"
aligner(df_trim_full2, target_sequence, "align_local2", result, output_path,lliteral, rliteral, 3,1)
#AA
####################
corr_frame=0
result="unaligned/Exp2_3p_mcherry_SP4.fasta"
output_html="aligned/AA/Exp2_3p_mcherry_SP4_AA.html"
out_csv="aligned/AA/Exp2_3p_mcherry_SP4_AA.csv"
df_aa=translate_nt_aa_csv(result,corr_frame, out_csv)
df_aa.columns
df_aa["Frame_corr:0|LYKMELDHMTRCTGGLHAYPAPRGGPAAKPNVILQIGKCRAEMLEHVRRTHRHLLTEVSKQVERELKGLHRSVGKLENNLDGYVPTGDSQRWKKSIKACLCRCQETIANLE"][1]
alt_frames=[0,1,2]
alt_frames.remove(corr_frame)
#first create a dict which has keys containing the translated ref and frame info and value containing a list of amplicons translated in
#same ref
ref_aa_cor=dict()
aa_ampls=[]
seq_info=[] #this will be merged later with the final dict that has the aligned AAs against the ref at certain frame
for record in SeqIO.parse(result, "fasta"):
    if "Ref" in record.description:
        #refs_aa_frames["Frame:" + str(alt_frame)]=str(Seq(record.seq[alt_frame:]).translate())
        ref_key="Frame_corr:" + str(corr_frame) +"|" +str(Seq(record.seq[corr_frame:]).translate())
    else:
        seq_info.append(record.description)
        aa_ampls.append(str(Seq(record.seq[corr_frame:]).translate()))
ref_aa_cor[ref_key]=aa_ampls
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
ref_aa_cor.update(ref_aa)
ref_x_alignment={}

ref_aa_cor.keys()
frame_ref="Frame_corr:0|LYKMELDHMTRCTGGLHAYPAPRGGPAAKPNVILQIGKCRAEMLEHVRRTHRHLLTEVSKQVERELKGLHRSVGKLENNLDGYVPTGDSQRWKKSIKACLCRCQETIANLE"
r="LYKMELDHMTRCTGGLHAYPAPRGGPAAKPNVILQIGKCRAEMLEHVRRTHRHLLTEVSKQVERELKGLHRSVGKLENNLDGYVPTGDSQRWKKSIKACLCRCQETIANLE"
s="PRHGRAVQDGAGPYDPAGGLHAYPA"
s="LRHGRAVQDGAGPYDHVLELEIAS*"
matches=SequenceMatcher(None, r,s)
matches
ref_x_alignment=dict()
for frame_ref in ref_aa_cor.keys():
    alignments_per_ref=[]
    for ampl in ref_aa_cor[frame_ref]:
        matches=SequenceMatcher(None, frame_ref.split("|")[1],ampl)
        seqs=[]
        #you use range_line so that when you fill the remnants from left side of the match, you wont keep adding from
        #beginning since in the end, we merge the seq list elements into the whole alignment of the amplicon against the ref
        range_line=0

        seqs=[]
        range_line=0
        for i in range(len(matches.get_matching_blocks())):
            match=matches.get_matching_blocks()[i]
            aa_seq=str(ref_aa_cor[frame_ref])[match.b:match.b+match.size]
            odds=["['", "[", "]","', '", ",", "',", "'", "',"]
            if aa_seq in odds:
                continue
            #if any(ext in aa_seq for ext in odds):
            else:
                seqs.append(len(frame_ref.split("|")[1][range_line:match.a])*"-"+aa_seq)
                range_line=match.a+match.size
            #if there are empty elements, remove them as they mess up the next stage when joining the elements into a string
        mer_seq=''.join(seqs)
        alignments_per_ref.append(mer_seq)
    ref_x_alignment[frame_ref]=alignments_per_ref

#go over each ref's aa seqs, match them against ref and create a dict that has key with ref and frame info and a list of aligned 
#AAs against the ref as values
for frame_ref in ref_aa_cor.keys():
    alignments_per_ref=[]
    for ampl in ref_aa_cor[frame_ref]:
        matches=SequenceMatcher(None, frame_ref.split("|")[1],ampl)
        seqs=[]
        #you use range_line so that when you fill the remnants from left side of the match, you wont keep adding from
        #beginning since in the end, we merge the seq list elements into the whole alignment of the amplicon against the ref
        range_line=0
        for i in range(len(matches.get_matching_blocks())):
            match=matches.get_matching_blocks()[i]
            seqs.append(len(frame_ref.split("|")[1][range_line:match.a])*"-"+str(ref_aa_cor[frame_ref])[match.b:match.b+match.size])
            range_line=match.a+match.size
        alignments_per_ref.append(''.join(seqs))
    ref_x_alignment[frame_ref]=alignments_per_ref
ref_x_alignment.keys()
ref_x_alignment["Frame_corr:0|LYKMELDHMTRCTGGLHAYPAPRGGPAAKPNVILQIGKCRAEMLEHVRRTHRHLLTEVSKQVERELKGLHRSVGKLENNLDGYVPTGDSQRWKKSIKACLCRCQETIANLE"]
df=pd.DataFrame.from_dict(ref_x_alignment)

df["Seq_info"]=seq_info
# shift column 'Name' to first position
first_column = df.pop('Seq_info')

# insert column using insert(position,column_name,
# first_column) function
df.insert(0, 'Seq_info', first_column)
df.to_csv(out_csv)

print("data frame generated and saved as csv at " + out_csv)



f = open(output_html,"w")

color_scheme={"RED": [["#FF0000"], ["D", "E"]], "BLUE":[["#6495ED"],["K", "R"]], "GREEN":[["#9FE2BF"], ["C", "V", "I", "L", "P", "F", "Y", "M", "W"]], "ORANGE":[["#FF7F50"],["G", "A", "S", "T"]], "MAGENTA":[["#DE3163"], ["N", "Q", "H"]], "BLACK": [["#000000"]]}

for info, a1,a2 in zip(df["Seq_info"], list(df.iloc[:,1]), [df.columns[1].split(":")[2]]*len(list(df.iloc[:,1]))):
    aa1_list=[]
    aa2_list=[]

    a1,a2=a1[0:round(len(a1)*0.6)], a2[0:round(len(a2)*0.6)]
    f.write("<html> \n <body> <p>"+ "_".join(info.split("_")[:1]) +"</p> <p>")
    #test_html.append('<!DOCTYPE html><html> <head><link rel="stylesheet" href="format.css"></head> <meta charset="utf-8"> <body> <p>'+ "_".join(info.split("_")[:2]) + "<br>" + frame +'</p> <p>')
#go over the individual AAs of the alignemnts
    #write the html colouring into a string var for each aa and append into a list, then join and make a massive string,
    # then after all AAs have been iterated over for the certain seq, then write into the html 
    for aa1,aa2, in zip(a1,a2):
        print("aa1: "+ aa1)
        print("aa2: "+ aa2)
        # try:
        #     aa2_list.append('<span style="color:'+ color_scheme[find_colour(aa2,color_scheme)][0][0]+ '">' + aa2 + '</span>')
        # except KeyError:
        #     print("Keyerror with aa2 being: " + aa2)
        # try:
        #     aa1_list.append('<span style="color:'+ color_scheme[find_colour(aa1,color_scheme)][0][0]+ '">' + aa1 + '</span>')
        # except KeyError:
        #     print("Keyerror with aa1 being: " + aa1)
        
        aa2_list.append('<span style="color:'+ color_scheme[find_colour(aa2,color_scheme)][0][0]+ '">' + aa2 + '</span>')
        aa1_list.append('<span style="color:'+ color_scheme[find_colour(aa1,color_scheme)][0][0]+ '">' + aa1 + '</span>')

        #print("============")
    coloured_ref="".join(aa2_list)

    coloured_seq="".join(aa1_list)

    f.write(coloured_ref +"<br>")
    f.write(coloured_seq)
    f.write("</p>")

f.write("</body></html>")
f.close()


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
lliteral = ' literal=GTGTCTCCGGTCCCCAAAAT'
rliteral = ' literal=GGGCGAGGAGGATAACATGG'
filterlitteral = 'CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAG'
filterlitteral='AGCCTGTTAACAGGCGCGCCACCATGGTGAGCAAGGGCGAGGAGGATAACATGG'
target_sequence = "GTGTCTCCGGTCCCCAAAATCCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAGGGGCGGGTCCTTGAGGGCACCCACGGGAGGGGAGCGAGTAGGCGCGGAAGGCGGGGCCTGCGGCAGGAGAGGGCGCGGGCGGGCTCTGGCGCGGAGCCTGGGCGCCGCCAATGGGAGCCAGGGCTCCACGAGCTGCCGCCCACGGGCCCCGCGCAGCATAAATAGCCGCTGGTGGCGGTTTCGGTGCAGAGCTCAAGCGAGTTCTCCCGCAGCCGCAGTCTCTGGGCCTCTCTAGCTTCAGCGGCGACGAGCCTGCCACACTCGCTAAGCTCCTCCGGCACCGCACACCTGCCACTGCCGCTGCAGCCGCCGGCTCTGCTCCCTTCCGGCTTCTGCCTCAGAGGAGTTCTTAGCCTGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACGTCATATGGTCCAGCTCGGgtgagcaagggcgaggaggataacatgg"
target_sequence="CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAGGGGCGGGTCCTTGAGGGCACCCACGGGAGGGGAGCGAGTAGGCGCGGAAGGCGGGGCCTGCGGCAGGAGAGGGCGCGGGCGGGCTCTGGCGCGGAGCCTGGGCGCCGCCAATGGGAGCCAGGGCTCCACGAGCTGCCGCCCACGGGCCCCGCGCAGCATAAATAGCCGCTGGTGGCGGTTTCGGTGCAGAGCTCAAGCGAGTTCTCCCGCAGCCGCAGTCTCTGGGCCTCTCTAGCTTCAGCGGCGACGAGCCTGCCACACTCGCTAAGCTCCTCCGGCACCGCACACCTGCCACTGCCGCTGCAGCCGCCGGCTCTGCTCCCTTCCGGCTTCTGCCTCAGAGGAGTTCTTAGCCTGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACGTCATATGGTCCAGCTCGGgtgagcaa"
target_sequence=target_sequence.upper()
read_fwd = True
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP4_5p"
#########
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral, rliteral, read_fwd)
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
df_aa=translate_nt_aa_csv(result,corr_frame, out_csv)
####################
