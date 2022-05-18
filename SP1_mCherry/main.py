from cmath import nan
from heapq import merge
import pandas as pd
import os 
os.getcwd()
os.chdir("/media/data/AtteR/projects/hiti/pipeline_output_reorg/hiti-arc-analysis")
from scripts_hiti import *
from scripts_hiti import import_reads_process_mini
from scripts_hiti import calculate_perc_sd2

os.chdir("/media/data/AtteR/projects/hiti/pipeline_output_reorg/hiti-arc-analysis/SP1_mCherry")


####################
#SP1 3' 
'''
SP1
3p
ctgtacaagccgaacgttcggagccgcagcaccgacgaccagATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAG
lliteral cggcggcatggacgag rliteral GTCTTCTACCGTCTGGAGAGG

'''
trimmed_df="/media/data/AtteR/projects/hiti/dfs/full_df_trim_mcherry_p3_HITI_SP1.csv"
transgene = 'mCherry'
assay_end = '3p'
read_fwd = True
lliteral=" literal=CGGCGGCATGGACGAGC"
rliteral=" literal=GTCTTCTACCGTCTGGAGAGG"

base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP1_3p"
filterlitteral = 'CTCCCTCCACACGTGCATCTCACGCTTGACCCAGCGCTCCAGGTTGGCGATGGT' #region prior to r2 primer

target_sequence="ctgtacaagccgaacgttcggagccgcagcaccgacgaccagATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAG"
target_sequence=target_sequence.upper()

#########
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral,rliteral,read_fwd)
df_trim_full2=calculate_perc_sd2(df_full)
result="unaligned/Exp2_3p_mcherry_mcherry_SP1.fasta"
save_fasta(result, df_trim_full2, target_sequence)
####################
csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
df_trim_full2.to_csv(csv_file)

#NT
output_path="aligned/"
#result="aligned/Exp2_3p_mcherry_SP1_local3.fasta"
#test_res=aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, -3,-1)
result=output_path + "Exp2_3p_mcherry_SP1_local2_prim.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local2", result, output_path, lliteral, rliteral, 3,1)
####################
#AA
####################
corr_frame=1
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/unaligned_seqs/Exp2_3p_mcherry_mcherry_SP1.fasta"
output_html="/media/data/AtteR/projects/hiti/pipeline_output_reorg/AA_aligned_html/Exp2_3p_mcherry_mcherry_SP1_AA.html"
translate_nt_aa_hiti2(result, corr_frame, output_html)
####################


####################
#SP1 5' 
'''
5p
gcagagctcaagcgagttctcccgcagccgcagtctctgggcctctctagcttcagcggcgacgagcctgccacactcgctaagctcctccggcaccgcacacctgccactgccgctgcagccgccggctctgctcccttccggcttctgcctcagaggagttcttagcctaggctaagaactcctccgcgccaccatggtgagcaa
rliteral gggcgaggaggataacatgg
'''
assay_end = '5p'
filterlitteral = 'CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAG'
rliteral = ' literal=GGGCGAGGAGGATAACATGG'
lliteral = ' literal=GTGTCTCCGGTCCCCAAAAT'

target_sequence= "CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAGGGGCGGGTCCTTGAGGGCACCCACGGGAGGGGAGCGAGTAGGCGCGGAAGGCGGGGCCTGCGGCAGGAGAGGGCGCGGGCGGGCTCTGGCGCGGAGCCTGGGCGCCGCCAATGGGAGCCAGGGCTCCACGAGCTGCCGCCCACGGGCCCCGCGCAGCATAAATAGCCGCTGGTGGCGGTTTCGGTGCAGAGCTCAAGCGAGTTCTCCCGCAGCCGCAGTCTCTGGGCCTCTCTAGCTTCAGCGGCGACGAGCCTGCCACACTCGCTAAGCTCCTCCGGCACCGCACACCTGCCACTGCCGCTGCAGCCGCCGGCTCTGCTCCCTTCCGGCTTCTGCCTCAGAGGAGTTCTTAGCCTaggctaagaactcctccgcgccaccatggtgagcaa"
target_sequence=target_sequence.upper()
read_fwd = True

base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP1_5p"
#########
#gives an error: no columns to parse from file when trying to import starcode file as a df and then calculate count and perc

df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral, rliteral, read_fwd)
df_trim_full2=calculate_perc_sd2(df_full)
result="unaligned/pipeline_output_reorg/fastas/Exp2_5p_mcherry_mcherry_SP1.fasta"
save_fasta(result, df_trim_full2, target_sequence)

csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
df_trim_full2.to_csv(csv_file)

#NT
####################
output_path="aligned/"
# result=output_path + "Exp2_5p_mcherry_mcherry_SP1_local3.fasta"
# test_res=aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, -3,-1)
result=output_path+"Exp2_5p_mcherry_mcherry_SP1_local2.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local2", result, output_path,lliteral, rliteral, 4,2)
####################

#atm if theres a mismatch with NTs in rows colours the first NT mismatch is coloured correctly as red but the next one is not!


#AA
####################
corr_frame=1
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp2_5p_mcherry_mcherry_SP1.fasta"
output_html="/media/data/AtteR/projects/hiti/pipeline_output_reorg/AA_aligned_html/Exp2_5p_mcherry_mcherry_SP1_AA.html"
translate_nt_aa_hiti2(result, corr_frame, output_html)
####################



#AA draft part, ignore.
from Bio import SeqIO
from difflib import Differ, SequenceMatcher
from Bio.Seq import Seq
all_frames=[0,1,2]
#out_of_frames.remove(corr_frame) dont take out any of the frames as the seqs after the scar may be in frame
refs_aa_frames={}
aa_and_perc={}

len(aa_and_perc)
for record in SeqIO.parse(result, "fasta"):
    for alt_frame in all_frames:
        if record.id=="0":
            refs_aa_frames["Frame:" + str(corr_frame)]=str(Seq(record.seq[corr_frame:]).translate())
        refs_aa_frames["Frame:" + str(alt_frame)]=str(Seq(record.seq[alt_frame:]).translate())
        aa_and_perc[">"+str(record.description) + "_transl.frame:" + str(alt_frame)]=str(Seq(record.seq[alt_frame:]).translate())

#you go over the ref seqs in different frames and align all the amplicons to them. save the alignment into the df's specific column. initially
#save into a list or such
ref_x_alignment={}
for frame_ref in refs_aa_frames.keys():
    alignments_per_ref=[]
    for ampl in aa_and_perc.keys():
        matches=SequenceMatcher(None,refs_aa_frames[frame_ref],aa_and_perc[ampl])
        seqs=[]
        #you use range_line so that when you fill the remnants from left side of the match, you wont keep adding from
        #beginning since in the end, we merge the seq list elements into the whole alignment of the amplicon against the ref
        range_line=0
        for i in range(len(matches.get_matching_blocks())):
            match=matches.get_matching_blocks()[i]
            seqs.append(len(refs_aa_frames[frame_ref][range_line:match.a])*"-"+str(aa_and_perc[ampl])[match.b:match.b+match.size])
            range_line=match.a+match.size
        alignments_per_ref.append(''.join(seqs))
    ref_x_alignment[frame_ref + "|Ref:" +refs_aa_frames[frame_ref]]=alignments_per_ref
seq_info={"Seq_info:":aa_and_perc.keys()}
keys=list(aa_and_perc.keys())

seq_info=["Seq_info"]+list(aa_and_perc.keys())
ref_x_alig_list=[]
for keys, values in ref_x_alignment.items():
    #print("".join(list((keys))[:]))
    #ref_x_alig_list.append([("".join(list((keys))))]+list(values))
    ref_x_alig_list.append([keys]+list(values))

df = pd.DataFrame(data= {seq_info[0]: seq_info[1:], ref_x_alig_list[0][0]:ref_x_alig_list[0][1:], ref_x_alig_list[1][0]:ref_x_alig_list[1][1:], ref_x_alig_list[2][0]:ref_x_alig_list[2][1:]})
df.columns

#to show alignments to ref seq that has been translated in frame

def find_colour(aa,color_scheme):
    col=[key for key,val in color_scheme.items() if any(aa in s for s in val)]
    out="".join(col)
    if out==None or not col:
        out="BLACK"
        return(out)
    else:
        return(out)

def aa_coloring(seq_infos, aligs1, aligs2, frame_info, output):
    #using lesk color code scheme:
    #color_scheme={"RED": [["#FF0000"], ["D", "E"]], "BLUE":[["#6495ED"],["K", "R"]], "GREEN":[["#9FE2BF"], ["C", "V", "I", "L", "P", "F", "Y", "M", "W"]], "ORANGE":[["#FF7F50"],["G", "A", "S", "T"]], "MAGENTA":[["#DE3163"], ["N", "Q", "H"]], "BLACK": [["#000000"],["-","|", "*"]]}
    color_scheme={"RED": [["#FF0000"], ["D", "E"]], "BLUE":[["#6495ED"],["K", "R"]], "GREEN":[["#9FE2BF"], ["C", "V", "I", "L", "P", "F", "Y", "M", "W"]], "ORANGE":[["#FF7F50"],["G", "A", "S", "T"]], "MAGENTA":[["#DE3163"], ["N", "Q", "H"]], "BLACK": [["#000000"], ["-","|", "*"]]}
    print(output)
    f = open(output,"w")
    for info, a1,a2, frame in zip(seq_infos, aligs1, aligs2,frame_info):
        aa1_list=[]
        aa2_list=[]

        print(frame)
        a1,a2=a1[0:round(len(a1)*0.6)], a2[0:round(len(a2)*0.6)]
        #f.write("<html> \n <body> <p>"+ "_".join(info.split("_")[:1]) +"</p> <p>")
        f.write('<!DOCTYPE html><html> <head><link rel="stylesheet" href="format.css"></head> <meta charset="utf-8"> <body> <p>'+ "_".join(info.split("_")[:2]) + "<br>" + frame +'</p> <p>')

    #go over the individual AAs of the alignemnts
        #write the html colouring into a string var for each aa and append into a list, then join and make a massive string,
        # then after all AAs have been iterated over for the certain seq, then write into the html 
        for aa1,aa2, in zip(a1,a2):
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
out="/media/data/AtteR/projects/hiti/pipeline_output_reorg/SP1_mCherry/aligned/coloured_alignment.html"


#now you have a df where first col contains the information about the specific amplicon, the cols 1-3 are
#named based on the ref codon frame + the entire translated AA seq

#Go over each alignment col by row, find the longest match and align these next to each other.
#then move onto the next col (ref seq translated using a different frame), find the longest consequtive seq
#in the amplicon and align this to the ref too. repeat with the last col as well.

#when aligning to the ref using e.g. ref1 (in frame) and amplicon 1, cut off the amplicon 1 after
#the longest conseq. seq has been found. Then take the same amplicon but when it has been aligned to the
#ref2 (out of frame), find the longest seq, remove the extra "---" from the beginning relative to the 
#amplicon 1 and then align to the ref. merge the two. 
seqinfos=[]
aligs_merged_seq=[]
aligs_merged_ref=[]
frame_info=[]

for ampl_row in range(round(len(df.index)*0.6)):
    aligs_merged_ref.append(df.columns[1].split(":")[-1] + "|" + df.columns[2].split(":")[-1] + "|" + df.columns[3].split(":")[-1])
    aligs_merged_seq.append(df.iloc[ampl_row,1] + "|" + df.iloc[ampl_row,2] + "|" + df.iloc[ampl_row,3])
    seqinfos.append(df.iloc[ampl_row,0])
    frame_info.append(df.columns[1].split("|")[0] + len(df.columns[1].split(":")[-1])*" " + df.columns[2].split("|")[0] + len(df.columns[2].split(":")[-1])*" " + df.columns[3].split("|")[0])

aligs_merged_seq[0]
aligs_merged_ref[0]
CGMAWTSCTSRTFGAAAPTTRWSW|LRHGMDELYKPNVRSRSTDD
-G-A-T----RT--AAAPT-R---|LRHG---L--P--RSR-T--
#using lesk color code scheme:
#color_scheme={"RED": [["#FF0000"], ["D", "E"]], "BLUE":[["#6495ED"],["K", "R"]], "GREEN":[["#9FE2BF"], ["C", "V", "I", "L", "P", "F", "Y", "M", "W"]], "ORANGE":[["#FF7F50"],["G", "A", "S", "T"]], "MAGENTA":[["#DE3163"], ["N", "Q", "H"]], "BLACK": [["#000000"],["-","|", "*"]]}
color_scheme={"RED": [["#FF0000"], ["D", "E"]], "BLUE":[["#6495ED"],["K", "R"]], "GREEN":[["#9FE2BF"], ["C", "V", "I", "L", "P", "F", "Y", "M", "W"]], "ORANGE":[["#FF7F50"],["G", "A", "S", "T"]], "MAGENTA":[["#DE3163"], ["N", "Q", "H"]], "BLACK": [["#000000"]]}
f = open(output,"w")
for info, a1,a2, frame in zip(seq_infos, aligs1, aligs2,frame_info):
    aa1_list=[]
    aa2_list=[]

    print(frame)
    a1,a2=a1[0:round(len(a1)*0.6)], a2[0:round(len(a2)*0.6)]
    #f.write("<html> \n <body> <p>"+ "_".join(info.split("_")[:1]) +"</p> <p>")
    f.write('<!DOCTYPE html><html> <head><link rel="stylesheet" href="format.css"></head> <meta charset="utf-8"> <body> <p>'+ "_".join(info.split("_")[:2]) + "<br>" + frame +'</p> <p>')

#go over the individual AAs of the alignemnts
    #write the html colouring into a string var for each aa and append into a list, then join and make a massive string,
    # then after all AAs have been iterated over for the certain seq, then write into the html 
    for aa1,aa2, in zip(a1,a2):
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

out="/media/data/AtteR/projects/hiti/pipeline_output_reorg/SP1_mCherry/aligned/coloured_alignment.html"
aa_coloring(seqinfos,aligs_merged_seq, aligs_merged_ref, frame_info, out)
