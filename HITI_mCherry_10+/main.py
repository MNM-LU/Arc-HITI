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
scripts_dir="/media/data/AtteR/projects/hiti/pipeline_output_reorg/hiti-arc-analysis"
os.chdir(scripts_dir)
from scripts_hiti import *
#sample directory (inside which all the subdirectories exist)
sample_dir=scripts_dir + "/HITI_mCherry_10+"
os.chdir(sample_dir)
#path to the analysis folder
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP4_3p"
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
df_trim_full=calculate_perc_sd(df_full)
result="unaligned/HITI_mCherry_3p_10+.fasta"
save_fasta(result, df_trim_full, target_sequence)

csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
df_trim_full.to_csv(csv_file)
full_df_trim=pd.read_csv(csv_file, index_col=[0])
#NT
####################
output_path="aligned/NT/"
result=output_path + "HITI_mCherry_3p_10+_prim.fasta"
aligner(df_trim_full, target_sequence, "align_local2", result, output_path,lliteral, rliteral, 3,1)


#AA
####################
corr_frame=0
result="unaligned/HITI_mCherry_3p_10+.fasta"
output_html="aligned/AA/HITI_mCherry_3p_10+_AA.html"
out_csv="aligned/AA/HITI_mCherry_3p_10+_AA.csv"
aa_fasta="unaligned/HITI_mCherry_3p_10+_AA.fasta"
#translate and save results as csv
#df_aa=translate_nt_aa_csv(result,corr_frame, out_csv)

translate_NT(result, corr_frame,direc, out_csv)

#need to attach the translated ref on top, then add the translated amplicons. give the function
#to visualise based on scar vs not. if scar is given, get the longest length of a match from the 
#amplicon in frame, get the length coordinates and start counting from the second alignment (longest match)
#after this and attach this next to the initial amplicon. do the same with the last frame (starting point being 
# the end of the correct matching amplicon's length)

def translation_new(df, output_html):
    color_scheme={"RED": [["#FF0000"], ["D", "E"]], "BLUE":[["#6495ED"],["K", "R"]], "GREEN":[["#9FE2BF"], ["C", "V", "I", "L", "P", "F", "Y", "M", "W"]], "ORANGE":[["#FF7F50"],["G", "A", "S", "T"]], "MAGENTA":[["#DE3163"], ["N", "Q", "H"]], "BLACK": [["#000000"]]}
    for i, ref_fr in enumerate(df.columns[1:], start=1):
        out_html_ed="/".join(output_html.split("/")[:-1]) + "/" + output_html.split("/")[-1].split(".")[0] + str(ref_fr.split("|")[0]) + ".html"
        f = open(out_html_ed,"w")
        f.write("<html> \n <body> <p>"+ "_".join(info.split("_")[:1]) +"</p>")

        print(out_html_ed)
        for info, a1,a2 in zip(df["Seq_info"], list(df.iloc[:,i]), [ref_fr.split("|")[1]]*len(list(df.iloc[:,1]))):
            aa1_list=[]
            aa2_list=[]
            #a1,a2=a1[0:round(len(a1)*0.6)], a2[0:round(len(a2)*0.6)]

            #get the hybrid alignment prior
            f.write("<p>")

            #test_html.append('<!DOCTYPE html><html> <head><link rel="stylesheet" href="format.css"></head> <meta charset="utf-8"> <body> <p>'+ "_".join(info.split("_")[:2]) + "<br>" + frame +'</p> <p>')
        #go over the individual AAs of the alignemnts
            #write the html colouring into a string var for each aa and append into a list, then join and make a massive string,
            # then after all AAs have been iterated over for the certain seq, then write into the html 
            for aa1,aa2, in zip(a1,a2):
                # print("aa1: "+ aa1)
                # print("aa2: "+ aa2)
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

def aa_coloring(seq_infos, aligs1, aligs2, frame_info, output):
    #using lesk color code scheme:
    #color_scheme={"RED": [["#FF0000"], ["D", "E"]], "BLUE":[["#6495ED"],["K", "R"]], "GREEN":[["#9FE2BF"], ["C", "V", "I", "L", "P", "F", "Y", "M", "W"]], "ORANGE":[["#FF7F50"],["G", "A", "S", "T"]], "MAGENTA":[["#DE3163"], ["N", "Q", "H"]], "BLACK": [["#000000"],["-","|", "*"]]}
    color_scheme={"RED": [["#FF0000"], ["D", "E"]], "BLUE":[["#6495ED"],["K", "R"]], "GREEN":[["#9FE2BF"], ["C", "V", "I", "L", "P", "F", "Y", "M", "W"]], "ORANGE":[["#FF7F50"],["G", "A", "S", "T"]], "MAGENTA":[["#DE3163"], ["N", "Q", "H"]], "BLACK": [["#000000"]]}
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
out="/media/data/AtteR/projects/hiti/AA_3p_mcherry_coloured_alignment.html"


#now you have a df where first col contains the information about the specific amplicon, the cols 1-3 are
#named based on the ref codon frame + the entire translated AA seq

#Go over each alignment col by row, find the longest match and align these next to each other.
#then move onto the next col (ref seq translated using a different frame), find the longest consequtive seq
#in the amplicon and align this to the ref too. repeat with the last col as well.

#when aligning to the ref using e.g. ref1 (in frame) and amplicon 1, cut off the amplicon 1 after
#the longest conseq. seq has been found. Then take the same amplicon but when it has been aligned to the
#ref2 (out of frame), find the longest seq, remove the extra "---" from the beginning relative to the 
#amplicon 1 and then align to the ref. merge the two. 

def visualise_aa_hybrid_alignments(df, output_html): #output as html
    seqinfos=[]
    aligs_merged_seq=[]
    aligs_merged_ref=[]
    frame_info=[]
    for ampl_row in range(round(len(df.index)*0.1)):
        match1 = SequenceMatcher(None, df.columns[1].split(":")[-1], df.iloc[ampl_row,1]).find_longest_match(0, len(df.columns[1].split(":")[-1]), 0, len(df.iloc[ampl_row,1]))
        end_of_match=match1.b+match1.size
        matched_ampl1= len(df.columns[1].split(":")[-1][:match1.a])*"-" + str(df.iloc[ampl_row,1][match1.b:]) + "|"
        seq_inframe1=matched_ampl1[0:end_of_match] + "|"
        ref_inframe1=df.columns[1].split(":")[-1][0:end_of_match] + "|"
        #get the matches with the other reading frames (i.e. the last 2 cols of the df)
        match2_a = SequenceMatcher(None, df.columns[2].split(":")[-1], df.iloc[ampl_row,2]).find_longest_match(0, len(df.columns[2].split(":")[-1]), 0, len(df.iloc[ampl_row,2]))
        match2_b = SequenceMatcher(None, df.columns[3].split(":")[-1], df.iloc[ampl_row,3]).find_longest_match(0, len(df.columns[3].split(":")[-1]), 0, len(df.iloc[ampl_row,3]))
        
        #modify this part so that you include AAs generated by all different frames
        if match2_a.size>match2_b.size:
            matched_ampl2= len(df.columns[2].split(":")[-1][:match2_a.a])*"-" + str(df.iloc[ampl_row,2][match2_a.b:]) + "|"
            ref_outframe2=df.columns[2].split(":")[-1][0:end_of_match] + "|"
            # seq_inframe2=matched_ampl2[0:end_of_match]
            seq_outframe=matched_ampl2[end_of_match:]
            frame_info.append("Frame:+1 " + ref_inframe1.index("|")*" " + "Frame:0")
        else:
            matched_ampl2= len(df.columns[3].split(":")[-1][:match2_b.a])*"-" + str(df.iloc[ampl_row,3][match2_b.b:]) + "|"
            ref_outframe2=df.columns[3].split(":")[-1][end_of_match:] + "|"
            seq_outframe=matched_ampl2[end_of_match:]
            # seq_inframe2=matched_ampl2[0:end_of_match]
            seq_outframe=matched_ampl2[end_of_match:]
            frame_info.append("Frame:+1 " + ref_inframe1.index("|")*" " + "Frame:+2")
            print("Frame:+1 " + ref_inframe1.index("|")*" " + "Frame:+2")
        merged_align_seq=seq_inframe1+seq_outframe

        merged_align_ref=ref_inframe1+ref_outframe2
        header_info=df.iloc[ampl_row,0]
        #add one more line which adds the frames

        seqinfos.append(df.iloc[ampl_row,0])
        #aligs_merged_seq.append(merged_align_seq[0:round(len(aligs_merged_seq)*0.6)])
        aligs_merged_seq.append(merged_align_seq)
        #aligs_merged_ref.append(merged_align_ref[0:round(len(merged_align_ref)*0.6)])
        aligs_merged_ref.append(merged_align_ref)
    aa_coloring(seqinfos,aligs_merged_seq,aligs_merged_ref, frame_info, output_html)

##############################












#5'
#############
#User configurable variables

#path to the analysis folder
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP4_5p"
#############


#SP4 5'
direc = "5p"
lliteral = ' literal=CCATGTTATCCTCCTCGCCC'
rliteral = ' literal=ATTTTGGGGACCGGAGACAC'
filterlitteral="CCATGTTATCCTCCTCGCCCTTGCTCACCCGAGCTGGACCATATGACGTCATATGGT"
target_sequence="tgctcacCCGAGCTGGACCATATGACGTCATATGGTCCAGCTCCATCTGGTCGTCGGTGCTGCGGCTCCGAACAGGCTAAGAACTCCTCTGAGGCAGAAGCCGGAAGGGAGCAGAGCCGGCGGCTGCAGCGGCAGTGGCAGGTGT"
target_sequence=target_sequence.upper()
read_fwd = True
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP4_5p"

############
#read preprocessing for each sample: trim, record read counts before and after trimming, cluster the reads 
#using starcode, calculate percentage, merge into the full dataframe containing read count-and percentage for
#each sample.

df_full=import_reads_process_mini(base_path, target_sequence,filterlitteral,lliteral,rliteral,read_fwd, direc)
df_trim_full=calculate_perc_sd(df_full)
result="unaligned/HITI_mCherry_5p_10+.fasta"
save_fasta(result, df_trim_full, target_sequence)

csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
df_trim_full.to_csv(csv_file)

#NT
####################
output_path="aligned/NT/"
result=output_path + "HITI_mCherry_5p_10+_prim.fasta"
aligner(df_trim_full, target_sequence, "align_local2", result, output_path, lliteral, rliteral,4,2)
####################

#AA
#484%3
####################
corr_frame=1
result="unaligned/HITI_mCherry_5p_10+.fasta"
output_html="aligned/AA/HITI_mCherry_5p_10+_AA.html"
out_csv="aligned/AA/HITI_mCherry_5p_10+_AA.csv"
translate_NT(result, corr_frame,direc, out_csv)

####################
