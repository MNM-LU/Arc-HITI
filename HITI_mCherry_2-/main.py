def align_trimmer(aligned_data,target_sequence):
    for id in aligned_data.keys():
        if len(aligned_data[id])==len(target_sequence):
            continue
        if len(aligned_data[id])>len(target_sequence):
            N_dashes=len(target_sequence)-len(aligned_data[id])
            aligned_data[id]=aligned_data[id][:N_dashes]
            print("Seq length larger than ref by " + str(N_dashes) + " ... \n After removal length: " + str(len(aligned_data[id][:N_dashes])))
        else:
            N_dashes=len(target_sequence)-len(aligned_data[id])
            aligned_data[id]=aligned_data[id]+ N_dashes*"-"
            print("Seq length smaller than ref by " + str(N_dashes) + " ... \n After addition length: " + str(len(aligned_data[id])))
    return(aligned_data)
def write_align(aligned_data, filename, target_sequence):
    with open(filename, "w") as handle:
        header=">0"+" Ref_seq"
        handle.write(header + "\n" + target_sequence + "\n")
        for seq_i in aligned_data.keys():
            handle.write(seq_i + "\n" + aligned_data[seq_i] + "\n")


def find_remove_duplicates(aligned_data_trim):
    id_keys_dict=dict()
    algs=list(aligned_data_trim.values())

    for i, id in enumerate(aligned_data_trim.keys()):
        id_keys_dict[i]=id

    #go over the indeces where the duplicates of a seq found, two at a time and sum them all up together
    def remove_multiples(index_pos, id_keys_dict):
        perc_sum_all=0
        var_sum_all=0
        sd_sum_all=0
        for i in index_pos:
            perc_sum_all+=float(id_keys_dict[i].split("_")[0].split(":")[1])
            var_sum_all+=float(id_keys_dict[i].split("_")[1].split(":")[1])
            sd_sum_all+=float(id_keys_dict[i].split("_")[2].split(":")[1])
        #merged_id="CluSeq:"+str(round(perc_sum_all,5))+ "_var:"+str(round(var_sum_all,5)) +"_sd:"+str(round(math.sqrt(var_sum_all),5))
        merged_id="CluSeq:"+str(round(perc_sum_all,5))+ "_sd:"+str(round(math.sqrt(var_sum_all),5))

        return(merged_id)

    #gets the indices of the matching seqs
    def get_indeces_matching_seqs(list_of_elems, element):
        index_pos_list = []
        index_pos = 0
        while True:
            try:
                # Search for item in list from indexPos to the end of list
                index_pos = list_of_elems.index(element, index_pos)
                # Add the index position in list
                index_pos_list.append(index_pos)
                index_pos += 1
            except ValueError as e:
                break
        return index_pos_list

    def remove_key_from_dict(index_pos,id_keys_dict,aligned_data_trim):
        for pos in index_pos:
            #id_keys_dict[pos] value corresponds to the key of the alignments 
            print(id_keys_dict[pos])
            print(aligned_data_trim[id_keys_dict[pos]])
            del aligned_data_trim[id_keys_dict[pos]]
        return(aligned_data_trim)
    import collections
    #merged values
    dups = [item for item, count in collections.Counter(algs).items() if count >= 2]
    merged_multiples=dict()
    for dupl_seq in dups:
        index_pos = get_indeces_matching_seqs(algs, dupl_seq)
        #take the seq from any of the matches
        #seq=algs[index_pos[0]]
        #print(dupl_seq[:40] + " found at pos " + str(index_pos))
        aligned_data_trim=remove_key_from_dict(index_pos,id_keys_dict,aligned_data_trim)
        #key includes the merged matches values and the value the seq
        merged_multiples[remove_multiples(index_pos, id_keys_dict)]=algs[index_pos[0]]
    #now we have removed the seqs that have duplicates and effectively merged them together into merged_multiples dict. This
    #is then added back to the original dict from which the duplicates had been removed.
#    merged_multiples.update(aligned_data_trim)
    aligned_data_trim.update(merged_multiples)

    return(aligned_data_trim)
#once you get the indices of the duplicates, get the keys based on these indices, merge the key values
#Now we have a script for detecting duplicates from the data dict, removing them via merging the values
#need to remove the seqs from the original dict using the index_pos approach
def add_primers_save(aligned_data_trim,filename, target_sequence,lliteral):
    id_f = 1
    ref="Ref"
    with open(filename, "w") as handle:
        whole_ref=lliteral.split("=")[1] + "-" + target_sequence
        seq_obj = SeqRecord(Seq(whole_ref), id=str(0), description=ref)
        count = SeqIO.write(seq_obj, handle, "fasta")
        for id, seq_prims in aligned_data_trim.items():
            #descr=id.split("_")[0] + "_" + id.split("_")[1]
            descr=id
            whole_seq=lliteral.split("=")[1] + "-" +seq_prims
            seq_obj = SeqRecord(Seq(whole_seq), id=str(id_f), description=descr)
            print(seq_obj)
            count = SeqIO.write(seq_obj, handle, "fasta")
            id_f+=1
    print("Saved!")

#takes into account whether the sample is 3p or 5p. if 5p, then reverse translate the primer and put on the other
#end of the translated main seq
def add_primers_save2(aligned_data_trim,filename, target_sequence,lliteral, dir, rev_complement):
    id_f = 1
    ref="Ref"
    old_chars = "ACGT"
    replace_chars = "TGCA"
    tab = str.maketrans(old_chars,replace_chars)
    if rev_complement=="True":
        target_sequence = str(target_sequence).translate(tab)[::-1]
    else:
        target_sequence=target_sequence
    with open(filename, "w") as handle:
        if '3' in dir:
            whole_ref=lliteral.split("=")[1] + "-" + target_sequence
            seq_obj = SeqRecord(Seq(whole_ref), id=str(0), description=ref)
            count = SeqIO.write(seq_obj, handle, "fasta")
            for id, seq_prims in aligned_data_trim.items():
                #descr=id.split("_")[0] + "_" + id.split("_")[1]
                descr=id
                whole_seq=lliteral.split("=")[1] + "-" +seq_prims
                seq_obj = SeqRecord(Seq(whole_seq), id=str(id_f), description=descr)
                print(seq_obj)
                count = SeqIO.write(seq_obj, handle, "fasta")
                id_f+=1
        else:
            whole_ref=target_sequence + "-" + lliteral.split("=")[1].translate(tab)[::-1] # reverse primer position as the 5p seq is sequenced in rev
            # whole_ref=target_sequence + "-" + lliteral.split("=")[1](tab)[::-1] # reverse primer position as the 5p seq is sequenced in rev
            seq_obj = SeqRecord(Seq(whole_ref), id=str(0), description=ref)
            count = SeqIO.write(seq_obj, handle, "fasta")
            for id, seq_prims in aligned_data_trim.items():
                #descr=id.split("_")[0] + "_" + id.split("_")[1]
                descr=id
                whole_seq=seq_prims+ "-" + lliteral.split("=")[1].translate(tab)[::-1]
                seq_obj = SeqRecord(Seq(whole_seq), id=str(id_f), description=descr)
                print(seq_obj)
                count = SeqIO.write(seq_obj, handle, "fasta")
                id_f+=1


    print("Saved!")



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
#sample directory (inside which all the subdirectories exist)
sample_dir=scripts_dir + "/HITI_mCherry_2-/" 
os.chdir(sample_dir)
#path to the analysis folder
base_path = '/media/data/AtteR/projects/hiti/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
export_path=sample_dir = "trimmed_data/"

#where the programs bbduk and starcode are found
# program_path="/media/data/AtteR/Attes_bin"
#############


#mCherry 3p
#############
transgene="mCherry"
read_fwd = True
direc="3p"
animal_list = [7, 8, 9, 10, 11, 12]
filterlitteral = 'CTCCCTCCACACGTGCATCTCACGCTTGACCCAGCGCTCCAGGTTGGCGATGGT'
filterlitteral='CGGCGGCATGGACGAGCTGTACAAGGtcggtgctgcggctccgCGGAGCC'
lliteral = ' literal=CGGCGGCATGGACGAG'
rliteral = ' literal=CATATGACCACCGG' 
#original, full length
target_sequence = "CTGTACAAGGtcggtgctgcggctccgCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGC"
target_sequence=target_sequence.upper()

#read preprocessing for each sample: trim, record read counts before and after trimming, cluster the reads 
#using starcode, calculate percentage, merge into the full dataframe containing read count-and percentage for
#each sample.
full_df=analyze_all(base_path, transgene, filterlitteral,lliteral,rliteral,export_path,read_fwd, animal_list, target_sequence, direc)
result="unaligned/HITI_mCherry_3p_1-.fasta"
full_df_trim=calculate_perc_sd(full_df, 3)

save_fasta(result, full_df_trim, target_sequence)
csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
full_df.to_csv(csv_file)
#full_df_trim=pd.read_csv(csv_file, index_col=[0])

#NT
####################
output_path="aligned/NT/"
result=output_path + "HITI_mCherry_3p_1-_prim.fasta"
test_res=aligner(full_df_trim, target_sequence, "align_local2", result, output_path, lliteral, rliteral, 3,1)
####################

#AA
####################
corr_frame=0
aa_primer_frame=1

result="unaligned/HITI_mCherry_3p_1-.fasta"
out_csv="aligned/AA/HITI_mCherry_3p_1-_AA.csv"
output_html="aligned/AA/HITI_mCherry_3p_1-_AA.html"

translate_NT(result, corr_frame,direc, out_csv, lliteral.split("=")[1],aa_primer_frame)


#mCherry 5p
############
transgene='mCherry'
assay_end='5p'
read_fwd=True
direc='5p'

animal_list = [1, 2, 3, 4, 5, 6] 
filterlitteral = 'CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAG'
filterlitteral = 'TTATCCTCCTCGCCCTTGCTCACCATGGTGGCGCGcctgttAACAGGCTAAG'
lliteral = ' literal=TTATCCTCCTCGCCC'
rliteral = ' literal=CCTCTGAGGCAGAA'

base_path = '/media/data/AtteR/projects/hiti/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
export_path = '/media/data/AtteR/projects/hiti/output/'
target_sequence = "TTGCTCACCATGGTGGCGCGCCTGTTAACAGGCTAAGAACTCCTCTGAGGCAGAAGCCGGAAGGGAGCAGAGCCGGCGGCTGCAGCGGCAGTGGCAGGTGTGCGGTGCCGGAGGAGCTTAGCGAGTGTGGCAGGCTCGTCGCCGCTGAAGCTAGAGAGGCCCAGAGACTGCGGCTGCGGGAGAACTCGCTTGAGCTCTGCACCGAAACCGCCACCAGCGGCTATTTATGCTGCGCGGGGCCCGTGGGCGGCAGCTCGTGGAGCCCTGGCTCCCATTGGCGGCGCCCAGGCTCCGCGCCAGAGCCCGCCCGCGCCCTCTCCTGCCGCAGGCCCCGCCTTCCGCGCCTACTCGCTCCCCTCCCGTGGGTGCCCTCAAGGACCCGCCCCTTCACAGCCCCGAGTGACTAATGTGCTCTGCTGCGCGCCTCCCACCGGGAGGGATTTTGGGGACCGGAGACAC"

animal_list = [1, 2, 3, 4, 5, 6]
target_sequence = target_sequence.upper()

############
#read preprocessing for each sample: trim, record read counts before and after trimming, cluster the reads 
#using starcode, calculate percentage, merge into the full dataframe containing read count-and percentage for
#each sample.

full_df=analyze_all(base_path, transgene, filterlitteral,lliteral,rliteral,export_path,read_fwd, animal_list, target_sequence, direc)
full_df_trim=calculate_perc_sd(full_df,3)
result="unaligned/HITI_mCherry_5p_1-.csv"
save_fasta(result, full_df_trim, target_sequence)
#########
csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
full_df_trim.to_csv(csv_file)
full_df_trim = pd.read_csv(csv_file, index_col=[0])

#NT
####################
output_path="aligned/NT/"
result=output_path+"HITI_mCherry_5p_2-_prim.fasta"

# If given the argument True, creates a rev compl
aligner(full_df_trim, target_sequence, "align_local2", result, output_path, lliteral, rliteral,3,1, "True")
####################


#AA
#439%3
# corr_frame=1
# result="unaligned/HITI_mCherry_5p_1-.fasta"
# out_csv="aligned/AA/HITI_mCherry_5p_1-_AA.csv"
# output_html="aligned/AA/HITI_mCherry_5p_1-_AA.html"
from Bio.pairwise2 import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
full_df = pd.read_csv(csv_file, index_col=[0])

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

#########################################################
# def aligner(full_df, target_sequence, align_method, filename, output_path, lliteral, rliteral, gop=3, gep=1, dir="3p", rev_complement="False"):
align_class = {"align_local": align_local,
        "align_local2": align_local2,
        "align_local3":align_local3,
        "align_global":align_global,
        "align_global2":align_global2}  
align_method = align_local3
id_f=1
aligner_init = align_class.get(str(align_method), None)  # Get the chosen class, or None if input is bad

#########################################################
align_method= "align_local2"
rev_complement="True"
gop=3
gep=1
dir="5p"

full_df.iloc[,0]
old_chars = "ACGT"
replace_chars = "TGCA"
tab = str.maketrans(old_chars,replace_chars)
headers=[]
aligned_seqs=[]
target_sequence
for seq_i in range(len(full_df.iloc[:,-1])):
    #yield iteratively the header of the certain seq and the corresponding seq
    #header=">"+ str(id_f)+"_CluSeq:" + str((round(full_df.iloc[seq_i,-4],5))) + "_var:"+str((round(full_df.iloc[seq_i,-2],5))) +"_sd:" + str((round(full_df.iloc[seq_i,-1],5)))
    header="CluSeq:" + str((round(full_df.iloc[seq_i,-2], 5))) + "_sd:" + str((round(full_df.iloc[seq_i,-1],5)))
    #seq_obj_align = aligner_init(full_df.iloc[seq_i,0], target_sequence, gop, gep).align()
    seq_obj_align = aligner_init(str(full_df.iloc[seq_i,0]), target_sequence, gop, gep).align()
    seq_obj_align = re.sub(r'[(\d|\s]', '', seq_obj_align) #remove digits from the string caused by the alignment and empty spaces from the start
    print(f'seq_obj_normal: {seq_obj_align}')
    headers.append(header)
    if rev_complement=="True":
        aligned_seqs.append(str(seq_obj_align).translate(tab)[::-1])
        print(f'seq_obj_revcompl: {str(seq_obj_align).translate(tab)[::-1]}')
    else:
        aligned_seqs.append(seq_obj_align)
#########################################################
a = Seq("TTGCTCACCATGGTGGCGCGCCTGT")
b = Seq("GAGGCGAAGAACTCCTCTTAGGCGGCAGC")
a.reverse_complement()
b.reverse_complement()

GCTGCCGCCTAAGAGGAGTTCTTCGCCTC

ACAGGCGCGCCACCATGGTGAGCAA
ACAGGCGCGCCACCATGGTGAGCAA
# translate_NT(result, corr_frame,direc, out_csv)
def add_primers_save2(aligned_data_trim,filename, target_sequence,lliteral, dir, rev_complement):
    id_f = 1
    ref="Ref"
    old_chars = "ACGT"
    replace_chars = "TGCA"
    tab = str.maketrans(old_chars,replace_chars)
    if rev_complement=="True":
        target_sequence = str(target_sequence).translate(tab)[::-1]
    else:
        target_sequence=target_sequence
    with open(filename, "w") as handle:
        if '3' in dir:
            whole_ref=lliteral.split("=")[1] + "-" + target_sequence
            seq_obj = SeqRecord(Seq(whole_ref), id=str(0), description=ref)
            count = SeqIO.write(seq_obj, handle, "fasta")
            for id, seq_prims in aligned_data_trim.items():
                #descr=id.split("_")[0] + "_" + id.split("_")[1]
                descr=id
                whole_seq=lliteral.split("=")[1] + "-" +seq_prims
                seq_obj = SeqRecord(Seq(whole_seq), id=str(id_f), description=descr)
                print(seq_obj)
                count = SeqIO.write(seq_obj, handle, "fasta")
                id_f+=1
        else:
            whole_ref=target_sequence + "-" + lliteral.split("=")[1].translate(tab)[::-1] # reverse primer position as the 5p seq is sequenced in rev
            # whole_ref=target_sequence + "-" + lliteral.split("=")[1](tab)[::-1] # reverse primer position as the 5p seq is sequenced in rev
            seq_obj = SeqRecord(Seq(whole_ref), id=str(0), description=ref)
            count = SeqIO.write(seq_obj, handle, "fasta")
            for id, seq_prims in aligned_data_trim.items():
                #descr=id.split("_")[0] + "_" + id.split("_")[1]
                descr=id
                whole_seq=seq_prims+ "-" + lliteral.split("=")[1].translate(tab)[::-1]
                seq_obj = SeqRecord(Seq(whole_seq), id=str(id_f), description=descr)
                print(seq_obj)
                count = SeqIO.write(seq_obj, handle, "fasta")
                id_f+=1


    print("Saved!")


class gen_aligned_data():
    # generate alignments by aligning each seq to the rev, either compl or not
    def __init__(self, full_df, aligner_init, target_sequence, gop, gep, rev_compl):
        self.full_df=full_df
        self.aligner_init=aligner_init
        self.target_sequence=target_sequence
        self.gop=gop
        self.gep=gep
        self.rev_compl=rev_compl    

    def align_data(self):
        old_chars = "ACGT"
        replace_chars = "TGCA"
        tab = str.maketrans(old_chars,replace_chars)
        headers=[]
        aligned_seqs=[]
        for seq_i in range(len(self.full_df.iloc[:,-1])):
            #yield iteratively the header of the certain seq and the corresponding seq
            #header=">"+ str(id_f)+"_CluSeq:" + str((round(full_df.iloc[seq_i,-4],5))) + "_var:"+str((round(full_df.iloc[seq_i,-2],5))) +"_sd:" + str((round(full_df.iloc[seq_i,-1],5)))
            header="CluSeq:" + str((round(self.full_df.iloc[seq_i,-2], 5))) + "_sd:" + str((round(self.full_df.iloc[seq_i,-1],5)))
            #seq_obj_align = aligner_init(full_df.iloc[seq_i,0], target_sequence, gop, gep).align()
            seq_obj_align = self.aligner_init(str(self.full_df.iloc[seq_i,0]), self.target_sequence, self.gop, self.gep).align()
            seq_obj_align = re.sub(r'[(\d|\s]', '', seq_obj_align) #remove digits from the string caused by the alignment and empty spaces from the start
            print(f'seq_obj_normal: {seq_obj_align}')
            headers.append(header)
            if self.rev_compl=="True":
                aligned_seqs.append(str(seq_obj_align).translate(tab)[::-1])
                print(f'seq_obj_revcompl: {str(seq_obj_align).translate(tab)[::-1]}')
            else:
                aligned_seqs.append(seq_obj_align)
        return(dict(zip(headers, aligned_seqs)))

align_method= "align_local2"
rev_complement="True"
gop=3
gep=1
dir="5p"

# def aligner(full_df, target_sequence, align_method, filename, output_path, lliteral, rliteral, gop=3, gep=1, dir="3p", rev_complement="False"):
align_class = {"align_local": align_local,
        "align_local2": align_local2,
        "align_local3":align_local3,
        "align_global":align_global,
        "align_global2":align_global2}  
id_f=1
aligner_init = align_class.get(str(align_method), None)  # Get the chosen class, or None if input is bad
generate_alignments = gen_aligned_data(full_df, aligner_init, target_sequence, gop, gep, rev_complement) # generate alignments by aligning each seq to the rev, either compl or not

aligned_data=generate_alignments.align_data()

list(aligned_data.keys())[0]
aligned_data['CluSeq:0.90016_sd:0.02255']

#align all the data, save into dict, then ensure that all the seqs are same length (take the longest seq). IF not, then add padding!
aligned_data_trim=align_trimmer(aligned_data, target_sequence)
#data_trim_nodupl=find_remove_duplicates(aligned_data_trim)
aligned_data_trim=reorganise_perc(aligned_data_trim)
#Add primers to both ends of the seq and save
#write_align(data_trim_nodupl, filename, target_sequence)

add_primers_save2(aligned_data_trim, filename, target_sequence, lliteral, dir, rev_complement)
#Generate a visual alignment file using mview
mview_file=output_path +filename.split("/")[-1].split(".")[-2] + ".html"
mview_command='/media/data/AtteR/Attes_bin/mview -in fasta -html head -css on -reference 1 -coloring identity ' + filename + '>' + mview_file
#mview_command='/media/data/AtteR/Attes_bin/mview -in fasta -html head -css on -reference 1 -coloring identity ' + filename + '>' + mview_file
call([mview_command], shell=True)
print("html file created as "+ mview_file)
return(aligned_data_trim)
#os.system('/media/data/AtteR/Attes_bin/mview -in fasta -html head -css on -coloring any {} > {}'.format(str(filename), str(mview_file))) 
#subprocess.run(['/media/data/AtteR/Attes_bin/mview', '-in fasta', '-html head', '-css on', '-coloring any', filename, '>', mview_file])
