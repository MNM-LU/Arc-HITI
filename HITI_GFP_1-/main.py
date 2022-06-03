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
sample_dir=scripts_dir + "/HITI_GFP_1-/"
os.chdir(sample_dir)
#path to the analysis folder
base_path = '/media/data/AtteR/projects/hiti/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
export_path=sample_dir + "trimmed_data/"

# #where the programs bbduk and starcode are found
# program_path="/media/data/AtteR/Attes_bin/"

#############

#############
#edited script for read preprocessing and clustering:
def trimRead_hiti(animal_nr,base_path,transgene,filterlitteral,lliteral,rliteral,read_fwd,direc):
    animal_nr = str(animal_nr)
    "Filters and trims the reads"
    search_path = base_path+animal_nr+'*'+transgene+'*'+direc+'*/'
    
    animal_p5_cat = tempfile.NamedTemporaryFile(suffix = '.fastq.gz').name
    animal_p7_cat = tempfile.NamedTemporaryFile(suffix = '.fastq.gz').name
    test_file_p5_out = tempfile.NamedTemporaryFile(suffix = '.fastq').name
    test_file_p7_out = tempfile.NamedTemporaryFile(suffix = '.fastq').name
    test_file_p5_filter = tempfile.NamedTemporaryFile(suffix = '.fastq').name
    
    if read_fwd:
        animal_p5 = glob.glob(search_path+'*R1*')
        animal_p7 = glob.glob(search_path+'*R2*')
    else:
        animal_p5 = glob.glob(search_path+'*R2*')
        animal_p7 = glob.glob(search_path+'*R1*')
    
    stats_out = "trim_data/"+animal_nr+ "_" + direc +'_stats-filter.txt'

    cat_p5= "cat "+" ".join(animal_p5)+" > "+animal_p5_cat
    call([cat_p5], shell=True)
    cat_p7= "cat "+" ".join(animal_p7)+" > "+animal_p7_cat
    call([cat_p7], shell=True)

    #stats_out = export_path+animal_nr+'_'+transgene+'_'+assay_end+'_stats-filter.txt'
    
    kmer = '20'
    hdist = '3'
    param=" k="+kmer+" hdist="+hdist+" rcomp=f skipr2=t threads=32 overwrite=true"
    
    call_sequence = "bbduk.sh in="+animal_p5_cat +" in2="+animal_p7_cat+" outm1="+ test_file_p5_out +" outm2="+test_file_p7_out+" literal="+filterlitteral+" stats="+stats_out + param
    call([call_sequence], shell=True)

    call_sequence = "bbduk.sh in="+test_file_p5_out+" out="+test_file_p5_filter+ " literal=AAAAAAAAA,CCCCCCCCC,GGGGGGGGG,TTTTTTTTT k=9 mm=f overwrite=true minlength=40"
    call([call_sequence], shell=True)
    test_file_p5_filter2 = tempfile.NamedTemporaryFile(suffix = '.fastq').name #when cutadapt applied

    cutadapt_call="cutadapt -g "+lliteral+ " --discard-untrimmed -o " + test_file_p5_filter2 + " " + test_file_p5_filter
    call([cutadapt_call], shell=True)
    # cutadapt_call="cutadapt -a "+rliteral+" -o " + test_file_p5_filter2 + " " + test_file_p5_filter2
    # call([cutadapt_call], shell=True)

    test_file_p5_out_starcode = tempfile.NamedTemporaryFile(suffix = '.tsv').name
    starcode_call= "starcode -i "+test_file_p5_filter2+" -t 32 -o "+test_file_p5_out_starcode
    call([starcode_call], shell=True)

    df=pd.read_csv(test_file_p5_out_starcode, sep='\t', header=None)
    # result="unaligned/Starcode_HITI_GFP_3p_" + animal_nr + "_.csv"
    # df.to_csv(result)

    df = df.rename(columns={0: 'sequence', 1:'count'})
    df = df.rename(columns={'count':animal_nr+'_count',})
    
    return df

def analyze_all(base_path, transgene, filterlitteral,lliteral,rliteral,export_path,read_fwd,animal_list, target_sequence, direc):
    complete_df = pd.DataFrame({'sequence': [target_sequence]})
    for animal in animal_list:
        df_this = trimRead_hiti(animal,base_path,transgene,filterlitteral,lliteral,rliteral,read_fwd,direc)
        complete_df = pd.merge(complete_df, df_this, on="sequence", how='outer')
    
    complete_df = complete_df.fillna(value=0)
    #complete_df['percent_sum'] = complete_df[perc_cols].sum(axis=1)
    #export_csv = export_path+transgene+'_'+assay_end+'.csv'
    #complete_df.to_csv(export_csv, index=False)
    print("Done!")
    return complete_df

#saves file as fasta and csv




#############
transgene = 'GFP'
assay_end = '3p'
animal_list = [20, 21, 22, 23, 24] 
read_fwd = True
filterlitteral='CAGCTCCATCTGGTCGTCGGT'
filterlitteral=filterlitteral.upper()
lliteral = ' literal=GTGGTCATATGGTCCAGCTCC'

target_sequence='ATCTGGTCGTCGGTGCTGCGGCTCCGcggagccgcagcaccgaCCTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGAACTCCAGCAGGACCATGTGATCGCGCTTCTCGTTGGGGTCTTTGCTCAGGGCGGACTGGGTGCTCAGGTAGTGGTTGTCGGGCAGCAGCACGGGGCCGTCGCCGATGGGGGTGTTCTGCTGGTAGTGGTCGGCGAGCTGCACGCTGCCGTCCTCGATGTTGTGGCGGATCTTGAAGTTCACCTTGATGCCGTTCTTCTGCTTGTCGGCCATGATATAGACGTTGTGGCTGTTGTAGTTGTACTCCAGCTTGTGCCCCAGGATGTTGCCGTCCTCCTTGAAGTCGATGCCCTTCAGCTCGATGCGGTTCACCAGGGTGTCGCCCTCGAACTTCACCTCGGCGCGGG'
filterlitteral='CAGCTCCATCTGGTCGTCGGT'
rliteral=''
target_sequence='ATCTGGTCGTCGGTGCTGCGGCTCCGcggagccgcagcaccgaCCTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGAACTCCAGCAGGACCATGTGATCGCGCTTCTCGTTGGGGTCTTTGCTCAGGGCGGACTGGGTGCTCAGGTAGTGGTTGTCGGGCAGCAGCACGGGGCCGTCGCCGATGGGGGTGTTCTGCTGGTAGTGGTCGGCGAGCTGCACGCTGCCGTCCTCGATGTTGTGGCGGATCTTGAAGTTCACCTTGATGCCGTTCTTCTGCTTGTCGGCCATGATATAGACGTTGTGGCTGTTGTAGTTGTACTCCAGCTTGTGCCCCAGGATGTTGCCGTCCTCCTTGAAGTCGATGCCCTTCAGCTCGATGCGGTTCACCAGGGTGTCGCCCTCGAACTTCACCTCGGCGCGGG'
target_sequence=target_sequence.upper()
direc="3p"
####################

#read preprocessing for each sample: trim, record read counts before and after trimming, cluster the reads 
#using starcode, calculate percentage, merge into the full dataframe containing read count-and percentage for
#each sample.
full_df=analyze_all(base_path, transgene, filterlitteral,lliteral,rliteral,export_path,read_fwd, animal_list, target_sequence, direc)
full_df_trim=calculate_perc_sd(full_df, 3)
result="unaligned/HITI_GFP_3p_1-.fasta"
csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
full_df_trim.to_csv(csv_file)
save_fasta(result, full_df_trim, target_sequence)
#full_df_trim=pd.read_csv(csv_file,  index_col=[0])

####################
#NT
#Perform alignments
####################
output_path="aligned/NT/"
result=output_path + "HITI_GFP_3p_1-_prim.fasta"
aligner(full_df_trim, target_sequence, "align_local2", result, output_path, lliteral, rliteral, 3,1)
####################


#AA
####################
corr_frame=2
result="unaligned/HITI_GFP_3p_1-.fasta"
out_csv="aligned/AA/HITI_GFP_3p_1-_AA.csv"
output_html="aligned/AA/HITI_GFP_3p_1-_AA.html"
translate_NT(result, corr_frame,direc, out_csv)

####################



def trimRead_hiti(animal_nr,base_path,transgene,filterlitteral,lliteral,rliteral,read_fwd,direc):
    animal_nr = str(animal_nr)
    "Filters and trims the reads"
    search_path = base_path+animal_nr+'*'+transgene+'*'+direc+'*/'
    
    animal_p5_cat = tempfile.NamedTemporaryFile(suffix = '.fastq.gz').name
    animal_p7_cat = tempfile.NamedTemporaryFile(suffix = '.fastq.gz').name
    test_file_p5_out = tempfile.NamedTemporaryFile(suffix = '.fastq').name
    test_file_p7_out = tempfile.NamedTemporaryFile(suffix = '.fastq').name
    test_file_p5_filter = tempfile.NamedTemporaryFile(suffix = '.fastq').name
    
    if read_fwd:
        animal_p5 = glob.glob(search_path+'*R1*')
        animal_p7 = glob.glob(search_path+'*R2*')
    else:
        animal_p5 = glob.glob(search_path+'*R2*')
        animal_p7 = glob.glob(search_path+'*R1*')
    
    stats_out = "trim_data/"+animal_nr+ "_" + direc +'_stats-filter.txt'

    cat_p5= "cat "+" ".join(animal_p5)+" > "+animal_p5_cat
    call([cat_p5], shell=True)
    cat_p7= "cat "+" ".join(animal_p7)+" > "+animal_p7_cat
    call([cat_p7], shell=True)

    #stats_out = export_path+animal_nr+'_'+transgene+'_'+assay_end+'_stats-filter.txt'
    
    kmer = '20'
    hdist = '3'
    param=" k="+kmer+" hdist="+hdist+" rcomp=f skipr2=t threads=32 overwrite=true"
    
    call_sequence = "bbduk.sh in="+animal_p7_cat+" in2="+animal_p5_cat+" outm1="+test_file_p7_out+" outm2="+test_file_p5_out+" literal="+filterlitteral+" stats="+stats_out + param
    call([call_sequence], shell=True)

    call_sequence = "bbduk.sh in="+test_file_p5_out+" out="+test_file_p5_filter+ " literal=AAAAAAAAA,CCCCCCCCC,GGGGGGGGG,TTTTTTTTT k=9 mm=f overwrite=true minlength=40"
    call([call_sequence], shell=True)
    test_file_p5_filter2 = tempfile.NamedTemporaryFile(suffix = '.fastq').name #when cutadapt applied

    cutadapt_call="cutadapt -g "+lliteral+ " --discard-untrimmed -o " + test_file_p5_filter2 + " " + test_file_p5_filter
    call([cutadapt_call], shell=True)
    # cutadapt_call="cutadapt -a "+rliteral+" -o " + test_file_p5_filter2 + " " + test_file_p5_filter2
    # call([cutadapt_call], shell=True)

    test_file_p5_out_starcode = tempfile.NamedTemporaryFile(suffix = '.tsv').name
    starcode_call= "starcode -i "+test_file_p5_filter2+" -t 32 -o "+test_file_p5_out_starcode
    call([starcode_call], shell=True)

    df=pd.read_csv(test_file_p5_out_starcode, sep='\t', header=None)
    # result="unaligned/Starcode_HITI_GFP_3p_" + animal_nr + "_.csv"
    # df.to_csv(result)

    df = df.rename(columns={0: 'sequence', 1:'count'})
    df = df.rename(columns={'count':animal_nr+'_count',})
    
    return df

#GFP 5p
#############
transgene = 'GFP'
read_fwd = True
direc="5p"
#(2561-2210)%3
animal_list = [13, 14, 15, 16, 17, 18] 
# filterlitteral = 'GCGCGGGTCTTGTAGTTGCCGTCGTCCTTGAAGAAGATGGTGCGCTCCTGGACG'

#original vars used
lliteral = ' literal=CCTCAGAGGAGTTCT'
rliteral = ' literal=GGCGAGGAGCTGTT'
base_path = '/media/data/AtteR/projects/hiti/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
export_path = '/media/data/AtteR/projects/hiti/pipeline_output_reorg/'
#target_sequence="GCGCGGGTCTTGTAGTTGCCGTCGTCCTTGAAGAAGATGGTGCGCTCCTGGACGTAGCCTTCGGGCATGGCGGACTTGAAGAAGTCGTGCTGCTTCATGTGGTCGGGGTAGCGGCTGAAGCACTGCACGCCGTAGGTCAGGGTGGTCACGAGGGTGGGCCAGGGCACGGGCAGCTTGCCGGTGGTGCAGATGAACTTCAGGGTCAGCTTGCCGTAGGTGGCATCGCCCTCGCCCTCGCCGGACACGCTGAACTTGTGGCCGTTTACGTCGCCGTCCAGCTCGACCAGGATGGGCACCACCCCGGTGAACAGCTCCTCGCCCTTGCTCACCATGGTGGCGCGcctgttAACAGGCTA"
target_sequence="TTAGCTTCTGCCTCAGAGGAGTTCTTAGCCTGTTAACAGGCGCGCCACCATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGC"
target_sequence=target_sequence.upper()
filterlitteral = 'GCGCGGGTCTTGTAGTTGCCGTCGTCCTTGAAGAAGATGGTGCGCTCCTGGACG'
lliteral = ' literal=CCTCAGAGGAGTTCT'
rliteral = ' literal=GGCGAGGAGCTGTT'

target_sequence = "TAGCCTGTTaacaggCGCGCCACCATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGC"
#read preprocessing for each sample: trim, record read counts before and after trimming, cluster the reads 
#using starcode, calculate percentage, merge into the full dataframe containing read count-and percentage for
#each sample.
full_df=analyze_all(base_path, transgene, filterlitteral,lliteral,rliteral,export_path,read_fwd, animal_list, target_sequence, direc)

full_df_trim=calculate_perc_sd(full_df,3)

result="unaligned/HITI_GFP_5p_1-.fasta"
csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
full_df_trim.to_csv(csv_file)
save_fasta(result, full_df_trim, target_sequence)
#full_df_trim=pd.read_csv(csv_file,  index_col=[0])
#########


#NT
####################
output_path="aligned/NT/"
result=output_path+"HITI_GFP_5p_1-_prim.fasta"
test_res=aligner(full_df_trim, target_sequence, "align_local2", result, output_path, lliteral, rliteral,3,1)
####################

#AA
####################
# corr_frame=1
# result="unaligned/HITI_GFP_5p_1-.fasta"
# out_csv="aligned/AA/HITI_GFP_5p_1-_AA.csv"
# output_html="aligned/AA/HITI_GFP_5p_1-_AA.html"

# translate_NT(result, corr_frame,direc, out_csv)
#############

