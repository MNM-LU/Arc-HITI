# Arc-HITI
User should have starcode, cutadapt, bbduk and mview installed on their computer and have them in the path variables. If mview does not work,
the amino acid alignments can be visualised by taking the fasta file inside specific sample's aligned/AA/fasta/ directory and uploading them to
https://alignmentviewer.org/. 

Inside Arc-HiTi folder one can find samples characterised by folder name starting
with HITI. Each HITI sample subfolder contains main.py python script for the analysis workflow of the given sample. The script calls on functions from scripts_main.py found in the root.The folder amplicon_files contains DNA sequences of the reference files with the whole fluorescent protein and the Arc gene along with information about the primer binding sites and primer templates. In the root folder there is also the following files:

**scripts_main.py** - Contains all the functions and classes used for the
analysis. The functions are called in the analysis workflows performed
individually in the main.py script found within each sample subfolder.
**conda_env.yml** - The required python packages and their versions. to install,
run the following command:
> conda env create -f conda_env.yml


