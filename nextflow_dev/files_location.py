import os
import re
import argparse

#Welcome to creation of Files_location.tsv file !
#Note: this is a python script, don't be scared
#To run it, type in a terminal:
# cd /location/to/script
# python3 files_location.py -d path/to/data/folder

## FUNCTIONS ##
def write_ONT_reads(subrepo, abs_path_dir, output_file):
    sample_type = "ONT"
    sample = subrepo
    sample_dir_path = os.path.join(abs_path_dir, subrepo)
    for nano_file in os.listdir(sample_dir_path):
        full_path = os.path.join(sample_dir_path, nano_file)
        output_file.write(full_path+"\t"+sample+"\t"+sample_type+"\n")

def write_Illumina_reads(subrepo, abs_path_dir, output_file):
    full_path = os.path.join(abs_path_dir,subrepo)
    # Get characters before first underscore to get sample name
    sample = subrepo.split('_')[0]
    # Remove 'bis' after sample name
    sample = re.sub(r'\D+$', '', sample)
    # Remove '-bis' after sample name
    if sample.count("-") > 1:
        sample = sample.rsplit('-',1)[0]
    sample_type = subrepo.split('_')[-1].split('.')[0]
    if sample_type != "R1" and sample_type != "R2":
        sample_type = subrepo.split('_')[-2]
    output_file.write(full_path+"\t"+sample+"\t"+sample_type+"\n")

def search_reads(dir_path, method, output_file):
    dir_samples = os.listdir(dir_path)
    for repo in dir_samples:
        abs_path_dir = os.path.join('.', dir_path, repo)
        subdir_samples = os.listdir(abs_path_dir)
        for subrepo in subdir_samples:
            if method == "Nanopore":
                write_ONT_reads(subrepo, abs_path_dir, output_file)
            elif method == "Illumina":
                write_Illumina_reads(subrepo, abs_path_dir, output_file)

def main(project_dir, output_name):
    #Step1 - output file creation
    output_path = os.path.join(project_dir, output_name)
    output_file=open(output_path, "w")
    output_file.write("full_path\tSample\ttype\n")
    #Step2- which input information do we want
    full_path = "" #path of file containing reads
    sample = "" #name of sample
    sample_type = "" #could be ONT-R1-R2
    #Step3- Illumina reads
    path_to_illumina = os.path.join(project_dir,"Illumina")
    search_reads(path_to_illumina, "Illumina", output_file)
    #Step4- Nanopore reads
    path_to_nanopore = os.path.join(project_dir,"Nanopore")
    search_reads(path_to_nanopore, "Nanopore", output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Creation of Files_location.tsv')
    parser.add_argument("-d", dest="path_to_data", required=True, help="Path to the folder containing Illumina and Nanopore data")
    parser.add_argument("-o", dest="output_name", required=False, default="Files_location.tsv", help="Name of output file (only tsv file extension)")
    args = parser.parse_args()

    path_to_data = args.path_to_data
    output_name = args.output_name

    main(path_to_data, output_name)

# THE END
