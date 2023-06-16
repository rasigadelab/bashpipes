import os
import re
import argparse

#Welcome to creation of Files_location.tsv file !
#Note: this is a python script, don't be scared
#To run it, type in a terminal:
# cd /location/to/script
# python3 files_location.py -d path/to/data/folder

## FUNCTIONS ##
def remove_spaces(filename):
    # Example case: "Sample_S18_R2.fastq (2).gz"
    # Change filename to differentiate it from the other file with same name
    file_prefix = filename.split('.')[0]+"_2"
    filename = '.'.join([file_prefix,"fastq.gz"])
    return filename

def write_ONT_reads(subrepo, abs_path_dir, output_file):
    sample_type = "ONT"
    sample = subrepo
    sample_dir_path = os.path.join(abs_path_dir, subrepo)
    for nano_file in os.listdir(sample_dir_path):
        full_path = os.path.join(sample_dir_path, nano_file)
        output_file.write(full_path+"\t"+sample+"\t"+sample_type+"\n")
    return sample

def write_Illumina_reads(subrepo, abs_path_dir, output_file):
    full_path = os.path.join(abs_path_dir,subrepo)
    # Remove whitespaces in filename
    if " " in subrepo:
        subrepo = remove_spaces(subrepo)
        new_path = os.path.join(abs_path_dir,subrepo)
        os.rename(full_path, new_path)
        full_path = new_path
    # Get characters before first underscore to get sample name
    sample = subrepo.split('_')[0]
    # Remove 'bis' after sample name
    sample = re.sub(r'\D+$', '', sample)
    # Remove '-bis' after sample name
    if sample.count("-") > 1:
        sample = sample.rsplit('-',1)[0]
    sample_type = subrepo.split('_')[-1].split('.')[0]
    # While R1/R2 info is not catched, continue to search in string where is it
    if sample_type != "R1" and sample_type != "R2":
        for i in range(2,len(subrepo.split('_'))):
            sample_type = subrepo.split('_')[-i]
            if sample_type == "R1" or sample_type == "R2":
                break 
    sample = sample.replace("Epitrack", "Epi")
    output_file.write(full_path+"\t"+sample+"\t"+sample_type+"\n")
    return sample

def search_reads(dir_path, method, output_file):
    samples_list = []
    dir_samples = os.listdir(dir_path)
    for repo in dir_samples:
        abs_path_dir = os.path.join('.', dir_path, repo)
        subdir_samples = os.listdir(abs_path_dir)
        for subrepo in subdir_samples:
            if method == "Nanopore":
                sample = write_ONT_reads(subrepo, abs_path_dir, output_file)
            elif method == "Illumina":
                sample = write_Illumina_reads(subrepo, abs_path_dir, output_file)
            if sample not in samples_list:
                samples_list.append(sample)
    return samples_list

def main(project_dir, output_dir, illumina_only):
    #Step1 - output file creation
    output_name = os.path.join(output_dir,"Files_location.tsv")
    output_file=open(output_name, "w")
    output_file.write("full_path\tSample\ttype\n")
    #Step2- which input information do we want
    full_path = "" #path of file containing reads
    sample = "" #name of sample
    sample_type = "" #could be ONT-R1-R2
    #Step3- Illumina reads
    path_to_illumina = os.path.join(project_dir,"Illumina")
    samples_illumina = search_reads(path_to_illumina, "Illumina", output_file)
    #Step4- Nanopore reads
    path_to_nanopore = os.path.join(project_dir,"Nanopore")
    samples_nanopore = search_reads(path_to_nanopore, "Nanopore", output_file)
    output_file.close()
    #Step5- Check for missing file
    if not illumina_only:
        missing_illumina = list(set(samples_nanopore) - set(samples_illumina))
        with open(os.path.join(output_dir,'missing_illumina_data.txt'), 'w') as f1:
            f1.write('\n'.join(map(str, sorted(missing_illumina))))
        missing_nanopore = list(set(samples_illumina) - set(samples_nanopore))
        with open(os.path.join(output_dir,'missing_nanopore_data.txt'), 'w') as f2:
            f2.write('\n'.join(map(str, sorted(missing_nanopore))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Creation of Files_location.tsv')
    parser.add_argument("-d", dest="path_to_data", required=True, help="Path to the folder containing Illumina and Nanopore data")
    parser.add_argument("-o", dest="output_dir", required=False, help="Path to output directory")
    parser.add_argument("--illumina-only", dest="illumina_only", required=False, default=False, action="store_true", help="Data type to analyze (boolean)")
    args = parser.parse_args()

    if args.output_dir is None:
        args.output_dir = args.path_to_data

    path_to_data = args.path_to_data
    output_dir = args.output_dir
    illumina_only = args.illumina_only

    main(path_to_data, output_dir, illumina_only)

# THE END
