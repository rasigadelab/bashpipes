import os
import re
import argparse

# Welcome to creation of replicons.tsv file !
# To run it, type in a terminal:
# cd /location/to/script
# python3 find_replicons.py -d path/to/data/folder

## FUNCTIONS ##

def main(project_dir, output_dir, chr_only):
    #Step1 - output file creation
    output_name = os.path.join(output_dir,"replicons.tsv")
    output_file=open(output_name, "w")
    output_file.write("fasta_file\tSample\treplicon\n")
    #Step2- which input information do we want
    fasta_file = "" #path of file containing assembly of the replicon
    sample = "" #name of sample
    replicon = "" #name of replicons
    #Step3- Searching assemblies in genomes folder
    replicons_list = {}
    dir_genomes = os.path.join(project_dir, "genomes")
    for sample_repo in os.listdir(dir_genomes):
        sample = sample_repo
        mob_recon_dir = os.path.join(dir_genomes, sample_repo, "mob_recon")
        mlst_dir = os.path.join(dir_genomes, sample_repo, "mlst")
        species = ""
        #Step3A- Checking sample specie
        for mlst_file in os.listdir(mlst_dir):
            if mlst_file == "mlst.tsv":
                mlst_path = os.path.join(mlst_dir, mlst_file)
                f = open(mlst_path, "r")
                species = f.readline().split('\t')[1]
                if species == "cronobacter":
                    species = "ecloacae"
        #Step3B- Checking which replicons are in the sample
        for out_file in os.listdir(mob_recon_dir):
            # Option to only analize chromosomes
            if chr_only:
                pattern = "chromosome.fasta"
            else:
                pattern = ".fasta"

            if pattern in out_file:
                fasta_file = os.path.join(mob_recon_dir, out_file)
                replicon = out_file.strip('.fasta')
                if replicon == "chromosome" and species != "-":
                    replicon = species
                if replicon not in replicons_list.keys():
                    replicons_list[replicon] = []
                replicons_list[replicon].append((sample, fasta_file))
    #Step4- Only writing replicons with at least 2 samples (because needed to compare at least 2 genomes)
    for rep in replicons_list:
        if len(replicons_list[rep]) >= 2:
            for sample_item in replicons_list[rep]:
                output_file.write(sample_item[1]+"\t"+sample_item[0]+"\t"+rep+"\n")
    output_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Creation of replicons.tsv')
    parser.add_argument("-d", dest="path_to_data", required=True, help="Path to the folder containing the genomes repository")
    parser.add_argument("-o", dest="output_dir", required=False, help="Path to output directory")
    parser.add_argument("--only-chromosome", dest="chr_only", required=False, default=False, action="store_true", help="Option to only take into account chromosomes")
    args = parser.parse_args()

    if args.output_dir is None:
        args.output_dir = args.path_to_data

    path_to_data = args.path_to_data
    output_dir = args.output_dir
    chr_only = args.chr_only

    main(path_to_data, output_dir, chr_only)

# THE END
