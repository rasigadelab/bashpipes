#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Title: files_location.py
# Description: Outputs a list containing name, type of read and path to FASTQ files for each sample analyzed.
# Author: Aurélie Fischer
# Research team: PHE3ID, CIRI, Hospices Civils de Lyon
# License: GNU Affero General Public License v3.0 (AGPL-3.0)
# SPDX-License-Identifier: AGPL-3.0-only
#
# Copyright (C) 2026 Aurélie Fischer

import os
import re
import argparse

# Welcome to creation of Files_location.tsv file !
# To run it, type in a terminal:
# cd /location/to/script
# python3 files_location.py -d path/to/data/folder [--illumina-only]

## FUNCTIONS ##
def remove_spaces(filename):
    """
    Change filename of reads from same sample with a duplicated name.

    Args:
        filename (str): Name of the file to correct.

    Returns:
        New filename with the replacement "_2" instead of " (2)".
    """
    # Example case: "Sample_S18_R2.fastq (2).gz" --> "Sample_S18_R2_2.fastq.gz"
    file_prefix = filename.split('.')[0]+"_2"
    filename = '.'.join([file_prefix,"fastq.gz"])
    return filename

def write_ONT_reads(subrepo, abs_path_dir, output_file):
    """
    Search for filepath of ONT reads of a sample and write it to TSV file.

    Args:
        subrepo (str): Repository containing ONT reads of sample.
        abs_path_dir (str): Absolute path of parent directory containing ONT sample reads.
        output_file (str): Path to output TSV file containing sample information.

    Returns:
        Repository containing ONT reads of sample (also name of processed sample).
    """
    sample_type = "ONT"
    sample = subrepo # example: "1-1"
    sample_dir_path = os.path.join(abs_path_dir, subrepo) # example: "/data/Nanopore/Batch1/1-1"
    # For each file inside sample ONT reads directory 
    for nano_file in os.listdir(sample_dir_path):
        # Write path to reads inside output file
        full_path = os.path.join(sample_dir_path, nano_file)
        output_file.write(full_path+"\t"+sample+"\t"+sample_type+"\n") # example: "/data/Nanopore/Batch1/1-1/1-1_azejbs.fastq.gz 1-1 ONT"
    return sample

def write_Illumina_reads(subrepo, abs_path_dir, output_file):
    """
    Search for filepath of Illumina reads of a sample and write it to TSV file.

    Args:
        subrepo (str): File containing Illumina reads of sample.
        abs_path_dir (str): Absolute path of parent directory containing Illumina sample reads.
        output_file (str): Path to output TSV file containing sample information.

    Returns:
        Repository containing Illumina reads of sample (also name of processed sample).
    """
    full_path = os.path.join(abs_path_dir,subrepo) # example: "/data/Illumina/Run/Epi-1_S18_R1.fastq.gz"
    # Remove whitespaces in filename
    if " " in subrepo:
        subrepo = remove_spaces(subrepo)
        new_path = os.path.join(abs_path_dir,subrepo)
        os.rename(full_path, new_path)
        full_path = new_path
    # Get characters before first underscore to get sample name 
    sample = subrepo.split('_')[0] # example: "Epi-1"
    # Remove 'bis' after sample name (Epi-1bis_S18_R1.fastq.gz --> get only Epi-1)
    sample = re.sub(r'\D+$', '', sample)
    # Remove '-bis' after sample name (Epi-1-bis_S18_R1.fastq.gz --> get only Epi-1)
    if sample.count("-") > 1:
        sample = sample.rsplit('-',1)[0]
    # Get type of reads, either R1 or R2 from filepath (Epi-1_S18_R1.fastq.gz --> get R1)
    sample_type = subrepo.split('_')[-1].split('.')[0]
    # While R1/R2 info is not caught, continue to search in string where is it (Epi-1_S18_R1_L001.fastq.gz --> get R1)
    if sample_type != "R1" and sample_type != "R2":
        for i in range(2,len(subrepo.split('_'))):
            sample_type = subrepo.split('_')[-i]
            if sample_type == "R1" or sample_type == "R2":
                break 
    sample = sample.replace("Epitrack", "Epi")
    # Write read info to output file [/data/Illumina/Run/Epi-1_S18_R1.fastq.gz Epi-1 R1]
    output_file.write(full_path+"\t"+sample+"\t"+sample_type+"\n")
    return sample

def search_reads(dir_path, method, output_file):
    """
    Search for path to reads, based on type of reads, inside a directory.

    Args:
        dir_path (str): Path to directory inside which there could be reads.
        method (str): Type of sequencing, either Nanopore or Illumina.
        output_file (str): Path to output TSV file containing sample information.

    Returns:
        List of samples with reads processed and found inside directory.
    """
    samples_list = []
    dir_samples = os.listdir(dir_path)
    # For each directory inside dir_path
    for repo in dir_samples:
        abs_path_dir = os.path.join('.', dir_path, repo)
        subdir_samples = os.listdir(abs_path_dir)
        for subrepo in subdir_samples:
            if method == "Nanopore":
                # Architecture = [/data/Nanopore/BatchX/Sample/reads.fastq.gz]
                # Here subrepo is a sample directory
                sample = write_ONT_reads(subrepo, abs_path_dir, output_file)
            elif method == "Illumina":
                # Architecture = [/data/Illumina/BatchX/Sample_reads.fastq.gz]
                # Here subrepo is a fastq.gz file
                sample = write_Illumina_reads(subrepo, abs_path_dir, output_file)
            if sample not in samples_list:
                samples_list.append(sample)
    return samples_list

def main(project_dir, output_dir, illumina_only):
    # Step1 - output file creation
    output_name = os.path.join(output_dir,"Files_location.tsv")
    output_file=open(output_name, "w")
    output_file.write("full_path\tSample\ttype\n")
    # Step2- which input information do we want
    full_path = "" # path of file containing reads
    sample = "" # name of sample
    sample_type = "" # could be ONT-R1-R2
    # Step3- Illumina reads
    path_to_illumina = os.path.join(project_dir,"Illumina")
    samples_illumina = search_reads(path_to_illumina, "Illumina", output_file)
    # Step4- Nanopore reads
    path_to_nanopore = os.path.join(project_dir,"Nanopore")
    samples_nanopore = search_reads(path_to_nanopore, "Nanopore", output_file)
    output_file.close()
    # Step5- Check for missing file (if there are Nanopore reads but no Illumina reads found for 1 sample, or the opposite)
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
