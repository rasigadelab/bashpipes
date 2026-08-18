#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Title: preparation.py
# Description: Outputs a list containing sample, gene, plasmid id, plasmid inc type and path to fasta for each plasmid analyzed.
# Author: Aurélie Fischer
# Research team: PHE3ID, CIRI, Hospices Civils de Lyon
# License: GNU Affero General Public License v3.0 (AGPL-3.0)
# SPDX-License-Identifier: AGPL-3.0-only
#
# Copyright (C) 2026 Aurélie Fischer

import os
import argparse
import pandas as pd

# Welcome to creation of plasmid_locations.tsv file !
# To run it, type in a terminal:
# cd /location/to/script
# python3 preparation.py -d path/to/data/folder

## FUNCTIONS ##
def plasmid_db(dir_path):
    """
    Gather information around plasmids carried by analyzed genomes.

    Args:
        dir_path (str): Path to the directory containing genomes folder.

    Returns:
        None
    """
    sample_plasmid_global = pd.DataFrame()
    # For each analysed genomes
    for sample in os.listdir(os.path.join(dir_path, "genomes")):
        # Read amrfinder.tsv, get the gene content.
        amrfinder_path = os.path.join(dir_path, "genomes", sample, "amrfinder/amrfinder.tsv")
        # Filter name of gene + name of contig.
        amrfinder_res = pd.read_csv(amrfinder_path, sep = "\t", header = 0, usecols = ['Contig id', 'Gene symbol'])
        # Read mob_suite contig_report.txt
        contig_path = os.path.join(dir_path, "genomes", sample, "mob_recon/contig_report.txt")
        # Get : contig and plasmid name associated
        contig_res = pd.read_csv(contig_path, sep = "\t", header = 0, usecols = ['molecule_type', 'primary_cluster_id', 'contig_id'])
        # Filter out chromosomal contigs
        contig_res = contig_res[contig_res['molecule_type'] != "chromosome"]
        # Build plasmid name 
        contig_res = contig_res.drop(columns=['molecule_type'])
        # Read mob_suite mobtyper_results.txt
        typer_path = os.path.join(dir_path, "genomes", sample, "mob_recon/mobtyper_results.txt")
        # Get plasmid name and plasmid inc type.
        typer_res = pd.read_csv(typer_path, sep = "\t", header = 0, usecols = ['primary_cluster_id', 'rep_type(s)'])
        # Merge all info in 1 big table [sample - plasmid name - gene - plasmid_inc]
        sample_plasmid_info = pd.merge(contig_res, amrfinder_res, left_on=['contig_id'], right_on=['Contig id'])
        sample_plasmid_info = pd.merge(sample_plasmid_info, typer_res, on=['primary_cluster_id'])
        # Filter only for genes of interest : blaOXA-48 ; blaVIM ; blaNDM
        sample_plasmid_info = sample_plasmid_info[
                                    (sample_plasmid_info['Gene symbol'] == "blaOXA-48") |
                                    (sample_plasmid_info['Gene symbol'].str.contains('blaVIM', na=False)) |
                                    (sample_plasmid_info['Gene symbol'].str.contains('blaNDM', na=False))]
        # Reorganize columns and add a column path to fasta
        sample_plasmid_info = sample_plasmid_info.drop(columns=['Contig id', 'contig_id'])
        sample_plasmid_info['sample'] = sample
        sample_plasmid_info['full_path'] = dir_path + "/genomes/" + sample + "/mob_recon/plasmid_"+sample_plasmid_info['primary_cluster_id'].astype(str)+".fasta"
        # Reordering columns 
        sample_plasmid_info = sample_plasmid_info.loc[:, ['sample', 'Gene symbol', 'primary_cluster_id', 'rep_type(s)', 'full_path']]
	# Replacing "/" by "-" for inc types
        sample_plasmid_info['rep_type(s)'] = sample_plasmid_info['rep_type(s)'].str.replace("/", "-", regex = False)
        # Bind to other results
        sample_plasmid_global = pd.concat([sample_plasmid_global, sample_plasmid_info], axis = 0)
    return sample_plasmid_global

def main(project_dir, output_dir):
    # Step1 - output file creation
    output_name = os.path.join(output_dir,"plasmid_locations.tsv")
    # Step2- search plasmid info inside genome results dir
    plasmid_info = plasmid_db(project_dir)
    # Step3- write to output file
    plasmid_info.to_csv(output_name, sep = "\t", header = False, index = False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Creation of plasmid_locations.tsv')
    parser.add_argument("-d", dest="path_to_data", required=True, help="Path to the folder containing samples genomes analysis")
    parser.add_argument("-o", dest="output_dir", required=False, help="Path to output directory")
    args = parser.parse_args()

    if args.output_dir is None:
        args.output_dir = args.path_to_data

    path_to_data = args.path_to_data
    output_dir = args.output_dir

    main(path_to_data, output_dir)

# THE END
