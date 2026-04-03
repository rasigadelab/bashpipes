#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Title: circos_plot.py
# Description: Outputs an annotated circular plot of aligned plasmid sequences.
# Author: Aurélie Fischer
# Research team: PHE3ID, CIRI, Hospices Civils de Lyon
# License: GNU Affero General Public License v3.0 (AGPL-3.0)
# SPDX-License-Identifier: AGPL-3.0-only
#
# Copyright (C) 2026 Aurélie Fischer

import re
import os
import argparse
import matplotlib.pyplot as plt # type: ignore


def scale_to_kb(path_to_circos):
    """
    Display correct scale on Circos graph by adjusting ticks number.

    Args:
        path_to_circos (str): Path to circos results directory.

    Returns:
        None
    """
    # Open ticks file to change
    ticks_conf = open(path_to_circos+"/data/ticks.conf", "r")
    # Create a tmp file that will contain all new lines
    tmp_conf = open(path_to_circos+"/data/tmp.conf", "w")
    for line in ticks_conf:
        if "label_multiplier" in line:
            # Multiply number of ticks [scale], 0.001 = 1 tick /1kb
            tmp_conf.write("label_multiplier = 0.001\n")
        else:
            tmp_conf.write(line)
    ticks_conf.close()
    tmp_conf.close()
    # Remove old ticks config file
    os.remove(path_to_circos+"/data/ticks.conf")
    # Replace it with tmp file freshly filled
    os.rename(path_to_circos+"/data/tmp.conf", path_to_circos+"/data/ticks.conf")

def scatter_for_mismatches(path_to_circos, track_corresp):
    """
    Add a scatterplot on histograms presenting mismatches

    Args:
        path_to_circos (str): Path to circos results directory.
        track_corresp (str): Track labels configuration (labels.txt).

    Returns:
        None
    """
    circos_conf = open(path_to_circos+"/circos.conf", "r")
    tmp_conf = open(path_to_circos+"/tmp.conf", "w")
    # First change filling of every histogram to black
    for line in circos_conf:
        tmp_conf.write(line.replace("fill_color = vlyellow", "fill_color = black"))
    circos_conf.close()
    tmp_conf.close()
    os.remove(path_to_circos+"/circos.conf")
    os.rename(path_to_circos+"/tmp.conf", path_to_circos+"/circos.conf")
    # Then, add a scatterplot on histograms so points are clearly indicated. 
    # This must be added between a "plots /plots" structure.
    # Get Epi-name and track combination
    tmp = track_corresp.split("\t")[-1].split(",") # List of [track0=Epi-name, track1=Epi-name, ...]
    epi_track = {name.split("=")[1]: name.split("=")[0] for name in tmp} # Dict of {Epi-name: track0, Epi-name: track1}
    epi_track.pop("ANY") # Remove ANY from list of tracks to be treated
    # For each sample track, create the scatterplot structure to add
    all_scatters = ""
    for sample in epi_track.keys():
        scatter_to_add = "<plot>\ntype = scatter\nthickness = 1\nfile = data/" + sample
        scatter_to_add += ".mismatches.txt\nr0 = eval(sprintf('%.3fr', conf(" + epi_track[sample] + "_pos) - "
        scatter_to_add += " conf(track_width)))\nr1 = eval(sprintf('%.3fr', conf(" + epi_track[sample] + "_pos)))\n"
        scatter_to_add += "stroke_color = black\nstroke_thickness = 2\nglyph_size = 12\n</plot>"
        all_scatters += scatter_to_add + "\n"
    circos_conf = open(path_to_circos+"/circos.conf", "r")
    tmp_conf = open(path_to_circos+"/tmp.conf", "w")
    # Add all scatterplots just before ending plots declaration in config file
    for line in circos_conf:
        if "</plots>" in line:
            tmp_conf.write(all_scatters)
        tmp_conf.write(line)
    circos_conf.close()
    tmp_conf.close()
    # Replace old config file by new configuration 
    os.remove(path_to_circos+"/circos.conf")
    os.rename(path_to_circos+"/tmp.conf", path_to_circos+"/circos.conf")


def generate_color_range(palette_name, n_colors):
    """
    Create a palette of a specific number of colors.

    Args:
        palette_name (str): Name of standard palette of colors (should be called by matplotlib colormaps function).
        n_colors (int): Number of colors to integrate to palette

    Returns:
        list: List of colors from chosen standard palette, in RGB format
    """
    # Load standard color palette
    cmap = plt.colormaps[palette_name]
    if n_colors == 1:
        colors_list = [cmap(1)]
    else:
        # If more than 1 color => get a color palette with the number of needed colors
        colors_list = [cmap(i / (n_colors - 1)) for i in range(n_colors)]
    # Convert into RGB colors 
    return [(int(r * 250),int(g * 250),int(b * 250)) for r, g, b, _ in colors_list]

def one_color_per_track(path_to_circos, assembly_epi):
    """
    Color each track with a specific color.

    Args:
        path_to_circos (str): Path to circos results directory.
        assembly_epi (dict): Dictionary association between old track names and sample names ({assembly1: Epi-15, assembly2: Epi-412, ...})

    Returns:
        None
    """
    # Generate color range
    palette_name = 'Spectral'  # ColorBrewer palette
    n_colors = len(assembly_epi.keys()) 
    colors = generate_color_range(palette_name, n_colors)
    # In a python dictionary, associate each track to 1 color of the palette {Epi-name: color}
    assembly_color = {list(assembly_epi.values())[i]: str(colors[i]).strip("()").replace(" ", "") for i in range(len(assembly_epi.keys()))}
    # For each sample inside Circos graph, modify Circos config file
    for sample in assembly_epi.values():
        sample_conf = path_to_circos+"/data/"+sample+".conf"
        tmp_conf = path_to_circos+"/data/tmp.conf"
        sample_file = open(sample_conf, "r")
        tmp_file = open(tmp_conf, "w")
        # Replace old color for the sample, by the color attributed in previous step to the sample
        for line in sample_file:
            tmp_file.write(re.sub("color=[a-z]+", "color="+assembly_color[sample], line))
        sample_file.close()
        tmp_file.close()
        os.remove(sample_conf)
        os.rename(tmp_conf, sample_conf)

def redefine_tick_labels(path_to_circos):
    """
    Resize labels inside Circos graph.

    Args:
        path_to_circos (str): Path to circos results directory.

    Returns:
        None
    """
    ## Tick labels
    ticks_file = open(path_to_circos+"/data/ticks.conf", "r")
    tmp_file = open(path_to_circos+"/data/tmp.conf", "w")
    # Change size of every labels to 32pt
    for line in ticks_file:
        tmp_file.write(re.sub("label_size = [0-9]+p", "label_size = 32p", line))
    tmp_file.close()
    ticks_file.close()
    # Remove old configuration file
    os.remove(path_to_circos+"/data/ticks.conf")
    # Replace it by the new one
    os.rename(path_to_circos+"/data/tmp.conf", path_to_circos+"/data/ticks.conf")

def remove_track_blocks(input_file, output_file, track_id, track_pos, nb_block):
    """
    Remove a track from Circos graph.
    A track is described in configuration file by a block of nb_block lines.
    The function remove every line of the block associated to track inside configuration file.

    Args:
        input_file (str): Path to original circos main configuration file.
        output_file (str): Path to new circos configuration file that will contain the changes.
        track_id (str): Track number to remove (format="conf(track#_pos)")
        track_pos (int): Position (0-based) of line in block that contains track name to remove
        nb_block (int): Number of lines to check as a block 

    Returns:
        None
    """
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        buffer = []  # To store the last "nb_block" lines
        for line in infile:
            buffer.append(line)

            # Keep buffer size at most "nb_block" to scan file "nb_block" lines at a time
            if len(buffer) > nb_block:
                # Write the first line of the block in output and remove it from buffer
                outfile.write(buffer.pop(0))

            # Check the "track_pos" line in buffer for problem_id
            if len(buffer) == nb_block and track_id in buffer[track_pos]:
                buffer.clear()  # Remove the block if problem_id found

        # Write any remaining lines that didn't fit into a full block
        for remaining_line in buffer:
            outfile.write(remaining_line)

def remove_any_track(path_to_circos, tracks_corresp):
    """
    Remove ANY track from Circos graph.

    Args:
        path_to_circos (str): Path to circos results directory.
        tracks_corresp (str): Configuration of labels.

    Returns:
        None
    """
    # ANY track displays gene density all along the reference genome. 
    # As it does not add interesting information for us, it would be better to remove this track.
    
    # Get the track number associated to ANY
    any_track = re.findall("track[0-9]+=ANY", tracks_corresp)[0].replace("=ANY", "")
    # Parts that needed to be removed in circos configuration file
    tmp_file = path_to_circos+"/tmp.conf"
    circos_conf = path_to_circos+"/circos.conf"
    # Remove the ANY heatmap by changing config in a tmp file
    remove_track_blocks(circos_conf, tmp_file, "conf("+any_track+"_pos)", 4, 7)
    # Remove old config file
    os.remove(circos_conf)
    # Replace it with temporary file
    os.rename(tmp_file, circos_conf)
    # Remove ANY label
    remove_track_blocks(circos_conf, tmp_file, "track_idx = "+any_track, 1, 4)
    os.remove(circos_conf)
    os.rename(tmp_file, circos_conf)


def loading_sample_names(path_to_circos):
    """
    Associate sample labels to tracks drawn on graph.
    Change track values for sample names inside Circos configuration file.

    Args:
        path_to_circos (str): Path to circos results directory.

    Returns:
        list: string of modified labels and dictionary of association between old track labels and accurate sample names.
    """
    # Labels are stored in _data/labels.txt_ file, where track names are listed.
    # Assembly names and their associated track names are stored in _legend.txt_.
    label_file = open(path_to_circos+"/data/labels.txt", "r")
    old_label = label_file.readlines()
    label_file.close()
    assembly_name = open(path_to_circos+"/legend.txt", "r").readlines()
    tmp = [name for name in assembly_name if 'assembly' in name]
    tmp = [name.strip() for name in tmp]
    # transform into {'assembly1': "Epi-54"}
    assembly_track = {name.split(" ")[0]: name.split(" ")[2] for name in tmp}
    # modify assembly1,... by Epi-name
    for assembly in assembly_track.keys():
        old_label[0] = old_label[0].replace(assembly, assembly_track[assembly])
    # remove ANY track
    final_label = re.sub(",track[0-9]+=ANY", "", old_label[0])
    # create new labels.txt config file with Epi-name and without ANY track
    new_label = open(path_to_circos+"/data/labels.txt", "w")
    new_label.write(final_label)
    new_label.close()
    # Return modified labels config (with ANY) and association between old track labels and Epi-names
    return [old_label[0], assembly_track]

def change_working_dir(path_to_circos):
    """
    Change name of working directory inside Circos configuration file.

    Args:
        path_to_circos (str): Path to circos results directory.

    Returns:
        None
    """
    # Circos configuration file indicates specifically the name of directory. 
    # Since we change the directory to process files, we have to modify the variable in code.
    # Nextflow work/aa/ajfkl12... --> analyses/blaVIM-1/IncL-M/visual/...

    # Open circos configuration file
    circos_conf_file = open(path_to_circos+"/circos.conf", "r")
    # Create tmp file that will contain new config
    tmp_file = open(path_to_circos+"/tmp.conf", "w")
    # Change dir path
    for line in circos_conf_file:
        tmp_file.write(re.sub("dir = [^\n]+", "dir = "+path_to_circos, line))
    circos_conf_file.close()
    tmp_file.close()
    # Remove old config file
    os.remove(path_to_circos+"/circos.conf")
    # Replace it by tmp file
    os.rename(path_to_circos+"/tmp.conf", path_to_circos+"/circos.conf")

def main(project_dir):
    print("Processing some changes to Circos configuration files.")
    print("First step, updating working directory.")
    change_working_dir(project_dir)
    print("Second step, editing sample labels.")
    out = loading_sample_names(project_dir)
    assembly_tracks = out[0] # String with new labels in config
    assembly_epi = out[1] # Dictionary with correspondance between old labels and sample names
    print("Third step, remove ANY track from plot.")
    remove_any_track(project_dir, assembly_tracks)
    print("Step four, getting bigger tick labels.")
    redefine_tick_labels(project_dir)
    print("Step five, define one colour for each track.")
    one_color_per_track(project_dir, assembly_epi)
    print("Step six, adding scatterplots to represent mismatches.")
    scatter_for_mismatches(project_dir, assembly_tracks)
    print("Step seven, adjusting scale on kb.")
    scale_to_kb(project_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Modifying Circos configuration files')
    parser.add_argument("-d", dest="path_to_data", required=True, help="Path to the folder containing Circos files")
    parser.add_argument("-o", dest="output_dir", required=False, help="Path to output directory")
    args = parser.parse_args()

    if args.output_dir is None:
        args.output_dir = args.path_to_data

    path_to_data = args.path_to_data
    output_dir = args.output_dir

    main(path_to_data, output_dir)