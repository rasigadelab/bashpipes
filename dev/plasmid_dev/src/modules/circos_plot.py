#!/usr/bin/python3
# coding: utf-8

import re
import os
import argparse
import matplotlib.pyplot as plt # type: ignore
from matplotlib import colormaps # type: ignore
import numpy as np # type: ignore




def scatter_for_mismatches(path_to_circos, track_corresp):
    # Scatterplot on histograms presenting mismatches
    # First change fill_color of every histogram to black
    circos_conf = open(path_to_circos+"/circos.conf", "r")
    tmp_conf = open(path_to_circos+"/tmp.conf", "w")
    for line in circos_conf:
        tmp_conf.write(line.replace("fill_color = vlyellow", "fill_color = black"))
    circos_conf.close()
    tmp_conf.close()
    os.remove(path_to_circos+"/circos.conf")
    os.rename(path_to_circos+"/tmp.conf", path_to_circos+"/circos.conf")
    # Then, add a scatterplot on histograms so points are clearly indicated. 
    # This must be added between "plots "/plots".
    # Get Epi-name and track combination
    tmp = track_corresp.split("\t")[-1].split(",")
    epi_track = {name.split("=")[1]: name.split("=")[0] for name in tmp}
    epi_track.pop("ANY")
    all_scatters = ""
    for sample in epi_track.keys():
        scatter_to_add = "<plot>\ntype = scatter\nthickness = 1\nfile = data/" + sample
        scatter_to_add += ".mismatches.txt\nr0 = eval(sprintf('%.3fr', conf(" + epi_track[sample] + "_pos) - "
        scatter_to_add += " conf(track_width)))\nr1 = eval(sprintf('%.3fr', conf(" + epi_track[sample] + "_pos)))\n"
        scatter_to_add += "stroke_color = black\nstroke_thickness = 2\nglyph_size = 12\n</plot>"
        all_scatters += scatter_to_add + "\n"
    circos_conf = open(path_to_circos+"/circos.conf", "r")
    tmp_conf = open(path_to_circos+"/tmp.conf", "w")
    for line in circos_conf:
        if "</plots>" in line:
            tmp_conf.write(all_scatters)
        tmp_conf.write(line)
    circos_conf.close()
    tmp_conf.close()
    os.remove(path_to_circos+"/circos.conf")
    os.rename(path_to_circos+"/tmp.conf", path_to_circos+"/circos.conf")


def generate_color_range(palette_name, n_colors):
    # Creation of a color palette associated to each sample
    cmap = plt.colormaps[palette_name]
    if n_colors == 1:
        colors_list = [cmap(1)]
    else:
        colors_list = [cmap(i / (n_colors - 1)) for i in range(n_colors)]
    return [(int(r * 250),int(g * 250),int(b * 250)) for r, g, b, _ in colors_list]

def one_color_per_track(path_to_circos, assembly_epi):
    # Generate color range
    palette_name = 'Spectral'  # ColorBrewer palette
    n_colors = len(assembly_epi.keys())
    colors = generate_color_range(palette_name, n_colors)
    # In a python dictionary
    assembly_color = {list(assembly_epi.values())[i]: str(colors[i]).strip("()").replace(" ", "") for i in range(len(assembly_epi.keys()))}
    # Replace the element color=corresponding color
    for sample in assembly_epi.values():
        sample_conf = path_to_circos+"/data/"+sample+".conf"
        tmp_conf = path_to_circos+"/data/tmp.conf"
        sample_file = open(sample_conf, "r")
        tmp_file = open(tmp_conf, "w")
        for line in sample_file:
            tmp_file.write(re.sub("color=[a-z]+", "color="+assembly_color[sample], line))
        sample_file.close()
        tmp_file.close()
        os.remove(sample_conf)
        os.rename(tmp_conf, sample_conf)

def redefine_tick_labels(path_to_circos):
    ## Tick labels
    ticks_file = open(path_to_circos+"/data/ticks.conf", "r")
    tmp_file = open(path_to_circos+"/data/tmp.conf", "w")

    for line in ticks_file:
        tmp_file.write(re.sub("label_size = [0-9]+p", "label_size = 32p", line))

    tmp_file.close()
    ticks_file.close()
    os.remove(path_to_circos+"/data/ticks.conf")
    os.rename(path_to_circos+"/data/tmp.conf", path_to_circos+"/data/ticks.conf")

def remove_track_blocks(input_file, output_file, track_id, track_pos, nb_block):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        buffer = []  # To store the last four lines
        for line in infile:
            buffer.append(line)

            # Keep buffer size at most 4
            if len(buffer) > nb_block:
                outfile.write(buffer.pop(0))

            # Check the second line in buffer for problem_id
            if len(buffer) == nb_block and track_id in buffer[track_pos]:
                buffer.clear()  # Remove the block if problem_id found

        # Write any remaining lines that didn't fit into a full block
        for remaining_line in buffer:
            outfile.write(remaining_line)

def remove_any_track(path_to_circos, tracks_corresp):
    # ANY track displays gene density all along the reference genome. 
    # As it does not add interesting information for us, it would be better to remove this track.
    any_track = re.findall("track[0-9]+=ANY", tracks_corresp)[0].replace("=ANY", "")
    # Parts that needed to be removed in circos configuration file
    tmp_file = path_to_circos+"/tmp.conf"
    circos_conf = path_to_circos+"/circos.conf"
    remove_track_blocks(circos_conf, tmp_file, "conf("+any_track+"_pos)", 4, 7)
    os.remove(circos_conf)
    os.rename(tmp_file, circos_conf)

    remove_track_blocks(circos_conf, tmp_file, "track_idx = "+any_track, 1, 4)
    os.remove(circos_conf)
    os.rename(tmp_file, circos_conf)


def loading_sample_names(path_to_circos):
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
    # change labels.txt
    for assembly in assembly_track.keys():
        old_label[0] = old_label[0].replace(assembly, assembly_track[assembly])
    # supprimer track ANY
    final_label = re.sub(",track[0-9]+=ANY", "", old_label[0])
    new_label = open(path_to_circos+"/data/labels.txt", "w")
    new_label.write(final_label)
    new_label.close()
    return [old_label[0], assembly_track]

def change_working_dir(path_to_circos):
    # Circos configuration file indicates specifically the name of directory. 
    # Since we change the directory to process files, we have to modify the variable in code.
    circos_conf_file = open(path_to_circos+"/circos.conf", "r")
    tmp_file = open(path_to_circos+"/tmp.conf", "w")
    for line in circos_conf_file:
        tmp_file.write(re.sub("dir = [^\n]+", "dir = "+path_to_circos, line))
    circos_conf_file.close()
    tmp_file.close()
    os.remove(path_to_circos+"/circos.conf")
    os.rename(path_to_circos+"/tmp.conf", path_to_circos+"/circos.conf")

def main(project_dir, output_dir):
    print("Processing some changes to Circos configuration files.")
    print("First step, updating working directory.")
    change_working_dir(project_dir)
    print("Second step, editing sample labels.")
    out = loading_sample_names(project_dir)
    assembly_tracks = out[0]
    assembly_epi = out[1]
    print("Third step, remove ANY track from plot.")
    remove_any_track(project_dir, assembly_tracks)
    print("Step four, getting bigger tick labels.")
    redefine_tick_labels(project_dir)
    print("Step five, define one colour for each track.")
    one_color_per_track(project_dir, assembly_epi)
    print("Step six, adding scatterplots to represent mismatches.")
    scatter_for_mismatches(project_dir, assembly_tracks)

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