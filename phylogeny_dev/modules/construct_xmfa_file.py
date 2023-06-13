from Bio import SeqIO
import argparse
import os

# GOAL : build .xmfa file for aln input for ClonalFrameML

def main(panaroo_dir, output_dir):
    core_genes = open(os.path.join(panaroo_dir,"core_alignment_filtered_header.embl"), "r")
    out_xmfa = open(os.path.join(output_dir,"core_genes_aln.xmfa"), "w")
    #Step1- Get all core genes labels in a list
    core_labels = []
    for line in core_genes:
        line = line.split()
        if line[0] == "FT" and len(line) == 2:            
            label = line[1].split("=")[1].replace('.aln','')
            if label not in core_labels:
                core_labels.append(label)

    #Step2- Get core genes alignment file and write it to XMFA file
    for core_gene in core_labels:
        aln_file = os.path.join(panaroo_dir,"aligned_gene_sequences",core_gene+".aln.fas")
        for seq_record in SeqIO.parse(aln_file, "fasta"):
            seq_id = seq_record.id
            new_seq_id = seq_id.split(";")[0]+" +"+core_gene
            seq = seq_record.seq
            out_xmfa.write(">"+new_seq_id+"\n")
            out_xmfa.write(str(seq)+"\n")
        out_xmfa.write("=\n")


        # with open(aln_file) as f:
        #     lines = f.readlines()
        #     out_xmfa.writelines(lines)
        #     out_xmfa.write("=\n")
            
    core_genes.close()
    out_xmfa.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Getting core reference genome FASTA')
    parser.add_argument("-d", dest="path_to_panaroo", required=True, help="Path to the folder containing Panaroo outputs")
    parser.add_argument("-o", dest="output_dir", required=False, help="Path to output directory")
    args = parser.parse_args()

    if args.output_dir is None:
        args.output_dir = args.path_to_panaroo

    path_to_panaroo = args.path_to_panaroo
    output_dir = args.output_dir

    main(path_to_panaroo, output_dir)

# THE END

# What does it look like ? Well what is painful is that each gene has to be separated.
#[Core gene 1 block]
#>Sample1
#ATCGCTV-FTFYH
#>Sample2
#jnlUHD
#....
#>SampleN
#kfmeh
#=
#[Core gene 2 block]
#>Sample1
#ATCGCTV-FTFYH
#>Sample2
#jnlUHD
#....
#>SampleN
#kfmeh
#Does eventually Panaroo have an option already implemented ? :D