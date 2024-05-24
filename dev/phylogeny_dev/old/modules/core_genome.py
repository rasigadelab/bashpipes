from Bio import SeqIO
import argparse
import os

def main(panaroo_dir, output_dir):
    core_genes = open(os.path.join(panaroo_dir,"core_alignment_filtered_header.embl"), "r")
    pan_genome = os.path.join(panaroo_dir, "pan_genome_reference.fa")
    out_core_genome = open(os.path.join(output_dir,"core_genome_reference.fa"), "w")
    #Step1- Get all core genes labels in a list
    core_labels = []
    for line in core_genes:
        line = line.split()
        if line[0] == "FT" and len(line) == 2:
            label = line[1].split("=")[1].replace('.aln','')
            if label not in core_labels:
                core_labels.append(label)

    #Step2- Get core genes ref sequence and write it to output file
    for seq_record in SeqIO.parse(pan_genome, "fasta"):
        gene_id = seq_record.id 
        if gene_id in core_labels:
            gene_seq = seq_record.seq
            out_core_genome.write(">"+gene_id+"\n")
            out_core_genome.write(str(gene_seq)+"\n")

    core_genes.close()
    out_core_genome.close()

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