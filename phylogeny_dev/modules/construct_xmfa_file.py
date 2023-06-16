from Bio import SeqIO
import argparse
import os
import pandas

# GOAL : build .xmfa file for aln input for ClonalFrameML

def main(prokka_dir, panaroo_dir, output_dir):
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

    #Step2- Open Panaroo metadata on genes
    panaroo_metadata = os.path.join(panaroo_dir, "gene_presence_absence.csv")
    metadata_df = pandas.read_csv(panaroo_metadata, sep = ',', index_col = "Gene")
    panaroo_gene_data = os.path.join(panaroo_dir, "gene_data.csv")
    gene_data_df = pandas.read_csv(panaroo_gene_data, sep = ',', index_col = "clustering_id")
    #Step2- Get core genes alignment file and write it to XMFA file
    for core_gene in core_labels:
        aln_file = os.path.join(panaroo_dir,"aligned_gene_sequences",core_gene+".aln.fas")
        for seq_record in SeqIO.parse(aln_file, "fasta"):
            # Work on having right IDs [>sample_name:start-end +gene_name]
            seq_id = seq_record.id
            # Correct if panaroo renamed gene by '_R_gene'
            seq_id = seq_id.replace('_R_', '')
            # Get gene start and end for sample
            sample = seq_id.split(";")[0]
            #Strat2 - step1 - Check the gene_id for the sample in panaroo_metadata
            gff_id = metadata_df.loc[core_gene, sample]
            # If it's a gene present in GFF file
            if "refound" not in gff_id:
                #Strat2 - step2 - Go to gff file and grep the line containing gff_id
                gff_path = os.path.join(prokka_dir, sample, "prokka", sample+'.gff')
                gff_file=open(gff_path, "r")
                seq_annotation = ""
                for line in gff_file:
                    if gff_id in line:
                        seq_annotation = line
                        
                gff_file.close()
                #Strat2 - step3 - grep start and end of gene sequence
                seq_start = seq_annotation.split('\t')[3]
                seq_end = seq_annotation.split('\t')[4]
            # If it's a gene refound by Panaroo
            else:
                #Strat2 - step2 - Go to 'gene_data.csv' file
                data = gene_data_df.loc[gff_id, "description"]
                data = data.split(";")
                for item in data:
                    if "location" in item:
                        data = item
                        break
                seq_start = data.split(":")[1].split("-")[0]
                seq_end = data.split(":")[1].split("-")[1]

            #Strat2 - final step :)
            new_seq_id = sample+":"+seq_start+"-"+seq_end+" +"+core_gene
            seq = seq_record.seq
            out_xmfa.write(">"+new_seq_id+"\n")
            out_xmfa.write(str(seq)+"\n")
        # Adding an "=" between each FASTA block, that's .XMFA formatting
        out_xmfa.write("=\n")


        # with open(aln_file) as f:
        #     lines = f.readlines()
        #     out_xmfa.writelines(lines)
        #     out_xmfa.write("=\n")
            
    core_genes.close()
    out_xmfa.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Getting core reference genome FASTA')
    parser.add_argument("-d", dest="path_to_prokka_dir", required=True, help="Path to the folder containing prokka results")
    parser.add_argument("-p", dest="path_to_panaroo_dir", required=True, help="Path to the folder containing panaroo results")
    parser.add_argument("-o", dest="output_dir", required=False, help="Path to output directory")
    args = parser.parse_args()

    if args.output_dir is None:
        args.output_dir = args.path_to_panaroo_dir

    path_to_prokka_dir = args.path_to_prokka_dir
    path_to_panaroo_dir = args.path_to_panaroo_dir
    output_dir = args.output_dir

    main(path_to_prokka_dir, path_to_panaroo_dir, output_dir)

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