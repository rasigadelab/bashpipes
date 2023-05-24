import os
import argparse
import ast 

# Create input.tab, list of samples to give to snippy.

def main(work_dir, samples_list, out_dir):
    genomes_dir = os.path.join(work_dir, "genomes")
    #samples_list = ast.literal_eval(samples_list)
    #Step 1- Check if there is already an output directory
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    #Step 2- Write input.tab file in output dir [Sample \t R1_path \t R2_path]
    out_file = open(os.path.join(out_dir, "input.tab"), "w")
    for sample in samples_list:
        R1 = os.path.join(genomes_dir, sample, sample+"_R1.fastq.gz")
        R2 = os.path.join(genomes_dir, sample, sample+"_R2.fastq.gz")
        out_file.write(sample+"\t"+R1+"\t"+R2+"\n")
    out_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Creating Snippy-multi "input.tab" file''')
    parser.add_argument("-d", dest="path_to_workdir", required=True, help="Path to the folder containing genomes and phylogeny folders")
    parser.add_argument("-l", dest="list_of_samples", required=True, nargs='+', help="List of a batch of samples owning the same replicon")
    parser.add_argument("-o", dest="path_to_outdir", required=True, help="Path to the folder that will contain output files")
    args = parser.parse_args()

    path_to_workdir = args.path_to_workdir
    list_of_samples = args.list_of_samples
    path_to_outdir = args.path_to_outdir

    main(path_to_workdir, list_of_samples, path_to_outdir)

# THE END
