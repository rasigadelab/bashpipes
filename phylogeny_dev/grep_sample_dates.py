import pandas
import os
import argparse

#GOAL : Creation of a TXT file containing Epi-1       \t 2023-02-04
#                                         Sample_name \t Sampling_date

## FUNCTIONS ##

def main(project_dir, output_dir):
    #Step1- Read Metadataeq_rasigade.xlsx
    excel_file = os.path.join("metadata.xlsx")
    #Get column 1 = SAMPLE_ID and column 4 = DATE_PRELEVEMENT
    df = pandas.read_excel(excel_file, usecols = ['SAMPLE_ID', 'DATE_PRELEVEMENT'])
    #Read replicons.tsv
    sample_file = os.path.join(project_dir,"replicons.tsv")
    df2 = pandas.read_csv(sample_file, sep = '\t', usecols = ['Sample', 'replicon'])
    df2 = df2.reset_index()  # make sure indexes pair with number of rows
    #Create a dictionary {replicon: samples_list}
    replicon_dict = {}
    for index, row in df2.iterrows():
        rep = row['replicon']
        sample = row['Sample']
        if rep not in replicon_dict.keys():
            replicon_dict[rep] = [sample]
        else:
            replicon_dict[rep].append(sample)
    #Write down to output file
    for repl in replicon_dict.keys():
        out_file = open(os.path.join(output_dir,repl+"_sampling_list.csv"), "w")
        out_file.write("strain,date\n")
        for sample in replicon_dict[repl]:
            sample_date = df.loc[df['SAMPLE_ID'] == sample, ['DATE_PRELEVEMENT']].values[0][0]
            out_file.write(sample+","+sample_date+"\n")
        out_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Creation of sampling list for each replicon')
    parser.add_argument("-d", dest="path_to_data", required=True, help="Path to the folder containing the genomes repository")
    parser.add_argument("-o", dest="output_dir", required=False, help="Path to output directory")
    args = parser.parse_args()

    if args.output_dir is None:
        args.output_dir = args.path_to_data

    path_to_data = args.path_to_data
    output_dir = args.output_dir

    main(path_to_data, output_dir)

# THE END



