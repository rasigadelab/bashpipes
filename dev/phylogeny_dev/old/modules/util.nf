process make_sample_dir{

  // Creating genomes folder, containing one folder per sample.
  // In each sample folder, one can find two or three files, named : sample_R1.fastq.gz, sample_R2.fastq.gz, sample_ONT.fastq.gz

  label 'lowCPU'
  storeDir params.result
  debug false
  tag "Merging files of sample $sample of type $type"

  input:
    tuple val(sample), val(type), path(list_of_filepaths)

  output:
    tuple val(sample), path("genomes/$sample/${sample}_${type}.fastq.gz")
  
  """
  mkdir -p -m 777 genomes/$sample
  cat $list_of_filepaths >> genomes/$sample/${sample}_${type}.fastq.gz
  """
}