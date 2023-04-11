process make_sample_dir{

  // Creating genomes folder, containing one folder per sample.
  // In each sample folder, one can find two or three files, named : sample_R1.fastq.gz, sample_R2.fastq.gz, sample_ONT.fastq.gz

  label 'highCPU'
  storeDir params.result
  debug false
  tag "Merging file $full_path.simpleName"

  input:
    tuple path(full_path), val(sample), val(type)

  output:
    tuple val(sample), path("genomes/$sample/${sample}_${type}.fastq.gz")
  
  """
  mkdir -p -m 777 genomes/$sample
  cat $full_path >> genomes/$sample/${sample}_${type}.fastq.gz
  """
}