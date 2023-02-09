process make_sample_dir{
  tag "Merging file $full_path.simpleName"
  publishDir "$sample_path", mode: 'copy'
  debug true

  input:
    tuple path(full_path), val(sample), val(type)
    val sample_path

  output:
    stdout
  
  """
  mkdir -p -m 777 $sample_path/genomes/$sample
  cat $full_path >> $sample_path/genomes/$sample/${sample}_${type}.fastq.gz
  """
}