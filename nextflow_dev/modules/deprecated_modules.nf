process create_files_location {
  //DEPRECATED, launched before nextflow pipeline
  //Worklflow command to build a channel from CSV : 
  //file_locations = Channel.fromPath(params.result+"/Files_location.tsv", checkIfExists:true).splitCsv(sep:'\t', header: true)

  tag "Creating Files_location.tsv"
  publishDir "$sample_path", mode: 'copy'
  debug true

  input:
    val sample_path

  output:
    path "Files_location.tsv"

  """
  python3 /home/afischer/script_to_try/files_location.py -d $sample_path -o .
  """
}

process merger_preparation {
  //DEPRECATED, has been replaced by process make_sample_dir()

  input:
    path x
  
  output:
    stdout

  """
  sudo Rscript /mnt/c/Users/admin/Documents/Projets/4.Pipeline_improvements/5.Prepare_reads/230104_nextflow_prepare_reads.R -d $x
  """
}