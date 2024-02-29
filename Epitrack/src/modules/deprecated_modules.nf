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

process annotate_prokka {
  // Tool: prokka. 
  // Genome annotation consists in locating genes on the genome and giving their function 

  label 'prokka'
  storeDir params.result
  debug false
  tag "Prokka on $sample"  

  when:
    params.annotate_prokka.todo == 1

  input:
    tuple val(sample), path(final_assembly), path(taxonomy_file)

  output:
    tuple val(sample), path("genomes/$sample/$final_assembly"), emit : final_assembly
    path("genomes/$sample/prokka/${sample}.err")
    path("genomes/$sample/prokka/${sample}.gff")
    path("genomes/$sample/prokka/${sample}.log")

  script:
  """
  OUT_DIR=genomes/$sample/prokka
  mkdir -p -m 777 \${OUT_DIR}

  # Extract genus and species names if available
  GENUS=\$(cut -d',' -f8 $taxonomy_file | tail -n 1)
  SPECIES=\$(cut -d',' -f9 $taxonomy_file | tail -n 1)
  
  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate prokka
  prokka $final_assembly --force ${params.annotate_prokka["genes"]} –-cpus $task.cpus --outdir \${OUT_DIR} ${params.annotate_prokka["mode"]} \
      --prefix $sample –-usegenus –-genus \$GENUS –-species \$SPECIES &> \${OUT_DIR}/prokka.log
  conda deactivate
  """
}