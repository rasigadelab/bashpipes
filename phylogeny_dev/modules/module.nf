process annotate_prokka {
  // Tool: prokka. 
  // Genome annotation consists in locating genes on the genome and giving their function 

  label 'lowCPU'
  storeDir params.result
  debug false
  tag "Prokka on $sample"  

  when:
    params.annotate_prokka.todo == 1

  input:
    tuple val(replicon), val(sample), path(fasta_file)

  output:
    tuple val(replicon), val(sample), path("phylogeny/$replicon/sequences/$sample/prokka/${sample}.gff"), emit : final_assembly
    path("phylogeny/$replicon/sequences/$sample/prokka/*.gff")
    path("phylogeny/$replicon/sequences/$sample/prokka/*.log")

  script:
  """
  OUT_DIR=phylogeny/$replicon/sequences/$sample/prokka
  mkdir -p -m 777 \${OUT_DIR}
  
  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate prokka
  prokka $fasta_file --force ${params.annotate_prokka["genes"]} –-cpus $task.cpus --outdir \${OUT_DIR} ${params.annotate_prokka["mode"]} \
    --prefix $sample &> \${OUT_DIR}/prokka.log
  conda deactivate
  """

  stub:
  """
  OUT_DIR=phylogeny/$replicon/sequences/$sample/prokka
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/${sample}.gff
  touch \${OUT_DIR}/prokka.log
  """  
}









