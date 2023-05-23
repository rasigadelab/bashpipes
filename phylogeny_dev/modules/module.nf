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
    tuple val(replicon), path("phylogeny/$replicon/sequences/$sample/prokka/${sample}.gff"), emit : annotation_file
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

process pan_genome_panaroo {
  // Tool: panaroo. 
  // Pan-genome analysis of a batch of sample enables to get the core genome of a specific replicon.
  // Core-genome gathers genes that are present in every sample of the batch.

  label 'highCPU'
  storeDir params.result
  debug false
  tag "Panaroo on $replicon"  

  when:
    params.pan_genome_panaroo.todo == 1
  
  input:
    tuple val(replicon), path(annotation_files)

  output:
    tuple val(replicon), path("phylogeny/$replicon/panaroo/core_gene_alignment_filtered.aln"), emit : core_genome_aln
    tuple val(replicon), path("phylogeny/$replicon/panaroo/pan_genome_reference.fa"), path("phylogeny/$replicon/panaroo/core_alignment_filtered_header.embl"), emit : core_genome_ref
    path("phylogeny/$replicon/panaroo/summary_statistics.txt")
    path("phylogeny/$replicon/panaroo/panaroo.log")
    path("phylogeny/$replicon/panaroo/core_genome_reference.fa")

  script:
  """
  OUT_DIR=phylogeny/$replicon/panaroo
  mkdir -p -m 777 \${OUT_DIR}
  
  #Step1- Running Panaroo
  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate panaroo
  panaroo -i $annotation_files -o \${OUT_DIR} --clean-mode ${params.pan_genome_panaroo["clean_mode"]} \
    -a ${params.pan_genome_panaroo["gene_alignment"]} --core_threshold ${params.pan_genome_panaroo["core_threshold"]} \
    -t $task.cpus &> \${OUT_DIR}/panaroo.log
  conda deactivate

  #Step2- Creating core genome reference FASTA file
  python3 ${params.nfpath}/modules/core_genome.py -d \${OUT_DIR} 
  """

  stub:
  """
  OUT_DIR=phylogeny/$replicon/panaroo
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/core_gene_alignment_filtered.aln
  touch \${OUT_DIR}/pan_genome_reference.fa
  touch \${OUT_DIR}/core_alignment_filtered_header.embl
  touch \${OUT_DIR}/summary_statistics.txt
  touch \${OUT_DIR}/panaroo.log
  touch \${OUT_DIR}/core_genome_reference.fa
  """  
}

process core_tree_iqtree {
  // Tool: iqtree. 
  // Tree construction of a specific replicon based on its core genome.

  label 'highCPU'
  storeDir params.result
  debug false
  tag "Iqtree on $replicon"  

  when:
    params.core_tree_iqtree.todo == 1
  
  input:
    tuple val(replicon), path(core_genome_aln)

  output:
    tuple val(replicon), path("phylogeny/$replicon/iqtree_after_panaroo/${replicon}.treefile")
    path("phylogeny/$replicon/iqtree_after_panaroo/iqtree.err")
    path("phylogeny/$replicon/iqtree_after_panaroo/iqtree.log")

  script:
  """
  OUT_DIR=phylogeny/$replicon/iqtree_after_panaroo
  mkdir -p -m 777 \${OUT_DIR}
  
  iqtree -s $core_genome_aln -m ${params.core_tree_iqtree["model"]} --prefix $replicon -T $task.cpus 1> \${OUT_DIR}/iqtree.log 2> \${OUT_DIR}/iqtree.err
  """

  stub:
  """
  OUT_DIR=phylogeny/$replicon/iqtree_after_panaroo
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/${replicon}.treefile
  touch \${OUT_DIR}/iqtree.err
  touch \${OUT_DIR}/iqtree.log
  """  
}

process core_snps_snippy {
  // Tool: snippy. 
  // Illumina FASTQ mapping on core_genome reference FASTA file.
  // Then, SNP calling on mapping results. 

  label 'highCPU'
  storeDir params.result
  debug false
  tag "Snippy on $replicon"  

  when:
    params.core_snps_snippy.todo == 1
  
  input:
    tuple val(replicon), path(pan_genome_ref), path(core_genes)

  output:
    tuple val(replicon), path("phylogeny/$replicon/snippy")

  script:
  """
  OUT_DIR=phylogeny/$replicon/snippy
  mkdir -p -m 777 \${OUT_DIR}

  #First step: getting core genome reference FASTA
  python3 core_genome.py -d phylogeny/$replicon/panaroo 

  #Second step: creating input-tab for snippy-multi
  echo $samples

  
  """

  stub:
  """
  OUT_DIR=phylogeny/$replicon/snippy
  mkdir -p -m 777 \${OUT_DIR}
  echo $samples
  """ 

}

  








