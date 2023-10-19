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
    tuple val(replicon), val(sample), path(fasta_file)

  output:
    tuple val(replicon), val(sample), path("phylogeny/$replicon/sequences/$sample/prokka/${sample}.gff"), emit : annotation_file
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

  label 'panaroo'
  storeDir params.result
  debug true
  tag "Panaroo on $replicon"  

  when:
    params.pan_genome_panaroo.todo == 1
  
  input:
    tuple val(replicon), val(samples), path(annotation_files)

  output:
    tuple val(replicon), path("phylogeny/$replicon/panaroo/core_gene_alignment_filtered.aln"), emit : core_genome_aln
    tuple val(replicon), val(samples), path("phylogeny/$replicon/panaroo/core_genome_reference.fa"), emit : core_genome_ref
    //tuple val(replicon), val(samples), path("phylogeny/$replicon/panaroo/core_genome_reference.fa"), path("phylogeny/$replicon/panaroo/core_genes_aln.xmfa"), emit : core_genome_ref
    path("phylogeny/$replicon/panaroo/summary_statistics.txt")
    path("phylogeny/$replicon/panaroo/panaroo.log")

  script:
  """
  OUT_DIR=phylogeny/$replicon/panaroo
  mkdir -p -m 777 \${OUT_DIR}
  
  #Step1- Running Panaroo
  panaroo -i $annotation_files -o \${OUT_DIR} --clean-mode ${params.pan_genome_panaroo["clean_mode"]} \
    -a ${params.pan_genome_panaroo["gene_alignment"]} --core_threshold ${params.pan_genome_panaroo["core_threshold"]} \
    -t $task.cpus &> \${OUT_DIR}/panaroo.log

  #Step2- Creating core genome reference FASTA file
  python3 ${params.nfpath}/modules/core_genome.py -d \${OUT_DIR} 
  
  #Step3- Creating .xmfa core genome alignment file
  #python3 ${params.nfpath}/modules/construct_xmfa_file.py -d ${params.result}/phylogeny/$replicon/sequences -p \${OUT_DIR}

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
  touch \${OUT_DIR}/core_genes_aln.xmfa
  """  
}

process core_tree_iqtree {
  // Tool: iqtree. 
  // Tree construction of a specific replicon based on its core genome.

  label 'lowCPU'
  storeDir params.result
  debug false
  tag "Iqtree on $replicon"  

  when:
    params.core_tree_iqtree.todo == 1
  
  input:
    tuple val(replicon), path(core_genome_aln)
    val(out_prefix)

  output:
    tuple val(replicon), path("phylogeny/$replicon/$out_prefix/${replicon}.treefile")
    path("phylogeny/$replicon/$out_prefix/*")
    // path("phylogeny/$replicon/$out_prefix/${replicon}.iqtree")
    // path("phylogeny/$replicon/$out_prefix/iqtree.err")
    // path("phylogeny/$replicon/$out_prefix/iqtree.log")
    // path("phylogeny/$replicon/$out_prefix/${replicon}.alninfo")

  script:
  """
  OUT_DIR=phylogeny/$replicon/$out_prefix
  mkdir -p -m 777 \${OUT_DIR}

  
  
  iqtree -s $core_genome_aln -m ${params.core_tree_iqtree["model"]} --prefix \${OUT_DIR}/$replicon -T $task.cpus -alninfo 1> \${OUT_DIR}/iqtree.log 2> \${OUT_DIR}/iqtree.err
  """

  stub:
  """
  OUT_DIR=phylogeny/$replicon/$out_prefix
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/${replicon}.treefile
  touch \${OUT_DIR}/iqtree.err
  touch \${OUT_DIR}/iqtree.log
  """  
}

process create_input_tab {
  // Tool: homemade script. 
  // Creating tsv file 'input.tab' for snippy usage (Sample, Fasta_R1, Fasta_R2)

  label 'lowCPU'
  storeDir params.result
  debug false
  tag "Input tab for $replicon"  

  when:
    params.core_snps_snippy.todo == 1
  
  input:
    tuple val(replicon), val(samples), path(core_genome_fasta)//, path(core_genes_aln)

  output:
    tuple val(replicon), path("$core_genome_fasta"), path("phylogeny/$replicon/snippy/input.tab"), emit : input_tab
    //tuple val(replicon), path("$core_genome_fasta"), path("phylogeny/$replicon/snippy/input.tab"), path("$core_genes_aln"), emit : input_tab

  script:
  """
  OUT_DIR=phylogeny/$replicon/snippy
  mkdir -p -m 777 \${OUT_DIR}
  
  #samples = [1, 2, 3] But I want : samples = 1 2 3
  SAMPLES_LIST=\$(echo "$samples" | tr -d "[]" | tr -d ",")

  python3 ${params.nfpath}/modules/snippy_multi_list.py -d $params.result -l \$SAMPLES_LIST -o \${OUT_DIR}
  """

  stub:
  """
  OUT_DIR=phylogeny/$replicon/snippy
  mkdir -p -m 777 \${OUT_DIR}
  #samples = [1, 2, 3] But I want : samples = 1 2 3
  SAMPLES_LIST=\$(echo "$samples" | tr -d "[]" | tr -d ",")
  python3 /home/afischer/11.Phylogeny_test/snippy_multi_list.py -d $params.result -l \$SAMPLES_LIST -o \${OUT_DIR}
  """

}

process core_snps_snippy {
  // Tool: snippy. 
  // Illumina FASTQ mapping on core_genome reference FASTA file.
  // Then, SNP calling on mapping results. 

  label 'snippy'
  storeDir params.result
  debug false
  tag "Snippy on $replicon" 

  when:
    params.core_snps_snippy.todo == 1
  
  input:
    tuple val(replicon), path(core_genome_ref), path(input_tab)//, path(core_genes_aln)

  output:
    tuple val(replicon), path("phylogeny/$replicon/snippy/core_without_ref.aln"), emit : core_snps_aln
    //tuple val(replicon), path("phylogeny/$replicon/snippy/core_without_ref.aln"), path("$core_genes_aln"), emit : core_snps_aln
    path("phylogeny/$replicon/snippy/*/snps.aligned.fa")
    path("phylogeny/$replicon/snippy/*/snps.log")
    path("phylogeny/$replicon/snippy/*/snps.subs.vcf")
    path("phylogeny/$replicon/snippy/core.tab")
    path("phylogeny/$replicon/snippy/core.txt")
    path("phylogeny/$replicon/snippy/core.vcf")
    path("phylogeny/$replicon/snippy/snippy.log")
    path("phylogeny/$replicon/snippy/snippy.err")

  script:
  """
  OUT_DIR=phylogeny/$replicon/snippy
  mkdir -p -m 777 \${OUT_DIR}

  snippy-multi $input_tab --ref $core_genome_ref --cpus $task.cpus --force --mincov 30 > \${OUT_DIR}/snippy_commands.sh
  sed -i -e "s|snippy-core --ref '|snippy-core --prefix \${OUT_DIR}/core --ref '|g" \${OUT_DIR}/snippy_commands.sh
  sh \${OUT_DIR}/snippy_commands.sh 1> \${OUT_DIR}/snippy.log 2> \${OUT_DIR}/snippy.err
  # Remove 2 last lines because it's the reference sequence and we don't want it in phylogeny
  head -n -2 \${OUT_DIR}/core.full.aln > \${OUT_DIR}/core_without_ref.aln
  """

  stub:
  """
  OUT_DIR=phylogeny/$replicon/snippy
  mkdir -p -m 777 \${OUT_DIR}/sample
  touch "\${OUT_DIR}/core.full.aln"
  touch "\${OUT_DIR}/core.aln"
  touch "\${OUT_DIR}/core_without_ref.aln"
  touch "\${OUT_DIR}/sample/snps.aligned.fa"
  touch "\${OUT_DIR}/sample/snps.log"
  touch "\${OUT_DIR}/sample/snps.subs.vcf"
  touch "\${OUT_DIR}/core.tab"
  touch "\${OUT_DIR}/core.txt"
  touch "\${OUT_DIR}/core.vcf"
  touch "\${OUT_DIR}/snippy.log"
  touch "\${OUT_DIR}/snippy.err"
  """ 

}

process snps_tree_iqtree {
  // Tool: iqtree. 
  // Tree construction of a specific replicon based on its core snps.

  label 'lowCPU'
  storeDir params.result
  debug false
  tag "Iqtree on $replicon"  

  when:
    params.snps_tree_iqtree.todo == 1
  
  input:
    tuple val(replicon), path(core_snps_aln)//, path(core_genes_aln)
    val(out_prefix)

  output:
    tuple val(replicon), path("phylogeny/$replicon/$out_prefix/${replicon}.treefile"), emit : treefiles
    path("phylogeny/$replicon/$out_prefix/*")
    // path("phylogeny/$replicon/$out_prefix/${replicon}.iqtree")
    // path("phylogeny/$replicon/$out_prefix/iqtree.err")
    // path("phylogeny/$replicon/$out_prefix/iqtree.log")
    // path("phylogeny/$replicon/$out_prefix/${replicon}.alninfo")
    // path("phylogeny/$replicon/$out_prefix/${replicon}.timetree.lsd")
    // path("phylogeny/$replicon/$out_prefix/${replicon}.timetree.nwk")
    // path("phylogeny/$replicon/$out_prefix/${replicon}.timetree.nex")

  script:
  """
  OUT_DIR=phylogeny/$replicon/$out_prefix
  mkdir -p -m 777 \${OUT_DIR}
  # TO DELETE AFTER
  sed -n '/>Reference/q;p' $core_snps_aln > \${OUT_DIR}/snps.aln
  iqtree -s \${OUT_DIR}/snps.aln -m ${params.snps_tree_iqtree["model"]} --prefix \${OUT_DIR}/$replicon -T $task.cpus -alninfo --date $params.result/${replicon}_sampling_list.txt 1> \${OUT_DIR}/iqtree.log 2> \${OUT_DIR}/iqtree.err
  
  """

  stub:
  """
  OUT_DIR=phylogeny/$replicon/$out_prefix
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/${replicon}.treefile
  touch \${OUT_DIR}/iqtree.err
  touch \${OUT_DIR}/iqtree.log
  """  
}

process rec_removal_clonalframeml {
  // Tool: ClonalFrameML. 
  // Recombination correction on core SNPs tree.

  label 'highCPU'
  storeDir params.result
  debug false
  tag "ClonalFrameML on $replicon"  

  when:
    params.rec_removal_clonalframeml.todo == 1
  
  input:
    tuple val(replicon), path(core_snps_aln), path(core_genes_aln), path(treefile)

  output:
    tuple val(replicon), path("phylogeny/$replicon/clonalframeml/${replicon}.labelled_tree.newick")
    path("phylogeny/$replicon/clonalframeml/clonalframeml.err")
    path("phylogeny/$replicon/clonalframeml/clonalframeml.log")

  script:
  """
  OUT_DIR=phylogeny/$replicon/clonalframeml
  mkdir -p -m 777 \${OUT_DIR}
  
  ClonalFrameML $treefile $core_genes_aln \${OUT_DIR}/$replicon -xmfa_file true 1> \${OUT_DIR}/clonalframeml.log 2> \${OUT_DIR}/clonalframeml.err
  """

  stub:
  """
  OUT_DIR=phylogeny/$replicon/clonalframeml
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/${replicon}.labelled_tree.newick
  touch \${OUT_DIR}/clonalframeml.err
  touch \${OUT_DIR}/clonalframeml.log
  """  

}




