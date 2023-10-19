process rename_fasta {
  // Tool: Unix command
  // Gather and rename FASTA files so they can be used in Nf without problem

  label 'lowCPU'
  storeDir params.result
  debug false
  tag "FASTA gathering for $replicon"  

  when:
    params.rename_fasta.todo == 1
  
  input:
    tuple val(replicon), val(sample), path(fasta_file)

  output:
    tuple val(replicon), val(sample), path("phylogeny/$replicon/sequences/${sample}.fasta"), emit : fasta_renamed

  script:
  """
  OUT_DIR=phylogeny/$replicon/sequences
  mkdir -p -m 777 \${OUT_DIR}

  cat $fasta_file > \${OUT_DIR}/${sample}.fasta
  """

  stub:
  """
  OUT_DIR=phylogeny/$replicon/sequences
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/${sample}.fasta
  """

}

process distance_matrix_mash {
  // Tool: mash
  // Producing tsv file with mash distance between samples. 

  label 'highCPU'
  storeDir params.result
  debug false
  tag "Mash distance matrix for $replicon"  

  when:
    params.distance_matrix_mash.todo == 1
  
  input:
    tuple val(replicon), val(samples), path(fasta_files)

  output:
    tuple val(replicon), val(samples), emit : replicons_ch
    path("phylogeny/$replicon/mash/mash_dist.tsv")
    path("phylogeny/$replicon/mash/mash_analysis.log")
    path("phylogeny/$replicon/mash/mini_clusters.tsv")
    path("phylogeny/$replicon/mash/no_clonal_samples.tsv")

  script:
  """
  OUT_DIR=phylogeny/$replicon/mash
  mkdir -p -m 777 \${OUT_DIR}

  mash triangle $fasta_files > \${OUT_DIR}/mash_dist.tsv
  Rscript ${params.nfpath}/modules/mash_analysis.R -d \${OUT_DIR} -th ${params.distance_matrix_mash["mash_threshold"]} > \${OUT_DIR}/mash_analysis.log
  """

  stub:
  """
  OUT_DIR=phylogeny/$replicon/mash
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/mash_dist.tsv
  touch \${OUT_DIR}/mash_analysis.log
  touch \${OUT_DIR}/mini_clusters.tsv
  touch \${OUT_DIR}/no_clonal_samples.tsv
  """
}

process duplicate_masker_repeatmasker {
  // Tool: RepeatMasker
  // Random choice of a reference sample in subgroup (CAN BE CHANGED LATER)
  // Mask duplicated regions of the reference genome

  label 'repeatmasker'
  storeDir params.result
  debug true
  tag "RepeatMasker for $replicon - $subgroup"  

  when:
    params.duplicate_masker_repeatmasker.todo == 1
  
  input:
    tuple val(samples), val(replicon), val(subgroup)

  output:
    tuple val(replicon), val(samples), val(subgroup), path("phylogeny/$replicon/minicluster_$subgroup/repeatmasker/*.bed"), emit : replicons_ch
    path("phylogeny/$replicon/minicluster_$subgroup/repeatmasker/repeatmasker.log")
    path("phylogeny/$replicon/minicluster_$subgroup/repeatmasker/repeatmasker.err")
    


  script:
  """
  OUT_DIR=phylogeny/$replicon/minicluster_$subgroup/repeatmasker
  mkdir -p -m 777 \${OUT_DIR}

  SAMPLES_LIST=\$(echo "$samples" | tr -d "[]" | tr -d ",")
  ARRAY=(\$SAMPLES_LIST)
  SIZE=\${#ARRAY[@]}
  RANDOM_INDEX=\$((\$RANDOM % \$SIZE))
  RANDOMLY_CHOSEN_REF=\${ARRAY[\$RANDOM_INDEX]}

  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate repeatmasker
  RepeatMasker -dir \${OUT_DIR} ${params.result}/phylogeny/$replicon/sequences/\${RANDOMLY_CHOSEN_REF}.fasta 1> \${OUT_DIR}/repeatmasker.log 2> \${OUT_DIR}/repeatmasker.err
  rmsk2bed < \${OUT_DIR}/\${RANDOMLY_CHOSEN_REF}.fasta.out | bedops --merge - > \${OUT_DIR}/masked_\${RANDOMLY_CHOSEN_REF}.bed
  conda deactivate
  """

  stub:
  """
  OUT_DIR=phylogeny/$replicon/minicluster_$subgroup/repeatmasker
  mkdir -p -m 777 \${OUT_DIR}

  SAMPLES_LIST=\$(echo "$samples" | tr -d "[]" | tr -d ",")
  ARRAY=(\$SAMPLES_LIST)
  SIZE=\${#ARRAY[@]}
  RANDOM_INDEX=\$((\$RANDOM % \$SIZE))
  RANDOMLY_CHOSEN_REF=\${ARRAY[\$RANDOM_INDEX]}

  touch \${OUTDIR}/repeatmasker.log
  touch \${OUTDIR}/repeatmasker.err
  touch \${OUTDIR}/masked_\${RANDOMLY_CHOSEN_REF}.bed
  """
}

process create_input_tab {
  // Tool: homemade script. 
  // Creating tsv file 'input.tab' for snippy usage (Sample, Fasta_R1, Fasta_R2)

  label 'lowCPU'
  storeDir params.result
  debug false
  tag "Input tab for $replicon - $subgroup"  

  when:
    params.core_snps_snippy.todo == 1
  
  input:
    tuple val(replicon), val(samples), val(subgroup), path(masked_regions)

  output:
    tuple val(replicon), val(subgroup), path("phylogeny/$replicon/minicluster_$subgroup/repeatmasker/*.bed"), path("phylogeny/$replicon/minicluster_$subgroup/snippy/input.tab"), emit : input_tab

  script:
  """
  OUT_DIR=phylogeny/$replicon/minicluster_$subgroup/snippy
  mkdir -p -m 777 \${OUT_DIR}
  
  #samples = [1, 2, 3] But I want : samples = 1 2 3
  SAMPLES_LIST=\$(echo "$samples" | tr -d "[]" | tr -d ",")

  python3 ${params.nfpath}/modules/snippy_multi_list.py -d $params.result -l \$SAMPLES_LIST -o \${OUT_DIR}
  """

  stub:
  """
  OUT_DIR=phylogeny/$replicon/minicluster_$subgroup/snippy
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
  tag "Snippy on $replicon - $subgroup" 

  when:
    params.core_snps_snippy.todo == 1
  
  input:
    tuple val(replicon), val(subgroup), path(masked_regions), path(input_tab)

  output:
    path("phylogeny/$replicon/minicluster_$subgroup/snippy/*/snps.aligned.fa")
    path("phylogeny/$replicon/minicluster_$subgroup/snippy/*/snps.log")
    path("phylogeny/$replicon/minicluster_$subgroup/snippy/*/snps.subs.vcf")
    path("phylogeny/$replicon/minicluster_$subgroup/snippy/core.tab")
    path("phylogeny/$replicon/minicluster_$subgroup/snippy/core.txt")
    path("phylogeny/$replicon/minicluster_$subgroup/snippy/core.vcf")
    path("phylogeny/$replicon/minicluster_$subgroup/snippy/snippy.log")
    path("phylogeny/$replicon/minicluster_$subgroup/snippy/snippy.err")
    path("phylogeny/$replicon/minicluster_$subgroup/snippy/snp_matrix.tsv")
    path("phylogeny/$replicon/minicluster_$subgroup/snippy/core_without_ref.aln")

  script:
  """
  OUT_DIR=phylogeny/$replicon/minicluster_$subgroup/snippy
  mkdir -p -m 777 \${OUT_DIR}

  REF_SAMPLE=\$(echo $masked_regions | cut -d '_' -f 2 | cut -d '.' -f 1)
  REF_PATH=${params.result}/phylogeny/$replicon/sequences/\${REF_SAMPLE}.fasta
  snippy-multi $input_tab --ref \${REF_PATH} --cpus $task.cpus --force --mincov 30 > \${OUT_DIR}/snippy_commands.sh
  sed -i -e "s|snippy-core --ref '|snippy-core --prefix \${OUT_DIR}/core --mask $masked_regions --ref '|g" \${OUT_DIR}/snippy_commands.sh
  sh \${OUT_DIR}/snippy_commands.sh 1> \${OUT_DIR}/snippy.log 2> \${OUT_DIR}/snippy.err
  
  snp-dists \${OUT_DIR}/core.full.aln > \${OUT_DIR}/snp_matrix.tsv

  # Remove last lines because it's the reference sequence and we don't want it in phylogeny
  sed -n '/>Reference/q;p' \${OUT_DIR}/core.full.aln > \${OUT_DIR}/core_without_ref.aln
  """

  stub:
  """
  OUT_DIR=phylogeny/$replicon/minicluster_$subgroup/snippy
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

process duplicate_masker_phylogeny {
  // Tool: RepeatMasker
  // Choice of reference genome: manually (cf params)
  // Mask duplicated regions of the reference genome

  label 'repeatmasker'
  storeDir params.result
  debug true
  tag "RepeatMasker for $replicon"  

  when:
    params.duplicate_masker_phylogeny.todo == 1
  
  input:
    tuple val(replicon), val(samples)

  output:
    tuple val(replicon), val(samples), path("phylogeny/$replicon/phylogenetic_tree/repeatmasker/*.bed"), emit : replicons_ch
    path("phylogeny/$replicon/phylogenetic_tree/repeatmasker/repeatmasker.log")
    path("phylogeny/$replicon/phylogenetic_tree/repeatmasker/repeatmasker.err")
    
  script:
  """
  OUT_DIR=phylogeny/$replicon/phylogenetic_tree/repeatmasker
  mkdir -p -m 777 \${OUT_DIR}

  REF=${params.duplicate_masker_phylogeny["reference_sample"]}

  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate repeatmasker
  RepeatMasker -dir \${OUT_DIR} ${params.result}/phylogeny/$replicon/sequences/\${REF}.fasta 1> \${OUT_DIR}/repeatmasker.log 2> \${OUT_DIR}/repeatmasker.err
  rmsk2bed < \${OUT_DIR}/\${REF}.fasta.out | bedops --merge - > \${OUT_DIR}/masked_\${REF}.bed
  conda deactivate
  """

  stub:
  """
  OUT_DIR=phylogeny/$replicon/phylogenetic_tree/repeatmasker
  mkdir -p -m 777 \${OUT_DIR}

  REF=${params.duplicate_masker_phylogeny["reference_sample"]}

  touch \${OUTDIR}/repeatmasker.log
  touch \${OUTDIR}/repeatmasker.err
  touch \${OUTDIR}/masked_\${REF}.bed
  """
}

process create_input_tab_phylogeny {
  // Tool: homemade script. 
  // Creating tsv file 'input.tab' for snippy usage (Sample, Fasta_R1, Fasta_R2)

  label 'lowCPU'
  storeDir params.result
  debug false
  tag "Input tab for $replicon"  

  when:
    params.core_snps_snippy_phylogeny.todo == 1
  
  input:
    tuple val(replicon), val(samples), path(masked_regions)

  output:
    tuple val(replicon), path("phylogeny/$replicon/phylogenetic_tree/repeatmasker/*.bed"), path("phylogeny/$replicon/phylogenetic_tree/snippy/input.tab"), emit : input_tab

  script:
  """
  OUT_DIR=phylogeny/$replicon/phylogenetic_tree/snippy
  mkdir -p -m 777 \${OUT_DIR}
  
  #samples = [1, 2, 3] But I want : samples = 1 2 3
  SAMPLES_LIST=\$(echo "$samples" | tr -d "[]" | tr -d ",")

  python3 ${params.nfpath}/modules/snippy_multi_list.py -d $params.result -l \$SAMPLES_LIST -o \${OUT_DIR}
  """

  stub:
  """
  OUT_DIR=phylogeny/$replicon/phylogenetic_tree/snippy
  mkdir -p -m 777 \${OUT_DIR}
  #samples = [1, 2, 3] But I want : samples = 1 2 3
  SAMPLES_LIST=\$(echo "$samples" | tr -d "[]" | tr -d ",")
  python3 /home/afischer/11.Phylogeny_test/snippy_multi_list.py -d $params.result -l \$SAMPLES_LIST -o \${OUT_DIR}
  """

}

process core_snps_snippy_phylogeny {
  // Tool: snippy. 
  // Illumina FASTQ mapping on core_genome reference FASTA file.
  // Then, SNP calling on mapping results. 

  label 'snippy'
  storeDir params.result
  debug false
  tag "Snippy on $replicon" 

  when:
    params.core_snps_snippy_phylogeny.todo == 1
  
  input:
    tuple val(replicon), path(masked_regions), path(input_tab)

  output:
    tuple val(replicon), path("phylogeny/$replicon/phylogenetic_tree/snippy/clean.full.aln"), emit : core_snps_aln
    path("phylogeny/$replicon/phylogenetic_tree/snippy/*/snps.aligned.fa")
    path("phylogeny/$replicon/phylogenetic_tree/snippy/*/snps.log")
    path("phylogeny/$replicon/phylogenetic_tree/snippy/*/snps.subs.vcf")
    path("phylogeny/$replicon/phylogenetic_tree/snippy/core.tab")
    path("phylogeny/$replicon/phylogenetic_tree/snippy/core.txt")
    path("phylogeny/$replicon/phylogenetic_tree/snippy/core.vcf")
    path("phylogeny/$replicon/phylogenetic_tree/snippy/snippy.log")
    path("phylogeny/$replicon/phylogenetic_tree/snippy/snippy.err")
    path("phylogeny/$replicon/phylogenetic_tree/snippy/snp_matrix.tsv")

  script:
  """
  OUT_DIR=phylogeny/$replicon/phylogenetic_tree/snippy
  mkdir -p -m 777 \${OUT_DIR}

  REF_SAMPLE=\$(echo $masked_regions | cut -d '_' -f 2 | cut -d '.' -f 1)
  REF_PATH=${params.result}/phylogeny/$replicon/sequences/\${REF_SAMPLE}.fasta
  snippy-multi $input_tab --ref \${REF_PATH} --cpus $task.cpus --force --mincov 30 > \${OUT_DIR}/snippy_commands.sh
  sed -i -e "s|snippy-core --ref '|snippy-core --prefix \${OUT_DIR}/core --mask $masked_regions --ref '|g" \${OUT_DIR}/snippy_commands.sh
  sh \${OUT_DIR}/snippy_commands.sh 1> \${OUT_DIR}/snippy.log 2> \${OUT_DIR}/snippy.err
  
  # Create SNP distance matrix between samples
  snp-dists \${OUT_DIR}/core.full.aln > \${OUT_DIR}/snp_matrix.tsv

  # Remove last lines because it's the reference sequence and we don't want it in phylogeny
  sed -n '/>Reference/q;p' \${OUT_DIR}/core.full.aln > \${OUT_DIR}/core_without_ref.aln

  # Clean SNP alignment with Snippy command (remove weird characters)
  snippy-clean_full_aln \${OUT_DIR}/core_without_ref.aln > \${OUT_DIR}/clean.full.aln
  """

  stub:
  """
  OUT_DIR=phylogeny/$replicon/phylogenetic_tree/snippy
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

  label 'highCPU'
  storeDir params.result
  debug false
  tag "Iqtree on $replicon"  

  when:
    params.snps_tree_iqtree.todo == 1
  
  input:
    tuple val(replicon), path(core_snps_aln)

  output:
    tuple val(replicon), path("phylogeny/$replicon/phylogenetic_tree/snippy/clean.full.aln"), path("phylogeny/$replicon/phylogenetic_tree/iqtree/${replicon}.treefile"), emit : treefile
    path("phylogeny/$replicon/phylogenetic_tree/iqtree/${replicon}.iqtree")
    path("phylogeny/$replicon/phylogenetic_tree/iqtree/iqtree.err")
    path("phylogeny/$replicon/phylogenetic_tree/iqtree/iqtree.log")
    path("phylogeny/$replicon/phylogenetic_tree/iqtree/${replicon}.alninfo")
    path("phylogeny/$replicon/phylogenetic_tree/iqtree/${replicon}.bionj")
    path("phylogeny/$replicon/phylogenetic_tree/iqtree/${replicon}.log")
    path("phylogeny/$replicon/phylogenetic_tree/iqtree/${replicon}.mldist")

  script:
  """
  OUT_DIR=phylogeny/$replicon/phylogenetic_tree/iqtree
  mkdir -p -m 777 \${OUT_DIR}

  iqtree -s $core_snps_aln -m ${params.snps_tree_iqtree["model"]} --prefix \${OUT_DIR}/$replicon -T $task.cpus -alninfo 1> \${OUT_DIR}/iqtree.log 2> \${OUT_DIR}/iqtree.err
  """

  stub:
  """
  OUT_DIR=phylogeny/$replicon/phylogenetic_tree/iqtree
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/${replicon}.treefile
  touch \${OUT_DIR}/iqtree.err
  touch \${OUT_DIR}/iqtree.log
  touch \${OUT_DIR}/${replicon}.iqtree
  touch \${OUT_DIR}/${replicon}.alninfo
  touch \${OUT_DIR}/${replicon}.bionj
  touch \${OUT_DIR}/${replicon}.log
  touch \${OUT_DIR}/${replicon}.mldist
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
    tuple val(replicon), path(core_snps_aln), path(treefile)

  output:
    tuple val(replicon), path("phylogeny/$replicon/phylogenetic_tree/snippy/clean.full.aln"), path("phylogeny/$replicon/phylogenetic_tree/clonalframeml/${replicon}.labelled_tree.newick"), emit : tree_without_rec
    path("phylogeny/$replicon/phylogenetic_tree/clonalframeml/clonalframeml.err")
    path("phylogeny/$replicon/phylogenetic_tree/clonalframeml/clonalframeml.log")
    path("phylogeny/$replicon/phylogenetic_tree/clonalframeml/${replicon}.em.txt")
    path("phylogeny/$replicon/phylogenetic_tree/clonalframeml/${replicon}.importation_status.txt")
    path("phylogeny/$replicon/phylogenetic_tree/clonalframeml/${replicon}.ML_sequence.fasta")
    path("phylogeny/$replicon/phylogenetic_tree/clonalframeml/${replicon}.position_cross_reference.txt")

  script:
  """
  OUT_DIR=phylogeny/$replicon/phylogenetic_tree/clonalframeml
  mkdir -p -m 777 \${OUT_DIR}
  
  ClonalFrameML $treefile $core_snps_aln \${OUT_DIR}/$replicon 1> \${OUT_DIR}/clonalframeml.log 2> \${OUT_DIR}/clonalframeml.err
  """

  stub:
  """
  OUT_DIR=phylogeny/$replicon/phylogenetic_tree/clonalframeml
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/${replicon}.labelled_tree.newick
  touch \${OUT_DIR}/clonalframeml.err
  touch \${OUT_DIR}/clonalframeml.log
  touch \${OUT_DIR}/${replicon}.em.txt
  touch \${OUT_DIR}/${replicon}.importation_status.txt
  touch \${OUT_DIR}/${replicon}.ML_sequence.fasta
  touch \${OUT_DIR}/${replicon}.position_cross_reference.txt
  """  

}

process dating_treetime {
  // Tool: Treetime. 
  // Phylogenetic dating of a tree.

  label 'highCPU'
  storeDir params.result
  debug false
  tag "Treetime on $replicon"  

  when:
    params.dating_treetime.todo == 1
  
  input:
    tuple val(replicon), path(core_snps_aln), path(treefile)

  output:
    path("phylogeny/$replicon/phylogenetic_tree/treetime/treetime.err")
    path("phylogeny/$replicon/phylogenetic_tree/treetime/treetime.log")
    path("phylogeny/$replicon/phylogenetic_tree/treetime/out/*")

  script:
  """
  OUT_DIR=phylogeny/$replicon/phylogenetic_tree/treetime
  mkdir -p -m 777 \${OUT_DIR}
  
  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate treetime
  treetime --aln $core_snps_aln --tree $treefile --dates $params.result/${replicon}_sampling_list.csv --outdir \${OUT_DIR}/out --plot-tree treetime.png 1> \${OUT_DIR}/treetime.log 2> \${OUT_DIR}/treetime.err
  conda deactivate
  """

  stub:
  """
  OUT_DIR=phylogeny/$replicon/phylogenetic_tree/treetime/out
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/treetime.err
  touch \${OUT_DIR}/treetime.log
  """  

}


