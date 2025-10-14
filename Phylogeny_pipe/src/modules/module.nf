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
}

process distance_matrix_mash {
  // Tool: mash
  // Producing tsv file with mash distance between samples. 

  label 'mash'
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

  script:
  """
  OUT_DIR=phylogeny/$replicon/mash
  mkdir -p -m 777 \${OUT_DIR}

  mash triangle $fasta_files > \${OUT_DIR}/mash_dist.tsv
  """
}

process distance_matrix_analysis {
  // Tool: R
  // Producing tsv file with cluster of samples based on MASH distance

  label 'r_analysis'
  storeDir params.result
  debug false
  tag "Mash distance clustering for $replicon"  

  when:
    params.distance_matrix_analysis.todo == 1
  
  input:
    tuple val(replicon), val(samples)

  output:
    tuple val(replicon), val(samples), emit : replicons_ch
    path("phylogeny/$replicon/mash/mash_analysis.log")
    path("phylogeny/$replicon/mash/mini_clusters.tsv")
    path("phylogeny/$replicon/mash/no_clonal_samples.tsv")
    path("phylogeny/$replicon/mash/reference_for_phylogeny.tsv")

  script:
  """
  OUT_DIR=phylogeny/$replicon/mash
  mkdir -p -m 777 \${OUT_DIR}

  Rscript ${params.nfpath}/modules/mash_analysis.R -d ${params.result}/\${OUT_DIR} -t ${params.distance_matrix_analysis["mash_threshold"]} > \${OUT_DIR}/mash_analysis.log
  """
}

process map_bowtie2 {
  // Tools: bowtie2 and samtools. 
  // Mapping step consists in creating an index of the genome on which reads will be mapped.
  // Then it aligns Illumina reads along draft assembly.
  // Then alignments are sorted and indexed.
  // Outputs a BAM and a BAI file per sample.

  label 'map_and_sort'
  storeDir params.result
  debug false
  tag "Mapping for $replicon - $subgroup"

  when:
    params.map_bowtie2.todo == 1

  input:
    tuple val(samples), val(replicon), val(subgroup), val(ref)

  output:
    tuple val(samples), val(replicon), val(subgroup), val(ref), path("phylogeny/$replicon/minicluster_$subgroup/polish/${ref}.sorted.bam"), emit : sorted_bam_files
    path("phylogeny/$replicon/minicluster_$subgroup/polish/${ref}.sorted.bam.bai")
    path("phylogeny/$replicon/minicluster_$subgroup/polish/bowtie2.index.log")
    path("phylogeny/$replicon/minicluster_$subgroup/polish/bowtie2.map.log")
    path("phylogeny/$replicon/minicluster_$subgroup/polish/samtools.index.log")
    path("phylogeny/$replicon/minicluster_$subgroup/polish/samtools.sort.log")

  script:
  memory = (task.memory =~ /([^\ ]+)(.+)/)[0][1]
  """
  OUT_DIR=phylogeny/$replicon/minicluster_$subgroup/polish
  SAMPLE_BAM=\${OUT_DIR}/${ref}.bam
  SAMPLE_BAM_SORTED=\${OUT_DIR}/${ref}.sorted.bam
  R1=$params.result/genomes/$ref/trimmomatic/${ref}_R1_paired.trimmed.fastq.gz
  R2=$params.result/genomes/$ref/trimmomatic/${ref}_R2_paired.trimmed.fastq.gz
  ASSEMBLY=${params.result}/phylogeny/$replicon/sequences/${ref}.fasta
  mkdir -p -m 777 \${OUT_DIR}

  bowtie2-build \${ASSEMBLY} \${OUT_DIR}/index &> \${OUT_DIR}/bowtie2.index.log
  bowtie2 -x \${OUT_DIR}/index -1 \$R1 -2 \$R2 -p $task.cpus 2>> \${OUT_DIR}/bowtie2.map.log | samtools view -bS > \${SAMPLE_BAM}
  
  samtools sort \${SAMPLE_BAM} -o \${SAMPLE_BAM_SORTED} -@ $task.cpus -m ${memory}G 2>> \${OUT_DIR}/samtools.sort.log
  samtools index \${SAMPLE_BAM_SORTED} -@ $task.cpus 2>> \${OUT_DIR}/samtools.index.log
  """
}

process polish_pilon {
  // Tool: pilon. 
  // Polishing step consists in correcting errors in draft assembly with Illumina reads alignment.

  label 'pilon'
  storeDir params.result
  debug false
  tag "Polishing for $replicon - $subgroup"

  when:
    params.polish_pilon.todo == 1

  input:
    tuple val(samples), val(replicon), val(subgroup), val(ref), path(sorted_bam)

  output:
    tuple val(samples), val(replicon), val(subgroup), val(ref), path("phylogeny/$replicon/minicluster_$subgroup/polish/${ref}_polished.fasta"), emit : polished_reference
    path("phylogeny/$replicon/minicluster_$subgroup/polish/pilon.log")
    path("phylogeny/$replicon/minicluster_$subgroup/polish/${ref}_polished.changes")

  script:
  memory = (task.memory =~ /([^\ ]+)(.+)/)[0][1]
  """
  OUT_DIR=phylogeny/$replicon/minicluster_$subgroup/polish
  SAMPLE_POLISHED=\${OUT_DIR}/${ref}_polished
  mkdir -p -m 777 \${OUT_DIR}
  ASSEMBLY=${params.result}/phylogeny/$replicon/sequences/${ref}.fasta
  SAMPLE_POLISHED=\${OUT_DIR}/${ref}_polished

  java -Xmx${memory}G -jar /pilon/pilon.jar --genome \${ASSEMBLY} --bam $sorted_bam ${params.polish_pilon["list_changes"]} --output \${SAMPLE_POLISHED} &> \${OUT_DIR}/pilon.log
  """
}

process duplicate_masker_repeatmasker {
  // Tool: RepeatMasker
  // Mask duplicated regions of the reference genome

  label 'repeatmasker'
  storeDir params.result
  debug true
  tag "RepeatMasker for $replicon - $subgroup"  

  when:
    params.duplicate_masker_repeatmasker.todo == 1
  
  input:
    tuple val(samples), val(replicon), val(subgroup), val(ref), path(polished_ref)

  output:
    tuple val(replicon), val(samples), val(subgroup), val(ref), path("phylogeny/$replicon/minicluster_$subgroup/repeatmasker/${ref}_polished.fasta.out"), emit : replicons_ch
    path("phylogeny/$replicon/minicluster_$subgroup/repeatmasker/repeatmasker.log")
    path("phylogeny/$replicon/minicluster_$subgroup/repeatmasker/repeatmasker.err")
    
  script:
  """
  OUT_DIR=phylogeny/$replicon/minicluster_$subgroup/repeatmasker
  mkdir -p -m 777 \${OUT_DIR}

  RepeatMasker -dir \${OUT_DIR} $polished_ref 1> \${OUT_DIR}/repeatmasker.log 2> \${OUT_DIR}/repeatmasker.err
  """
}

process duplicate_masker_bedops {
  // Tool: Bedops
  // Write duplicated regions of the reference genome into BED file

  label 'bedops'
  storeDir params.result
  debug true
  tag "Rmsk2bed for $replicon - $subgroup"  

  when:
    params.duplicate_masker_bedops.todo == 1
  
  input:
    tuple val(replicon), val(samples), val(subgroup), val(ref), path(repeat_regions)

  output:
    tuple val(replicon), val(samples), val(subgroup), path("phylogeny/$replicon/minicluster_$subgroup/repeatmasker/masked_${ref}.bed"), emit : replicons_ch
    
  script:
  """
  OUT_DIR=phylogeny/$replicon/minicluster_$subgroup/repeatmasker
  mkdir -p -m 777 \${OUT_DIR}

  rmsk2bed < $repeat_regions | bedops --merge - > \${OUT_DIR}/masked_${ref}.bed #\${OUT_DIR}/${ref}_polished.fasta.out
  """
}

process create_input_tab {
  // Tool: homemade script. 
  // Creating tsv file 'input.tab' for snippy usage (Sample, Fasta_R1, Fasta_R2)

  label 'python_tab_creation'
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
    tuple val(replicon), val(subgroup), path("phylogeny/$replicon/minicluster_$subgroup/snippy/clean.full.aln"), emit : clean_core_snps
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
    path("phylogeny/$replicon/minicluster_$subgroup/snippy/clean.full.withref.aln")
    path("phylogeny/$replicon/minicluster_$subgroup/snippy/core.ref.fa")

  script:
  """
  OUT_DIR=phylogeny/$replicon/minicluster_$subgroup/snippy
  mkdir -p -m 777 \${OUT_DIR}

  REF_SAMPLE=\$(echo $masked_regions | cut -d '_' -f 2 | cut -d '.' -f 1)
  REF_PATH=${params.result}/phylogeny/$replicon/minicluster_$subgroup/polish/\${REF_SAMPLE}_polished.fasta
  snippy-multi $input_tab --ref \${REF_PATH} --cpus $task.cpus --force --mincov ${params.core_snps_snippy["mincov"]} > \${OUT_DIR}/snippy_commands.sh
  sed -i -e "s|snippy-core --ref '|snippy-core --prefix \${OUT_DIR}/core --mask $masked_regions --ref '|g" \${OUT_DIR}/snippy_commands.sh
  sh \${OUT_DIR}/snippy_commands.sh 1> \${OUT_DIR}/snippy.log 2> \${OUT_DIR}/snippy.err
  
  snp-dists \${OUT_DIR}/core.aln > \${OUT_DIR}/snp_matrix.tsv

  # Remove last lines because it's the reference sequence and we don't want it in phylogeny
  sed -n '/>Reference/q;p' \${OUT_DIR}/core.full.aln > \${OUT_DIR}/core_without_ref.aln

  # Clean SNP alignment with Snippy command (remove weird characters)
  snippy-clean_full_aln \${OUT_DIR}/core_without_ref.aln > \${OUT_DIR}/clean.full.aln
  snippy-clean_full_aln \${OUT_DIR}/core.full.aln > \${OUT_DIR}/clean.full.withref.aln
  """
}

process vc_recombination_analysis_gubbins {
  // Tool: Gubbins
  // Run Iqtree based on SNP alignment and predict recombination blocks

  label 'gubbins'
  storeDir params.result
  debug false
  tag "Gubbins for $replicon - $subgroup"  

  when:
    params.vc_recombination_analysis_gubbins.todo == 1
  
  input:
    tuple val(replicon), val(subgroup), path(core_snps_aln)

  output:
    tuple val(replicon), val(subgroup), path("phylogeny/$replicon/minicluster_$subgroup/gubbins/${replicon}_core_masked.strict_core_after_masking.fasta"), emit : aln_without_rec
    path("phylogeny/$replicon/minicluster_$subgroup/gubbins/${replicon}.final_tree.tre")
    path("phylogeny/$replicon/minicluster_$subgroup/gubbins/${replicon}.log")
    path("phylogeny/$replicon/minicluster_$subgroup/gubbins/gubbins.err")
    path("phylogeny/$replicon/minicluster_$subgroup/gubbins/gubbins.log")
    path("phylogeny/$replicon/minicluster_$subgroup/gubbins/${replicon}.per_branch_statistics.csv")
    path("phylogeny/$replicon/minicluster_$subgroup/gubbins/${replicon}.recombination_predictions.gff")
    path("phylogeny/$replicon/minicluster_$subgroup/gubbins/*")
    path("phylogeny/$replicon/minicluster_$subgroup/gubbins/clean.recombination_masked.aln")

  script:
  """
  OUT_DIR=phylogeny/$replicon/minicluster_$subgroup/gubbins
  mkdir -p -m 777 \${OUT_DIR}

  # Count number of samples analyzed
  NB_SAMPLES=\$(grep -o ">" $core_snps_aln | wc -l)
  # If only 2 samples
  if [ \$NB_SAMPLES = 2 ]
  then
    run_gubbins.py --pairwise $core_snps_aln --threads $task.cpus --filter-percentage 50.0 \
                   --prefix \${OUT_DIR}/$replicon 1> \${OUT_DIR}/gubbins.log 2> \${OUT_DIR}/gubbins.err
  else
    run_gubbins.py --tree-builder iqtree --first-tree-builder iqtree-fast --model GTR --model-fitter raxml $core_snps_aln \
                  --threads $task.cpus --filter-percentage 50.0 --prefix \${OUT_DIR}/$replicon 1> \${OUT_DIR}/gubbins.log 2> \${OUT_DIR}/gubbins.err
  fi
  # Masking recombination areas in input alignment
  mask_gubbins_aln.py --aln $core_snps_aln --gff \${OUT_DIR}/${replicon}.recombination_predictions.gff --out \${OUT_DIR}/clean.recombination_masked.aln
  
  # Filter only core SNPs
  SNIPPY_DIR=$params.result/phylogeny/$replicon/minicluster_$subgroup/snippy
  REFERENCE="\${SNIPPY_DIR}/core.ref.fa"
  CLEAN_WITH_REF="\${SNIPPY_DIR}/clean.full.withref.aln"
  SNIPPY_CORE_VCF="\${SNIPPY_DIR}/core.vcf"
  GUBBINS_MASKED_ALN="\${OUT_DIR}/clean.recombination_masked.aln"
  OUT_PREFIX="\${OUT_DIR}/${replicon}_core_masked"

  python3 ${params.nfpath}/modules/extract_core_snps_after_gubbins.py \
          --reference \$REFERENCE \
          --clean_with_ref \$CLEAN_WITH_REF \
          --snippy_core_vcf \$SNIPPY_CORE_VCF \
          --gubbins_masked_aln \$GUBBINS_MASKED_ALN \
          --out_prefix \$OUT_PREFIX
  """  
}

process vc_final_snp_matrix {
  // Tool: snp-dists
  // Output SNP matrix based on alignment file without recombination

  label 'snp_dists'
  storeDir params.result
  debug false
  tag "SNP matrix for $replicon - $subgroup"  

  when:
    params.vc_final_snp_matrix.todo == 1
  
  input:
    tuple val(replicon), val(subgroup), path(core_snps_aln)

  output:
    path("phylogeny/$replicon/minicluster_$subgroup/gubbins/core_snp_matrix.tsv")

  script:
  """
  OUT_DIR=phylogeny/$replicon/minicluster_$subgroup/gubbins
  mkdir -p -m 777 \${OUT_DIR}

  snp-dists $core_snps_aln > \${OUT_DIR}/core_snp_matrix.tsv
  """  
}

process ref_phylogeny {
  // Tool: RepeatMasker
  // Choice of reference genome: manually (cf params)
  // Mask duplicated regions of the reference genome

  label 'lowCPU'
  storeDir params.result
  debug true
  tag "Reference for $replicon"  

  when:
    params.ref_phylogeny.todo == 1
  
  input:
    tuple val(replicon), val(samples), val(ref)

  output:
    tuple val(replicon), val(samples), path("phylogeny/$replicon/phylogenetic_tree/reference/masked_${ref}.bed"), emit : replicons_ch
    path("phylogeny/$replicon/phylogenetic_tree/reference/${ref}_polished.fasta")
    
  script:
  """
  OUT_DIR=phylogeny/$replicon/phylogenetic_tree/reference
  mkdir -p -m 777 \${OUT_DIR}

  FASTA_REF=\$(find ${params.result} -type f -iname "${ref}_polished.fasta")
  MASK_REF=\$(find ${params.result} -type f -iname "masked_${ref}.bed")
  cp \${FASTA_REF} \${OUT_DIR}
  cp \${MASK_REF} \${OUT_DIR}
  """
}

process create_input_tab_phylogeny {
  // Tool: homemade script. 
  // Creating tsv file 'input.tab' for snippy usage (Sample, Fasta_R1, Fasta_R2)

  label 'python_tab_creation'
  storeDir params.result
  debug false
  tag "Input tab for $replicon"  

  when:
    params.core_snps_snippy_phylogeny.todo == 1
  
  input:
    tuple val(replicon), val(samples), path(mask_file)

  output:
    tuple val(replicon), path("phylogeny/$replicon/phylogenetic_tree/reference/*.bed"), path("phylogeny/$replicon/phylogenetic_tree/snippy/input.tab"), emit : input_tab

  script:
  """
  OUT_DIR=phylogeny/$replicon/phylogenetic_tree/snippy
  mkdir -p -m 777 \${OUT_DIR}
  
  #samples = [1, 2, 3] But I want : samples = 1 2 3
  SAMPLES_LIST=\$(echo "$samples" | tr -d "[]" | tr -d ",")

  python3 ${params.nfpath}/modules/snippy_multi_list.py -d $params.result -l \$SAMPLES_LIST -o \${OUT_DIR}
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
  REF_PATH=${params.result}/phylogeny/$replicon/phylogenetic_tree/reference/\${REF_SAMPLE}_polished.fasta
  snippy-multi $input_tab --ref \${REF_PATH} --cpus $task.cpus --force --mincov ${params.core_snps_snippy_phylogeny["mincov"]} > \${OUT_DIR}/snippy_commands.sh
  sed -i -e "s|snippy-core --ref '|snippy-core --prefix \${OUT_DIR}/core --mask $masked_regions --ref '|g" \${OUT_DIR}/snippy_commands.sh
  sh \${OUT_DIR}/snippy_commands.sh 1> \${OUT_DIR}/snippy.log 2> \${OUT_DIR}/snippy.err
  
  # Create SNP distance matrix between samples
  snp-dists \${OUT_DIR}/core.full.aln > \${OUT_DIR}/snp_matrix.tsv

  # Remove last lines because it's the reference sequence and we don't want it in phylogeny
  sed -n '/>Reference/q;p' \${OUT_DIR}/core.full.aln > \${OUT_DIR}/core_without_ref.aln

  # Clean SNP alignment with Snippy command (remove weird characters)
  snippy-clean_full_aln \${OUT_DIR}/core_without_ref.aln > \${OUT_DIR}/clean.full.aln
  """
}

process snps_tree_iqtree {
  // Tool: iqtree. 
  // Tree construction of a specific replicon based on its core snps.

  label 'iqtree'
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
}

process rec_removal_clonalframeml {
  // Tool: ClonalFrameML. 
  // Recombination correction on core SNPs tree.

  label 'clonalframeml'
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
}

process dating_treetime {
  // Tool: Treetime. 
  // Phylogenetic dating of a tree.

  label 'treetime'
  storeDir params.result
  debug false
  tag "Treetime on $replicon"  

  when:
    params.dating_treetime.todo == 1
  
  input:
    tuple val(replicon), path(core_snps_aln), path(treefile)

  output:
    tuple val(replicon), path("phylogeny/$replicon/phylogenetic_tree/snippy/clean.full.aln"), path("phylogeny/$replicon/phylogenetic_tree/clonalframeml/${replicon}.labelled_tree.newick"), emit : snippy_aln
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
}

process recombination_analysis_gubbins {
  // Tool: Gubbins
  // Run Iqtree based on SNP alignment and predict recombination blocks

  label 'gubbins'
  storeDir params.result
  debug false
  tag "Gubbins for $replicon"  

  when:
    params.recombination_analysis_gubbins.todo == 1
  
  input:
    tuple val(replicon), path(core_snps_aln), path(treefile)

  output:
    tuple val(replicon), path("phylogeny/$replicon/phylogenetic_tree/gubbins/clean.recombination_masked.aln"), emit : aln_without_rec
    path("phylogeny/$replicon/phylogenetic_tree/gubbins/${replicon}.final_tree.tre")
    path("phylogeny/$replicon/phylogenetic_tree/gubbins/${replicon}.log")
    path("phylogeny/$replicon/phylogenetic_tree/gubbins/gubbins.err")
    path("phylogeny/$replicon/phylogenetic_tree/gubbins/gubbins.log")
    path("phylogeny/$replicon/phylogenetic_tree/gubbins/${replicon}.per_branch_statistics.csv")
    path("phylogeny/$replicon/phylogenetic_tree/gubbins/${replicon}.recombination_predictions.gff")
    path("phylogeny/$replicon/phylogenetic_tree/gubbins/*")

  script:
  """
  OUT_DIR=phylogeny/$replicon/phylogenetic_tree/gubbins
  mkdir -p -m 777 \${OUT_DIR}

  run_gubbins.py --tree-builder iqtree --first-tree-builder iqtree-fast --model GTR $core_snps_aln \
          --threads $task.cpus --filter-percentage 50.0 --prefix \${OUT_DIR}/$replicon 1> \${OUT_DIR}/gubbins.log 2> \${OUT_DIR}/gubbins.err
  
  # Masking recombination areas in input alignment
  mask_gubbins_aln.py --aln $core_snps_aln --gff \${OUT_DIR}/${replicon}.recombination_predictions.gff --out \${OUT_DIR}/clean.recombination_masked.aln
  """  
}


