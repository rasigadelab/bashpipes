process assembly_flye {

  // Tool: Flye
  // De novo assembler for Oxford Nanopore long reads. 

  label 'denovo'
  storeDir (params.result)
  debug false
  tag "FLYE on $ont_reads.simpleName"

  when:
    params.assembly_flye.todo == 1

  input:
    tuple val(sample), path(ont_reads)
  
  output:
    tuple val(sample), path("genomes/$sample/flye/assembly.fasta"), emit : draft_assembly
    path("genomes/$sample/flye/flye.log")
    path("genomes/$sample/flye/flye.err")
    path("genomes/$sample/flye/assembly_info.txt")
    path("genomes/$sample/${sample}_assembly_raw.fasta")

  script:
  """
  OUT_DIR=genomes/$sample/flye
  mkdir -p -m 777 \${OUT_DIR}

  flye ${params.assembly_flye["ont_type"]} $ont_reads -o \${OUT_DIR} --threads $task.cpus 1> \${OUT_DIR}/flye.log 2> \${OUT_DIR}/flye.err 
  cp \${OUT_DIR}/assembly.fasta genomes/$sample/${sample}_assembly_raw.fasta

  echo "Completed Flye process for sample $sample"

  rm -rf work
  """

  stub:
  """
  OUT_DIR=genomes/$sample/flye
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/assembly.fasta
  touch \${OUT_DIR}/flye.log
  touch \${OUT_DIR}/flye.err
  touch \${OUT_DIR}/assembly_info.txt
  echo "Completed Flye process for sample $sample"

  rm -rf work
  """  
}

process assembly_spades {

  // Tool: SPAdes. 
  // De novo assembler for Illumina short reads.

  label 'denovo'
  storeDir params.result
  debug false
  tag "SPAdes on $sample"

  when:
    params.assembly_spades.todo == 1

  input:
    tuple val(sample), path(R1), path(R2)
  
  output:
    tuple val(sample), path("genomes/$sample/spades/contigs.fasta"), emit : draft_assembly
    path("genomes/$sample/spades/spades_1.log")
    path("genomes/$sample/spades/spades.err")
    path("genomes/$sample/${sample}_assembly_raw.fasta")

  script:
  """
  OUT_DIR=genomes/$sample/spades
  mkdir -p -m 777 \${OUT_DIR}

  spades.py -1 $R1 -2 $R2 ${params.assembly_spades["sample_type"]} -o \${OUT_DIR} -t $task.cpus 1> \${OUT_DIR}/spades_1.log 2> \${OUT_DIR}/spades.err
  cp \${OUT_DIR}/contigs.fasta genomes/$sample/${sample}_assembly_raw.fasta

  echo "Completed SPAdes process for sample $sample"

  rm -rf work
  """

  stub:
  """
  OUT_DIR=genomes/$sample/spades
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/contigs.fasta
  touch \${OUT_DIR}/spades_1.log
  touch \${OUT_DIR}/spades.err
  echo "Completed SPAdes process for sample $sample"

  rm -rf work
  """ 
}

process map_bowtie2 {
  // Tools: bowtie2 and samtools. 
  // Mapping step consists in creating an index of the genome on which reads will be mapped.
  // Then it aligns Illumina reads along draft assembly.
  // Then alignments are sorted and indexed.
  // Outputs a BAM and a BAI file per sample.

  label 'denovo'
  storeDir params.result
  debug false
  tag "Mapping on $sample"

  when:
    params.map_bowtie2.todo == 1

  input:
    tuple val(sample), path(draft_assembly), path(R1), path(R2)

  output:
    tuple val(sample), path("genomes/$sample/polish/${sample}.sorted.bam"), emit : sorted_bam_files
    path("genomes/$sample/polish/${sample}.sorted.bam.bai")
    path("genomes/$sample/polish/*.log")


  script:
  """
  OUT_DIR=genomes/$sample/polish
  SAMPLE_BAM=\${OUT_DIR}/${sample}.bam
  SAMPLE_BAM_SORTED=\${OUT_DIR}/${sample}.sorted.bam
  mkdir -p -m 777 \${OUT_DIR}

  bowtie2-build $draft_assembly \${OUT_DIR}/index &> \${OUT_DIR}/bowtie2.index.log
  bowtie2 -x \${OUT_DIR}/index -1 $R1 -2 $R2 -p $task.cpus 2>> \${OUT_DIR}/bowtie2.map.log | samtools view -bS - > \${SAMPLE_BAM}
  samtools sort \${SAMPLE_BAM} -o \${SAMPLE_BAM_SORTED} 2>> \${OUT_DIR}/samtools.sort.log
  samtools index \${SAMPLE_BAM_SORTED} 2>> \${OUT_DIR}/samtools.index.log

  rm -rf work
  """

  stub:
  """
  OUT_DIR=genomes/$sample/polish
  SAMPLE_BAM=\${OUT_DIR}/${sample}.bam
  SAMPLE_BAM_SORTED=\${OUT_DIR}/${sample}.sorted.bam
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/index.1.bt2
  touch \${OUT_DIR}/index.2.bt2
  touch \${OUT_DIR}/bowtie2.index.log
  touch \${SAMPLE_BAM}
  touch \${OUT_DIR}/bowtie2.map.log
  touch \${SAMPLE_BAM_SORTED}
  touch \${OUT_DIR}/samtools.sort.log
  touch \${SAMPLE_BAM_SORTED}.bai
  touch \${OUT_DIR}/samtools.index.log

  rm -rf work
  """ 
}

process polish_pilon {
  // Tool: pilon. 
  // Polishing step consists in correcting errors in draft assembly with Illumina reads alignment.

  label 'denovo'
  storeDir params.result
  debug false
  tag "Polishing of $sample"

  when:
    params.polish_pilon.todo == 1

  input:
    tuple val(sample), path(draft_assembly), path(sorted_bam)

  output:
    tuple val(sample), path("genomes/$sample/polish/${sample}_polished.fasta"), emit : polished_assembly
    path("genomes/$sample/polish/pilon.log")
    path("genomes/$sample/${sample}_polished.fasta")

  script:
  memory = (task.memory =~ /([^\ ]+)(.+)/)[0][1]
  """
  OUT_DIR=genomes/$sample/polish
  SAMPLE_POLISHED=\${OUT_DIR}/${sample}_polished
  mkdir -p -m 777 \${OUT_DIR}

  pilon -Xmx${memory}G --genome $draft_assembly --bam $sorted_bam ${params.polish_pilon["list_changes"]} --output \${SAMPLE_POLISHED} &> \${OUT_DIR}/pilon.log
  cp \${SAMPLE_POLISHED}.fasta genomes/$sample/${sample}_polished.fasta

  echo "Completed Pilon process for sample $sample"
  
  rm -rf work
  """

  stub:
  """
  OUT_DIR=genomes/$sample/polish
  SAMPLE_POLISHED=\${OUT_DIR}/${sample}_polished
  mkdir -p -m 777 \${OUT_DIR}
  touch \${SAMPLE_POLISHED}.fasta
  touch \${OUT_DIR}/pilon.log

  echo "Completed Pilon process for sample $sample"

  rm -rf work
  """

}

process qc_quast {
  // Tool: quast. 
  // Assembly QC step consists in checking quality of de novo Illumina assembly.

  label 'denovo'
  storeDir params.result
  debug false
  tag "Quast QC on $sample"

  when:
    params.qc_quast.todo == 1

  input:
    tuple val(sample), path(draft_assembly)

  output:
    path("genomes/$sample/quast/*")

  script:
  """
  OUT_DIR=genomes/$sample/quast
  mkdir -p -m 777 \${OUT_DIR}

  quast -o \${OUT_DIR} $draft_assembly 1> \${OUT_DIR}/quast.log 2> \${OUT_DIR}/quast.err

  rm -rf work
  """

  stub:
  """
  OUT_DIR=genomes/$sample/quast
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/quast.log
  touch \${OUT_DIR}/quast.err

  rm -rf work
  """  

}

process fixstart_circlator {
  // Tool: circlator. 
  // Circularization step consists in re-aligning contig sequences to origin of replication.

  label 'denovo'
  storeDir params.result
  debug false
  tag "Circlator on $sample"

  when:
    params.fixstart_circlator.todo == 1

  input:
    tuple val(sample), path(denovo_assembly)

  output:
    tuple val(sample), path("genomes/$sample/circlator/${sample}_realigned.fasta"), emit : realigned_assembly
    path("genomes/$sample/circlator/*.log")
    path("genomes/$sample/${sample}_realigned.fasta")
 
  script:
  """
  OUT_DIR=genomes/$sample/circlator
  mkdir -p -m 777 \${OUT_DIR}

  circlator fixstart $denovo_assembly \${OUT_DIR}/${sample}_realigned
  cp \${OUT_DIR}/${sample}_realigned.fasta genomes/$sample/${sample}_realigned.fasta

  rm -rf work
  """

  stub:
  """
  OUT_DIR=genomes/$sample/circlator
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/${sample}_realigned.fasta
  touch \${OUT_DIR}/circlator.log
  rm -rf work
  """  
}












