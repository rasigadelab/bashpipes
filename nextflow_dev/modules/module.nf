process assembly_flye {

  // Tool: Flye
  // De novo assembler for Oxford Nanopore long reads. 

  label 'denovo'
  storeDir (params.result)
  debug true
  tag "FLYE on $ont_reads.simpleName"

  when:
    params.assembly_flye.todo == 1

  input:
    tuple val(sample), path(ont_reads)
  
  output:
    tuple val(sample), path("genomes/$sample/flye/assembly.fasta"), emit : draft_assembly
    path("genomes/$sample/flye/assembly.fasta")
    path("genomes/$sample/flye/flye.log")
    path("genomes/$sample/flye/flye.err")
    path("genomes/$sample/flye/assembly_info.txt")

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
  debug true
  tag "SPAdes on $sample"

  when:
    params.assembly_spades.todo == 1

  input:
    tuple val(sample), path(R1), path(R2)
  
  output:
    tuple val(sample), path("genomes/$sample/spades/contigs.fasta"), emit : draft_assembly
    path("genomes/$sample/spades/spades_1.log")
    path("genomes/$sample/spades/spades.err")

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
























