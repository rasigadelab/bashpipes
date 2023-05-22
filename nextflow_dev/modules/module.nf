process quality_fastqc {
  // Tool: fastqc
  // Quality control on FASTQ reads. 

  label 'lowCPU'
  storeDir (params.result)
  debug false
  tag "Fastqc on $sample"

  when:
    params.quality_fastqc.todo == 1

  input:
    tuple val(sample), path(R1), path(R2)
  
  output:
    tuple val(sample), path("genomes/$sample/$R1"), path("genomes/$sample/$R2"), emit : illumina_reads
    path("genomes/$sample/fastqc/*")

  script:
  """
  OUT_DIR=genomes/$sample/fastqc
  mkdir -p -m 777 \${OUT_DIR}

  fastqc $R1 $R2 -o \${OUT_DIR} -t $task.cpus 1> \${OUT_DIR}/fastqc.log 2> \${OUT_DIR}/fastqc.err 
  
  """

  stub:
  """
  OUT_DIR=genomes/$sample/fastqc
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/fastqc.err
  touch \${OUT_DIR}/fastqc.log
  """ 

}

process trim_trimmomatic {
  // Tool: trimmomatic
  // Trimming adapters and filtering reads of bad quality

  label 'trimmomatic'
  storeDir (params.result)
  debug false
  tag "Trimmomatic on $sample"  

  when:
    params.trim_trimmomatic.todo == 1

  input:
    tuple val(sample), path(R1), path(R2)
  
  output:
    tuple val(sample), path("genomes/$sample/trimmomatic/${sample}_R1_paired.trimmed.fastq.gz"), path("genomes/$sample/trimmomatic/${sample}_R2_paired.trimmed.fastq.gz"), path("genomes/$sample/trimmomatic/${sample}_R1_unpaired.trimmed.fastq.gz"), path("genomes/$sample/trimmomatic/${sample}_R2_unpaired.trimmed.fastq.gz"), emit : illumina_trimmed
    path("genomes/$sample/trimmomatic/*")

  script:
  """
  OUT_DIR=genomes/$sample/trimmomatic
  TRIM_R1=\${OUT_DIR}/$sample"_R1_paired.trimmed.fastq.gz"
  TRIM_R1_unpaired=\${OUT_DIR}/$sample"_R1_unpaired.trimmed.fastq.gz"
  TRIM_R2=\${OUT_DIR}/$sample"_R2_paired.trimmed.fastq.gz"
  TRIM_R2_unpaired=\${OUT_DIR}/$sample"_R2_unpaired.trimmed.fastq.gz"
  PARAMS="ILLUMINACLIP:${params.trim_trimmomatic["adapter"]}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
  mkdir -p -m 777 \${OUT_DIR}

  java -jar ${params.trim_trimmomatic["Trimmomatic"]} PE -validatePairs -threads $task.cpus $R1 $R2 \$TRIM_R1 \$TRIM_R1_unpaired \$TRIM_R2 \$TRIM_R2_unpaired \$PARAMS \
    1> \${OUT_DIR}/trimmomatic.log 2> \${OUT_DIR}/trimmomatic.err
  
  """

  stub:
  """
  OUT_DIR=genomes/$sample/trimmomatic
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/$sample"_R1_paired.trimmed.fastq.gz"
  touch \${OUT_DIR}/$sample"_R1_unpaired.trimmed.fastq.gz"
  touch \${OUT_DIR}/$sample"_R2_paired.trimmed.fastq.gz"
  touch \${OUT_DIR}/$sample"_R2_unpaired.trimmed.fastq.gz"
  touch \${OUT_DIR}/trimmomatic.err
  touch \${OUT_DIR}/trimmomatic.log
  """     

}

process assembly_flye {

  // Tool: Flye
  // De novo assembler for Oxford Nanopore long reads. 

  label 'flye'
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
  """

  stub:
  """
  OUT_DIR=genomes/$sample/flye
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/assembly.fasta
  touch genomes/$sample/${sample}_assembly_raw.fasta
  touch \${OUT_DIR}/flye.log
  touch \${OUT_DIR}/flye.err
  touch \${OUT_DIR}/assembly_info.txt
  echo "Completed Flye process for sample $sample"
  """  
}

process assembly_spades {

  // Tool: SPAdes. 
  // De novo assembler for Illumina short reads.

  label 'spades'
  storeDir params.result
  debug false
  tag "SPAdes on $sample"

  when:
    params.assembly_spades.todo == 1

  input:
    tuple val(sample), path(R1), path(R2), path(R1_UNPAIRED), path(R2_UNPAIRED)
  //  tuple val(sample), path(R1), path(R2)
  
  output:
    tuple val(sample), path("genomes/$sample/spades/contigs.fasta"), emit : draft_assembly
    path("genomes/$sample/spades/spades_1.log")
    path("genomes/$sample/spades/spades.err")
    path("genomes/$sample/${sample}_assembly_raw.fasta")

  script:
  memory = (task.memory =~ /([^\ ]+)(.+)/)[0][1]
  """
  OUT_DIR=genomes/$sample/spades
  mkdir -p -m 777 \${OUT_DIR}

  #If unpaired files are empty
  UNPAIRED_R1_SIZE=\$(wc -c $R1_UNPAIRED | awk '{print \$1}')
  UNPAIRED_R2_SIZE=\$(wc -c $R2_UNPAIRED | awk '{print \$1}')
  if [ \${UNPAIRED_R1_SIZE} -le 1000 ] || [ \${UNPAIRED_R2_SIZE} -le 1000 ]; then
      spades.py -1 $R1 -2 $R2 ${params.assembly_spades["sample_type"]} -o \${OUT_DIR} -t $task.cpus -m ${memory} 1> \${OUT_DIR}/spades_1.log 2> \${OUT_DIR}/spades.err
  else
      spades.py -1 $R1 -2 $R2 --s1 $R1_UNPAIRED --s2 $R2_UNPAIRED ${params.assembly_spades["sample_type"]} -o \${OUT_DIR} -t $task.cpus -m ${memory} 1> \${OUT_DIR}/spades_1.log 2> \${OUT_DIR}/spades.err
  fi

  cp \${OUT_DIR}/contigs.fasta genomes/$sample/${sample}_assembly_raw.fasta

  echo "Completed SPAdes process for sample $sample"
  """
  /*
  If unpaired files are empty
  UNPAIRED_R1_SIZE=\$(wc -c $R1_UNPAIRED | awk '{print \$1}')
  UNPAIRED_R2_SIZE=\$(wc -c $R2_UNPAIRED | awk '{print \$1}')
  if [ \${UNPAIRED_R1_SIZE} -le 1000 ] || [ \${UNPAIRED_R2_SIZE} -le 1000 ]; then
      spades.py -1 $R1 -2 $R2 ${params.assembly_spades["sample_type"]} -o \${OUT_DIR} -t $task.cpus -m ${memory} 1> \${OUT_DIR}/spades_1.log 2> \${OUT_DIR}/spades.err
  else
      spades.py -1 $R1 -2 $R2 --s1 $R1_UNPAIRED --s2 $R2_UNPAIRED ${params.assembly_spades["sample_type"]} -o \${OUT_DIR} -t $task.cpus -m ${memory} 1> \${OUT_DIR}/spades_1.log 2> \${OUT_DIR}/spades.err
  fi
  */

  stub:
  """
  OUT_DIR=genomes/$sample/spades
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/contigs.fasta
  touch genomes/$sample/${sample}_assembly_raw.fasta
  touch \${OUT_DIR}/spades_1.log
  touch \${OUT_DIR}/spades.err
  echo "Completed SPAdes process for sample $sample"
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
  tag "Mapping on $sample"

  when:
    params.map_bowtie2.todo == 1

  input:
    tuple val(sample), path(draft_assembly), path(R1), path(R2), path(R1_UNPAIRED), path(R2_UNPAIRED)

  output:
    tuple val(sample), path("genomes/$sample/polish/${sample}.sorted.bam"), emit : sorted_bam_files
    path("genomes/$sample/polish/${sample}.sorted.bam.bai")
    path("genomes/$sample/polish/*.log")


  script:
  memory = (task.memory =~ /([^\ ]+)(.+)/)[0][1]
  """
  OUT_DIR=genomes/$sample/polish
  SAMPLE_BAM=\${OUT_DIR}/${sample}.bam
  SAMPLE_BAM_SORTED=\${OUT_DIR}/${sample}.sorted.bam
  mkdir -p -m 777 \${OUT_DIR}

  bowtie2-build $draft_assembly \${OUT_DIR}/index &> \${OUT_DIR}/bowtie2.index.log
  bowtie2 -x \${OUT_DIR}/index -1 $R1 -2 $R2 -p $task.cpus 2>> \${OUT_DIR}/bowtie2.map.log | samtools view -bS > \${SAMPLE_BAM}
  
  samtools sort \${SAMPLE_BAM} -o \${SAMPLE_BAM_SORTED} -@ $task.cpus -m ${memory}G 2>> \${OUT_DIR}/samtools.sort.log
  samtools index \${SAMPLE_BAM_SORTED} -@ $task.cpus 2>> \${OUT_DIR}/samtools.index.log
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
  """ 
}

process polish_pilon {
  // Tool: pilon. 
  // Polishing step consists in correcting errors in draft assembly with Illumina reads alignment.

  label 'highCPU'
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
    path("genomes/$sample/polish/${sample}_polished.changes")

  script:
  memory = (task.memory =~ /([^\ ]+)(.+)/)[0][1]
  """
  OUT_DIR=genomes/$sample/polish
  SAMPLE_POLISHED=\${OUT_DIR}/${sample}_polished
  mkdir -p -m 777 \${OUT_DIR}

  pilon -Xmx${memory}G --genome $draft_assembly --bam $sorted_bam ${params.polish_pilon["list_changes"]} --output \${SAMPLE_POLISHED} &> \${OUT_DIR}/pilon.log
  cp \${SAMPLE_POLISHED}.fasta genomes/$sample/${sample}_polished.fasta

  echo "Completed Pilon process for sample $sample"
  """

  stub:
  """
  OUT_DIR=genomes/$sample/polish
  SAMPLE_POLISHED=\${OUT_DIR}/${sample}_polished
  mkdir -p -m 777 \${OUT_DIR}
  touch \${SAMPLE_POLISHED}.fasta
  touch genomes/$sample/${sample}_polished.fasta
  touch \${OUT_DIR}/pilon.log

  echo "Completed Pilon process for sample $sample"
  """

}

process qc_quast {
  // Tool: quast. 
  // Assembly QC step consists in checking quality of de novo Illumina assembly.

  label 'highCPU'
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
  """

  stub:
  """
  OUT_DIR=genomes/$sample/quast
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/quast.log
  touch \${OUT_DIR}/quast.err
  """  

}

process fixstart_circlator {
  // Tool: circlator. 
  // Circularization step consists in re-aligning contig sequences to origin of replication.

  label 'lowCPU'
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
  """

  stub:
  """
  OUT_DIR=genomes/$sample/circlator
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/${sample}_realigned.fasta
  touch genomes/$sample/${sample}_realigned.fasta
  touch \${OUT_DIR}/circlator.log
  """  
}

process mlst_sequence_typing {
  // Tool: mlst. 
  // Sequence typing consists in writing ST annotation for each isolate by checking the allele version of 7 house-keeping genes.

  label 'highCPU'
  storeDir params.result
  debug false
  tag "MLST on $sample"

  when:
    params.mlst_sequence_typing.todo == 1

  input:
    tuple val(sample), path(final_assembly)

  output:
    tuple val(sample), path("genomes/$sample/$final_assembly"), emit : final_assembly
    path("genomes/$sample/mlst/*")

  script:
  """
  OUT_DIR=genomes/$sample/mlst
  mkdir -p -m 777 \${OUT_DIR}

  mlst $final_assembly --threads $task.cpus > \${OUT_DIR}/mlst.tsv 2>> \${OUT_DIR}/mlst.log
  """

  stub:
  """
  OUT_DIR=genomes/$sample/mlst
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/mlst.tsv
  touch \${OUT_DIR}/mlst.log
  """  
 
}

process classify_sourmash {
  // Tool: sourmash. 
  // Taxon classification consists in describing right taxonomy for each isolate. 
  // Outputs specifically genus and species in this pipeline.

  label 'highCPU'
  storeDir params.result
  debug false
  tag "Sourmash on $sample"

  when:
    params.mlst_sequence_typing.todo == 1

  input:
    tuple val(sample), path(final_assembly)

  output:
    tuple val(sample), path("genomes/$sample/sourmash/sourmash.csv"), emit : sample_taxonomy
    path("genomes/$sample/sourmash/*")

  script:
  """
  OUT_DIR=genomes/$sample/sourmash
  OUT_SIG_FILE=\${OUT_DIR}/${sample}_realigned.fasta.sig
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/sourmash.csv
  sourmash sketch dna -o \${OUT_SIG_FILE} -p scaled=${params.classify_sourmash["scale"]},k=${params.classify_sourmash["k"]} $final_assembly &> \${OUT_DIR}/sourmash.log
  sourmash lca classify --query \${OUT_SIG_FILE} --db ${params.classify_sourmash["db"]} > \${OUT_DIR}/sourmash.csv 2>> \${OUT_DIR}/sourmash.log
  """

  stub:
  """
  OUT_DIR=genomes/$sample/sourmash
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/sourmash.csv
  echo "genus,species" > \${OUT_DIR}/sourmash.csv
  echo "Klebsiella,pneumoniae" >> \${OUT_DIR}/sourmash.csv
  touch \${OUT_DIR}/sourmash.log
  """  

}

process amr_typer_amrfinder {
  // Tool: AMRFinder. 
  // AMR genes typing consists in listing AMR genes present in the genome + some other genes of interest like biocide, stress response or virulence genes. 
  // Behave differently depending on genus and species given by Sourmash.

  label 'highCPU'
  storeDir params.result
  debug false
  tag "AMRFinder on $sample"

  when:
    params.amr_typer_amrfinder.todo == 1

  input:
    tuple val(sample), path(final_assembly), path(taxonomy_file)

  output:
    tuple val(sample), path("genomes/$sample/$final_assembly"), path("genomes/$sample/sourmash/$taxonomy_file"), emit : final_assembly
    path("genomes/$sample/amrfinder/*")
  
  script:
  """
  OUT_DIR=genomes/$sample/amrfinder
  mkdir -p -m 777 \${OUT_DIR}

  # Extract genus and species names if available
  GENUS=\$(cut -d',' -f8 $taxonomy_file | tail -n 1)
  SPECIES=\$(cut -d',' -f9 $taxonomy_file | tail -n 1)
  # Set Organism flag if genus/species are specified
  ORGANISM_FLAG=""
  case "\$GENUS" in
      "Escherichia" ) ORGANISM_FLAG="--organism Escherichia";;
      *) ORGANISM_FLAG=""
  esac

  amrfinder --plus --nucleotide $final_assembly --threads $task.cpus \${ORGANISM_FLAG} > \${OUT_DIR}/amrfinder.tsv 2> \${OUT_DIR}/amrfinder.log
  """

  stub:
  """
  OUT_DIR=genomes/$sample/amrfinder
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/amrfinder.tsv
  touch \${OUT_DIR}/amrfinder.log
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
    path("genomes/$sample/prokka/*")

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

  stub:
  """
  OUT_DIR=genomes/$sample/prokka
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/prokka.gff
  touch \${OUT_DIR}/prokka.log
  """  
}

process mge_mob_recon {
  // Tool: Mob_recon. 
  // MGE Analysis consists in finding and annotating plasmids, transposons and other mobile genetic elements encountered in the genome. 

  label 'highCPU'
  storeDir params.result
  debug true
  tag "Mob_Recon on $sample" 

  when:
    params.mge_mob_recon.todo == 1

  input:
    tuple val(sample), path(final_assembly)

  output:
    path("genomes/$sample/mob_recon/*")

  script:
  """
  OUT_DIR=genomes/$sample/mob_recon
  mkdir -p -m 777 \${OUT_DIR}
  
  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate mob_suite2
  mob_recon -n $task.cpus --force --infile $final_assembly --outdir \${OUT_DIR} 1> \${OUT_DIR}/mob_recon.log 2> \${OUT_DIR}/mob_recon.err
  conda deactivate
  """

  stub:
  """
  OUT_DIR=genomes/$sample/mob_recon
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/mob_recon.err
  touch \${OUT_DIR}/mob_recon.log
  """  

}









