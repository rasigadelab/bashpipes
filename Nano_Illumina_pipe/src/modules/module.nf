process quality_fastp {
  // Tool: fastp
  // Quality control on FASTQ reads. 

  label 'fastp'
  storeDir (params.result)
  debug false
  tag "Fastp on $sample"

  when:
    params.quality_fastp.todo == 1

  input:
    tuple val(sample), path(R1), path(R2)
  
  output:
    tuple val(sample), path("genomes/$sample/$R1"), path("genomes/$sample/$R2"), emit : illumina_reads
    path("genomes/$sample/fastp/fastp.json")
    path("genomes/$sample/fastp/fastp.html")
    path("genomes/$sample/fastp/fastp.err")
    path("genomes/$sample/fastp/fastp.log")

  script:
  """
  OUT_DIR=genomes/$sample/fastp
  mkdir -p -m 777 \${OUT_DIR}

  fastp -i $R1 -I $R2 -w $task.cpus --json \${OUT_DIR}/fastp.json --html \${OUT_DIR}/fastp.html \
  1> \${OUT_DIR}/fastp.log 2> \${OUT_DIR}/fastp.err
  """
}

process trimming_porechop {
  // Tool: porechop
  // Trimming of ONT FASTQ reads. 

  label 'porechop'
  storeDir (params.result)
  debug false
  tag "Porechop on $sample"

  when:
    params.trimming_porechop.todo == 1

  input:
    tuple val(sample), path(ont_reads)
  
  output:
    tuple val(sample), path("genomes/$sample/porechop/${sample}_ONT_trimmed.fastq.gz"), emit : trimmed_ont_reads
    path("genomes/$sample/porechop/porechop.err")
    path("genomes/$sample/porechop/porechop.log")

  script:
  """
  OUT_DIR=genomes/$sample/porechop
  mkdir -p -m 777 \${OUT_DIR}

  porechop -i $ont_reads -o \${OUT_DIR}/${sample}_ONT_trimmed.fastq.gz 1> \${OUT_DIR}/porechop.log 2> \${OUT_DIR}/porechop.err
  """
}

process stats_nanoplot{
  
  // Computing statistics on Nanopore FastQ and plotting some plots related.

  label 'nanoplot'
  storeDir params.result
  debug false
  tag "Computing ONT stats of sample $sample"

  when:
    params.stats_nanoplot.todo == 1 

  input:
    tuple val(sample), path(ont_reads)

  output:
    tuple val(sample), path("genomes/$sample/filtlong/${sample}_filtered_ONT.fastq.gz"), emit : nanopore_reads
    path("genomes/$sample/nanoplot/*.log")
    path("genomes/$sample/nanoplot/NanoPlot-report.html")
    path("genomes/$sample/nanoplot/NanoStats.txt")
    
  
  """
  OUT_DIR=genomes/$sample/nanoplot
  mkdir -p -m 777 \${OUT_DIR}
  
  NanoPlot -t $task.cpus --fastq $ont_reads -o \${OUT_DIR} --tsv_stats --N50 -f png 1> \${OUT_DIR}/nanoplot.log 2> \${OUT_DIR}/nanoplot.err
  """
}

process trim_trimmomatic {
  // Tool: trimmomatic
  // Trimming adapters and filtering reads of bad quality

  label 'trim'
  storeDir (params.result)
  debug false
  tag "Trimmomatic on $sample"  

  when:
    params.trim_trimmomatic.todo == 1

  input:
    tuple val(sample), path(R1), path(R2)
  
  output:
    tuple val(sample), path("genomes/$sample/trimmomatic/${sample}_R1_paired.trimmed.fastq.gz"), path("genomes/$sample/trimmomatic/${sample}_R2_paired.trimmed.fastq.gz"), path("genomes/$sample/trimmomatic/${sample}_R1_unpaired.trimmed.fastq.gz"), path("genomes/$sample/trimmomatic/${sample}_R2_unpaired.trimmed.fastq.gz"), emit : illumina_trimmed
    path("genomes/$sample/trimmomatic/${sample}_R1_paired.trimmed.fastq.gz")
    path("genomes/$sample/trimmomatic/${sample}_R2_paired.trimmed.fastq.gz")
    path("genomes/$sample/trimmomatic/${sample}_R1_unpaired.trimmed.fastq.gz")
    path("genomes/$sample/trimmomatic/${sample}_R2_unpaired.trimmed.fastq.gz")
    path("genomes/$sample/trimmomatic/trimmomatic.err")
    path("genomes/$sample/trimmomatic/trimmomatic.log")

  script:
  """
  OUT_DIR=genomes/$sample/trimmomatic
  TRIM_R1=\${OUT_DIR}/${sample}_R1_paired.trimmed.fastq.gz
  TRIM_R1_unpaired=\${OUT_DIR}/${sample}_R1_unpaired.trimmed.fastq.gz
  TRIM_R2=\${OUT_DIR}/${sample}_R2_paired.trimmed.fastq.gz
  TRIM_R2_unpaired=\${OUT_DIR}/${sample}_R2_unpaired.trimmed.fastq.gz
  PARAMS="ILLUMINACLIP:${params.trim_trimmomatic["adapter"]}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
  mkdir -p -m 777 \${OUT_DIR}

  trimmomatic PE -validatePairs -threads $task.cpus $R1 $R2 \$TRIM_R1 \$TRIM_R1_unpaired \$TRIM_R2 \$TRIM_R2_unpaired \$PARAMS \
    1> \${OUT_DIR}/trimmomatic.log 2> \${OUT_DIR}/trimmomatic.err
  
  """
}

process filter_filtlong {
  // Tool: Filtlong
  // Filtering of Nanopore reads.  

  label 'filtlong'
  storeDir (params.result)
  debug false
  tag "Filtlong on $ont_reads.simpleName"

  when:
    params.filter_filtlong.todo == 1

  input:
    tuple val(sample), path(ont_reads)
  
  output:
    tuple val(sample), path("genomes/$sample/filtlong/${sample}_filtered_ONT.fastq.gz"), emit : filtered_nanopore
    path("genomes/$sample/filtlong/filtlong.err")

  script:
  """
  OUT_DIR=genomes/$sample/filtlong
  mkdir -p -m 777 \${OUT_DIR}

  filtlong --min_length ${params.filter_filtlong["min_read_length"]} --keep_percent ${params.filter_filtlong["keep_percent"]} $ont_reads | gzip > \${OUT_DIR}/${sample}_filtered_ONT.fastq.gz 2> \${OUT_DIR}/filtlong.err
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
    path("genomes/$sample/flye/*")
    // path("genomes/$sample/flye/flye.log")
    // path("genomes/$sample/flye/flye.err")
    // path("genomes/$sample/flye/assembly_info.txt")
    // path("genomes/$sample/flye/*.gfa")
    // path("genomes/$sample/${sample}_assembly_raw.fasta")

  script:
  """
  OUT_DIR=genomes/$sample/flye
  mkdir -p -m 777 \${OUT_DIR}

  flye ${params.assembly_flye["ont_type"]} $ont_reads -o \${OUT_DIR} --threads $task.cpus --meta 1> \${OUT_DIR}/flye.log 2> \${OUT_DIR}/flye.err 
  cp \${OUT_DIR}/assembly.fasta genomes/$sample/${sample}_assembly_raw.fasta

  echo "Completed Flye process for sample $sample"
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
    path("genomes/$sample/polish/bowtie2.index.log")
    path("genomes/$sample/polish/bowtie2.map.log")
    path("genomes/$sample/polish/samtools.index.log")
    path("genomes/$sample/polish/samtools.sort.log")

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
}

process polish_pilon {
  // Tool: pilon. 
  // Polishing step consists in correcting errors in draft assembly with Illumina reads alignment.

  label 'pilon'
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

  java -Xmx${memory}G -jar /pilon/pilon.jar --genome $draft_assembly --bam $sorted_bam ${params.polish_pilon["list_changes"]} --output \${SAMPLE_POLISHED} &> \${OUT_DIR}/pilon.log
  cp \${SAMPLE_POLISHED}.fasta genomes/$sample/${sample}_polished.fasta
  
  echo "Completed Pilon process for sample $sample"
  """
}

process fixstart_circlator {
  // Tool: circlator. 
  // Circularization step consists in re-aligning contig sequences to origin of replication.

  label 'circlator'
  storeDir params.result
  debug false
  tag "Circlator on $sample"

  when:
    params.fixstart_circlator.todo == 1

  input:
    tuple val(sample), path(denovo_assembly)

  output:
    tuple val(sample), path("genomes/$sample/circlator/${sample}_realigned.fasta"), emit : realigned_assembly
    path("genomes/$sample/circlator/${sample}_realigned.log")
    path("genomes/$sample/circlator/${sample}_realigned.detailed.log")
    path("genomes/$sample/${sample}_realigned.fasta")
 
  script:
  """
  OUT_DIR=genomes/$sample/circlator
  mkdir -p -m 777 \${OUT_DIR}

  circlator fixstart $denovo_assembly \${OUT_DIR}/${sample}_realigned
  cp \${OUT_DIR}/${sample}_realigned.fasta genomes/$sample/${sample}_realigned.fasta
  """
}

process mlst_sequence_typing {
  // Tool: mlst. 
  // Sequence typing consists in writing ST annotation for each isolate by checking the allele version of 7 house-keeping genes.

  label 'mlst'
  storeDir params.result
  debug false
  tag "MLST on $sample"

  when:
    params.mlst_sequence_typing.todo == 1

  input:
    tuple val(sample), path(final_assembly)

  output:
    tuple val(sample), path("genomes/$sample/$final_assembly"), emit : final_assembly
    path("genomes/$sample/mlst/mlst.log")
    path("genomes/$sample/mlst/mlst.tsv")

  script:
  """
  OUT_DIR=genomes/$sample/mlst
  mkdir -p -m 777 \${OUT_DIR}

  mlst $final_assembly --threads $task.cpus > \${OUT_DIR}/mlst.tsv 2>> \${OUT_DIR}/mlst.log
  """
}

process classify_sourmash {
  // Tool: sourmash. 
  // Taxon classification consists in describing right taxonomy for each isolate. 
  // Outputs specifically genus and species in this pipeline.

  label 'sourmash'
  storeDir params.result
  debug false
  tag "Sourmash on $sample"

  when:
    params.classify_sourmash.todo == 1

  input:
    tuple val(sample), path(final_assembly)

  output:
    tuple val(sample), path("genomes/$sample/sourmash/sourmash.csv"), emit : sample_taxonomy
    path("genomes/$sample/sourmash/${sample}_realigned.fasta.sig")
    path("genomes/$sample/sourmash/sourmash.csv")
    path("genomes/$sample/sourmash/sourmash.log")

  script:
  """
  OUT_DIR=genomes/$sample/sourmash
  OUT_SIG_FILE=\${OUT_DIR}/${sample}_realigned.fasta.sig
  mkdir -p -m 777 \${OUT_DIR}
  touch \${OUT_DIR}/sourmash.csv
  sourmash sketch dna -o \${OUT_SIG_FILE} -p scaled=${params.classify_sourmash["scale"]},k=${params.classify_sourmash["k"]} $final_assembly &> \${OUT_DIR}/sourmash.log
  sourmash lca classify --query \${OUT_SIG_FILE} --db ${params.classify_sourmash["db"]} > \${OUT_DIR}/sourmash.csv 2>> \${OUT_DIR}/sourmash.log
  """
}

process amr_typer_amrfinder {
  // Tool: AMRFinder. 
  // AMR genes typing consists in listing AMR genes present in the genome + some other genes of interest like biocide, stress response or virulence genes. 
  // Behave differently depending on genus and species given by Sourmash.

  label 'amrfinder'
  storeDir params.result
  debug false
  tag "AMRFinder on $sample"

  when:
    params.amr_typer_amrfinder.todo == 1

  input:
    tuple val(sample), path(final_assembly), path(taxonomy_file)

  output:
    tuple val(sample), path("genomes/$sample/$final_assembly"), path("genomes/$sample/sourmash/$taxonomy_file"), emit : final_assembly
    path("genomes/$sample/amrfinder/amrfinder.log")
    path("genomes/$sample/amrfinder/amrfinder.tsv")
  
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
}

process annotate_bakta {
  // Tool: bakta. 
  // Genome annotation consists in locating genes on the genome and giving their function 

  label 'bakta'
  storeDir params.result
  debug false
  tag "Bakta on $sample"  

  when:
    params.annotate_bakta.todo == 1

  input:
    tuple val(sample), path(final_assembly), path(taxonomy_file)

  output:
    tuple val(sample), path("genomes/$sample/$final_assembly"), emit : final_assembly
    path("genomes/$sample/bakta/bakta.err")
    path("genomes/$sample/bakta/${sample}.gff3")
    path("genomes/$sample/bakta/${sample}.tsv")
    path("genomes/$sample/bakta/${sample}.png")
    path("genomes/$sample/bakta/${sample}.log")
    path("genomes/$sample/bakta/bakta.log")

  script:
  """
  OUT_DIR=genomes/$sample/bakta
  mkdir -p -m 777 \${OUT_DIR}

  # Extract genus and species names if available
  GENUS=\$(cut -d',' -f8 $taxonomy_file | tail -n 1)
  SPECIES=\$(cut -d',' -f9 $taxonomy_file | tail -n 1)
  
  bakta --force --prefix $sample --threads $task.cpus --output \${OUT_DIR} --keep-contig-headers --db ${params.annotate_bakta["db"]} $final_assembly 1> \${OUT_DIR}/bakta.log 2> \${OUT_DIR}/bakta.err
  """
}

process mge_mob_recon {
  // Tool: Mob_recon. 
  // MGE Analysis consists in finding and annotating plasmids, transposons and other mobile genetic elements encountered in the genome. 

  label 'mob_recon'
  storeDir params.result
  debug false
  tag "Mob_Recon on $sample" 

  when:
    params.mge_mob_recon.todo == 1

  input:
    tuple val(sample), path(final_assembly)

  output:
    tuple val(sample), path("genomes/$sample/mob_recon/contig_report.txt"), emit : samples_list
    path("genomes/$sample/mob_recon/*.fasta")
    path("genomes/$sample/mob_recon/mge.report.txt")
    path("genomes/$sample/mob_recon/*.txt")
    path("genomes/$sample/mob_recon/mob_recon.log")
    path("genomes/$sample/mob_recon/mob_recon.err")

  script:
  """
  OUT_DIR=genomes/$sample/mob_recon
  mkdir -p -m 777 \${OUT_DIR}
  
  mob_recon -n $task.cpus --force --infile $final_assembly --outdir \${OUT_DIR} 1> mob_recon.log 2> mob_recon.err
  mv mob_recon.log \${OUT_DIR}
  mv mob_recon.err \${OUT_DIR}
  """
}








