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

process resync_bbmap {
  // Tool: bbmap repair.sh
  // Resynchronizing R1/R2 reads inside fastq files

  label 'resync'
  storeDir (params.result)
  debug false
  tag "R1/R2 resync on $sample"

  when:
    params.resync_bbmap.todo == 1

  input:
    tuple val(sample), path(R1), path(R2), path(R1_UNPAIRED), path(R2_UNPAIRED)
  
  output:
    tuple val(sample), path("genomes/$sample/trimmomatic/${sample}_R1.trimmed.resync.fastq.gz"), path("genomes/$sample/${sample}_R2.trimmed.resync.fastq.gz"), path("genomes/$sample/trimmomatic/${sample}_R1_unpaired.trimmed.fastq.gz"), path("genomes/$sample/trimmomatic/${sample}_R2_unpaired.trimmed.fastq.gz"), emit : illumina_resync
    path("genomes/$sample/fastp/fastp.json")
    path("genomes/$sample/fastp/fastp.html")
    path("genomes/$sample/trimmomatic/resync.err")
    path("genomes/$sample/trimmomatic/resync.log")

  script:
  """
  OUT_DIR=genomes/$sample/trimmomatic
  mkdir -p -m 777 \${OUT_DIR}

  repair.sh in1=$R1 in2=$R2 out1=$\{OUT_DIR}/${sample}_R1.trimmed.resync.fastq.gz \
   out2=\${OUT_DIR}/${sample}_R2.trimmed.resync.fastq.gz outs=\${OUT_DIR}/singletons.fastq.gz repair ow=t \
   1> \${OUT_DIR}/resync.log 2> \${OUT_DIR}/resync.err
  """
}

process assembly_spades {

  // Tool: SPAdes. 
  // De novo assembler for Illumina short reads.
  // Filter for only contigs > 500bp

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
    path("genomes/$sample/spades/assembly_graph_with_scaffolds.gfa")
    path("genomes/$sample/spades/contigs.paths")

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

  echo "Completed SPAdes process for sample $sample"
  """
}

process filter_contigs_bbmap {

  // Tool: BBmap reformat.sh. 
  // Filter for only contigs > 500bp in final assembly

  label 'bbmap'
  storeDir params.result
  debug false
  tag "Contigs filtering on $sample"

  when:
    params.filter_contigs_bbmap.todo == 1

  input:
    tuple val(sample), path(raw_assembly)
  
  output:
    tuple val(sample), path("genomes/$sample/spades/contigs_filtered.fasta"), emit : draft_assembly
    path("genomes/$sample/spades/bbmap.log")
    path("genomes/$sample/spades/bbmap.err")
    path("genomes/$sample/${sample}_assembly_raw.fasta")

  script:
  memory = (task.memory =~ /([^\ ]+)(.+)/)[0][1]
  """
  OUT_DIR=genomes/$sample/spades
  mkdir -p -m 777 \${OUT_DIR}

  reformat.sh in=$raw_assembly out=\${OUT_DIR}/contigs_filtered.fasta minlength=500 1> \${OUT_DIR}/bbmap.log 2> \${OUT_DIR}/bbmap.err
  cp \${OUT_DIR}/contigs_filtered.fasta genomes/$sample/${sample}_assembly_raw.fasta

  """
}

process qc_quast {
  // Tool: quast. 
  // Assembly QC step consists in checking quality of de novo Illumina assembly.

  label 'quast'
  storeDir params.result
  debug false
  tag "Quast QC on $sample"

  when:
    params.qc_quast.todo == 1

  input:
    tuple val(sample), path(draft_assembly)

  output:
    tuple val(sample), path("genomes/$sample/spades/$draft_assembly"), emit : draft_assembly
    path("genomes/$sample/quast/basic_stats/*")
    path("genomes/$sample/quast/icarus_viewers/*")
    path("genomes/$sample/quast/icarus.html")
    path("genomes/$sample/quast/quast.err")
    path("genomes/$sample/quast/quast.log")
    path("genomes/$sample/quast/report.html")
    path("genomes/$sample/quast/report.pdf")
    path("genomes/$sample/quast/report.tsv")
    path("genomes/$sample/quast/transposed_report.tsv")
    path("genomes/$sample/quast/transposed_report.txt")

  script:
  """
  OUT_DIR=genomes/$sample/quast
  mkdir -p -m 777 \${OUT_DIR}

  quast.py -o \${OUT_DIR} $draft_assembly 1> \${OUT_DIR}/quast.log 2> \${OUT_DIR}/quast.err
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
  TMP_DIR=\${OUT_DIR}/tmp
  mkdir -p -m 777 \${TMP_DIR}
  
  export TMPDIR=\${TMP_DIR}
  bakta --force --prefix $sample --threads $task.cpus --output \${OUT_DIR} --skip-trna --keep-contig-headers --skip-crispr --db ${params.annotate_bakta["db"]} --tmp-dir \$TMP_DIR $final_assembly 1> \${OUT_DIR}/bakta.log 2> \${OUT_DIR}/bakta.err
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








