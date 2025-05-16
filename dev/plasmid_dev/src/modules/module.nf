process gather_plasmid_seq {
  // Tool: bash
  // Gathering of plasmid sequences. 

  label 'gather'
  storeDir (params.result)
  debug false
  tag "Gather plasmid seq of $sample for resistance $resistance"

  when:
    params.gather_plasmid_seq.todo == 1

  input:
    tuple val(sample), val(resistance), val(plasmid_id), val(inc), path(seq)
  
  output:
    tuple val(resistance), val(inc), val(sample), emit : plasmid_seq
    path("analyses/$resistance/$inc/plasmids/${sample}.fasta")

  script:
  """
  OUT_DIR=analyses/$resistance/$inc/plasmids
  mkdir -p -m 777 \${OUT_DIR}

  cat ${params.result}/genomes/$sample/mob_recon/$seq > \${OUT_DIR}/${sample}.fasta
  """
}

process align_to_reference {
  // Tool: quast
  // Align plasmid sequences to plasmid reference

  label 'map'
  storeDir (params.result)
  debug false
  tag "Align to reference for $inc/$resistance plasmids"  

  when:
    params.align_to_reference.todo == 1

  input:
    tuple val(resistance), val(inc), val(sample), val(batch)
  
  output:
    tuple val(resistance), val(inc), val(sample), emit : plasmid_seq
    path("analyses/$resistance/$inc/quast_results_$batch/*")


  script:
  """
  TIMESTAMP=\$(date +%Y%m%d_%H%M%S)
  OUT_DIR=analyses/$resistance/$inc/quast_results_$batch
  SEQ_DIR=${params.result}/analyses/$resistance/$inc/plasmids
  mkdir -p -m 777 \${OUT_DIR}

  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate bactopia
  # Command is same but reference adapts to resistance and plasmid inc type
  REFERENCE="${params.nfpath}/ref_db"
  SAMPLES=\$(echo $sample | tr -d '[],')
  SAMPLES_FASTA=\$(for s in \$SAMPLES; do echo -n "\${SEQ_DIR}/\${s}.fasta "; done)
  if [ $resistance = "blaOXA-48" -a $inc = "IncL-M" ]
  then
    REFERENCE=\$REFERENCE/${params.align_to_reference["reference_oxa48_incLM"]}
    quast -r \$REFERENCE.fasta -g \$REFERENCE.gff3 --circos -o \${OUT_DIR} \$SAMPLES_FASTA
  elif [ $resistance = "blaVIM-1" -a $inc = "IncN" ]
  then
    REFERENCE=\$REFERENCE/${params.align_to_reference["reference_vim1_incN"]}
    quast -r \$REFERENCE.fasta -g \$REFERENCE.gff3 --circos -o \${OUT_DIR} \$SAMPLES_FASTA
  elif [ $resistance = "blaVIM-1" -a $inc = "IncL-M" ]
  then
    REFERENCE=\$REFERENCE/${params.align_to_reference["reference_vim1_incLM"]}
    quast -r \$REFERENCE.fasta -g \$REFERENCE.gff3 --circos -o \${OUT_DIR} \$SAMPLES_FASTA
  elif [ $resistance = "blaVIM-4" -a $inc = "IncHI2A,rep_cluster_1088" ]
  then
    REFERENCE=\$REFERENCE/${params.align_to_reference["reference_vim4_incHI2A"]}
    quast -r \$REFERENCE.fasta -g \$REFERENCE.gff3 --circos -o \${OUT_DIR} \$SAMPLES_FASTA
  elif [ $resistance = "blaOXA-48" -a $inc = "IncC" ]
  then
    REFERENCE=\$REFERENCE/${params.align_to_reference["reference_oxa48_incC"]}
    quast -r \$REFERENCE.fasta -g \$REFERENCE.gff3 --circos -o \${OUT_DIR} \$SAMPLES_FASTA
  elif [ $resistance = "blaOXA-48" -a $inc = "IncL-M,IncR" ]
  then
    REFERENCE=\$REFERENCE/${params.align_to_reference["reference_oxa48_incLMR"]}
    quast -r \$REFERENCE.fasta -g \$REFERENCE.gff3 --circos -o \${OUT_DIR} \$SAMPLES_FASTA
  elif [ $resistance = "blaOXA-48" -a $inc = "-" ]
  then
    REFERENCE=\$REFERENCE/${params.align_to_reference["reference_oxa48_noInc"]}
    quast -r \$REFERENCE.fasta -g \$REFERENCE.gff3 --circos -o \${OUT_DIR} \$SAMPLES_FASTA
  elif [ $resistance = "blaNDM-1" -a $inc = "IncN" ]
  then
    REFERENCE=\$REFERENCE/${params.align_to_reference["reference_ndm1_incN"]}
    quast -r \$REFERENCE.fasta -g \$REFERENCE.gff3 --circos -o \${OUT_DIR} \$SAMPLES_FASTA
  elif [ $resistance = "blaNDM-1" -a $inc = "IncFIB" ]
  then
    REFERENCE=\$REFERENCE/${params.align_to_reference["reference_ndm1_incFIB"]}
    quast -r \$REFERENCE.fasta -g \$REFERENCE.gff3 --circos -o \${OUT_DIR} \$SAMPLES_FASTA
  elif [ $resistance = "blaNDM-1" -a $inc = "IncHI1A" ]
  then
    REFERENCE=\$REFERENCE/${params.align_to_reference["reference_ndm1_incHI1A"]}
    quast -r \$REFERENCE.fasta -g \$REFERENCE.gff3 --circos -o \${OUT_DIR} \$SAMPLES_FASTA
  elif [ $resistance = "blaNDM-1" -a $inc = "IncC" ]
  then
    REFERENCE=\$REFERENCE/${params.align_to_reference["reference_ndm1_incC"]}
    quast -r \$REFERENCE.fasta -g \$REFERENCE.gff3 --circos -o \${OUT_DIR} \$SAMPLES_FASTA
  elif [ $resistance = "blaNDM-1" -a $inc = "-" ]
  then
    REFERENCE=\$REFERENCE/${params.align_to_reference["reference_ndm1_noInc"]}
    quast -r \$REFERENCE.fasta -g \$REFERENCE.gff3 --circos -o \${OUT_DIR} \$SAMPLES_FASTA
  elif [ $resistance = "blaNDM-1" -a $inc = "IncC,IncN,IncR" ]
  then
    REFERENCE=\$REFERENCE/${params.align_to_reference["reference_ndm1_incCNR"]}
    quast -r \$REFERENCE.fasta -g \$REFERENCE.gff3 --circos -o \${OUT_DIR} \$SAMPLES_FASTA
  else
    echo "No comparison requested" > \${OUT_DIR}/here.txt
  fi
  conda deactivate
  """
}

process visualize_circos {
  // Tool: Circos
  // Refine Circos maps created by Quast at previous step

  label 'circos'
  storeDir params.result
  debug false
  tag "Visual of $resistance/$inc"

  when:
    params.visualize_circos.todo == 1

  input:
    tuple val(resistance), val(inc), val(sample)
  
  output:
    tuple val(resistance), val(inc), val(sample), emit : plasmids_info
    path("analyses/$resistance/$inc/visual/*")
    
  script:
  """
  OUT_DIR=analyses/$resistance/$inc/visual
  RESULTS_DIR=\$(ls $params.result/analyses/$resistance/$inc)
  QUAST=\$(echo "\$RESULTS_DIR" | grep "^quast")

  for quast_folder in \$QUAST; do
    mkdir -p -m 777 \${OUT_DIR}
    echo "Processing \$quast_folder" >> \${OUT_DIR}/out.txt
    CIRCOS_FOLDER=$params.result/analyses/$resistance/$inc/\${quast_folder}/circos
    if [ -d \$CIRCOS_FOLDER ]; 
    then
      mkdir -p -m 777 \${OUT_DIR}/\${quast_folder}/
      cp -r \$CIRCOS_FOLDER "\${OUT_DIR}/\${quast_folder}/"
    else
      echo "Circos folder does not exist. Check if a comparison has been requested." >> \${OUT_DIR}/out.txt
    fi
  done

  """
}

process change_circos_config {

  // Tool: Python script

  label 'circos_config'
  storeDir params.result
  debug false
  tag "Change circos config $resistance/$inc"

  when:
    params.change_circos_config.todo == 1

  input:
    tuple val(resistance), val(inc), val(sample)
  
  output:
    path("analyses/$resistance/$inc/visual/changes.txt")
    path("analyses/$resistance/$inc/visual/*")


  script:
  """
  OUT_DIR=analyses/$resistance/$inc/visual
  mkdir -p -m 777 \${OUT_DIR}
  RESULTS_DIR=\$(ls $params.result/analyses/$resistance/$inc)
  QUAST=\$(echo "\$RESULTS_DIR" | grep "^quast")

  echo "Changing configuration files." >> \${OUT_DIR}/changes.txt
  for quast_folder in \$QUAST; do
    CIRCOS_FOLDER=$params.result/analyses/$resistance/$inc/visual/\${quast_folder}/circos
    if [ -d \$CIRCOS_FOLDER ]; 
    then
      echo "Processing results in \$quast_folder" >> \${OUT_DIR}/changes.txt
      python3 $params.nfpath/modules/circos_plot.py -d \$CIRCOS_FOLDER
      source ~/miniconda3/etc/profile.d/conda.sh
      conda activate circos
      circos -conf \${CIRCOS_FOLDER}/circos.conf
      conda deactivate
    else
      echo "No results to process for \$quast_folder" >> \${OUT_DIR}/changes.txt
    fi
  done
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








