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
