/*
 * NOTE FOR LATER
 * Command to run the script
 * nextflow run bactopia_v0.nf -c bactopia.config 
 * One can add -resume option 
 *
 * To change process parameters => bactopia.config file
 * (CPUs, RAM)
*/

params.datadir = '/mnt/c/Users/admin/Documents/Projets/5.pipeline_maintenance'
params.ont = "$params.datadir/genomes/*/*_ONT.fastq.gz"
//params.reads = "$projectDir/genomes/*/*_{R1,R2}.fastq.gz"

samples_dir = file(params.ont).parent

process create_files_location {
  input:
    path x

  output:
    path "5.pipeline_maintenance/Files_location.tsv"

  """
  python3 /mnt/c/Users/admin/Documents/Projets/4.Pipeline_improvements/1.Files_location_script/files_location.py -d $x 
  """
}

process merger_preparation {
  //A TESTER POUR AFFICHER SUDO MDP debug true

  input:
    path x
  
  output:
    stdout

  """
  sudo Rscript /mnt/c/Users/admin/Documents/Projets/4.Pipeline_improvements/5.Prepare_reads/230104_nextflow_prepare_reads.R -d $x
  """
}

process draft_assembly {
  tag "FLYE on $reads.simpleName"
  publishDir "$sample_path", mode: 'copy'
  debug true

  input:
    path(reads)
    val sample_path
  
  output:
    path("flye/assembly.fasta")
    path("flye/flye.log")
    path("flye/flye.err")
    path("flye/assembly_info.txt")
  
  """
  SAMPLE_DIR=$sample_path
  SAMPLE="\${SAMPLE_DIR##*/}"
  echo "Completed Flye process for sample \$SAMPLE"
  mkdir -p "flye"
  flye --nano-raw $reads -o ./flye --threads $task.cpus 1> flye/flye.log 2> flye/flye.err 
  cp ./flye/assembly.fasta $sample_path/\${SAMPLE}_assembly_raw.fasta
  """
  //flye --nano-raw $reads -o "$sample_path/flye" --threads $task.cpus 1> flye.log 2> flye.err && 
  //cp $sample_path/flye/assembly.fasta $sample_path/\${SAMPLE}_assembly_raw.fasta
  
}


workflow {
  //create_files_location(params.datadir)
  //merger_preparation(params.datadir)
  ont_ch = Channel.fromPath(params.ont)
  samples_ch = Channel.fromList(samples_dir)
  draft_assembly(ont_ch, samples_ch)
  //illumina_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
 
}