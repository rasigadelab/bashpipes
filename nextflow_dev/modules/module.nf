process assembler_flye {
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
  
}
























