version 1.0

task sam_to_bam {

  input {
    String    sra_id
    File      samfile
  }

  command {
    date | tee DATE
    samtools --version | head -n1 | tee VERSION


    samtools view -S -b ${samfile}>${sra_id}.bam
    samtools sort ${sra_id}.bam -o ${sra_id}.sorted.bam
    samtools index ${sra_id}.sorted.bam
  }

  output {
    File    bamfile="${sra_id}.bam"
    File	sorted_bam="${sra_id}.sorted.bam"
    File	indexed_bam="${sra_id}.sorted.bam.bai"
    String     date          = read_string("DATE")
    String     version       = read_string("VERSION")
  }

  runtime {
    docker:       "staphb/samtools:1.12"
    memory:       "8 GB"
    cpu:          4
    disks:        "local-disk 100 SSD"
    preemptible:  1
  }
}
