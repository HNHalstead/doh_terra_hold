version 1.0

task bowtie2_se {
  input {
    String  sra_id
    File    read1_trim
    File    read2_trim
    File    reference_seq
  }

  command {
    date | tee DATE
    bowtie2 --version | head -n1 | tee VERSION

    bowtie2-build ${reference_seq} mumps_ref
    bowtie2 -x mumps_ref -U ${read1_trim},${read2_trim} -S ${sra_id}.sam --local
  }

  output {
    File    samfile="${sra_id}.sam"
    String     date          = read_string("DATE")
    String     version       = read_string("VERSION")
  }

  runtime {
    docker:       "quay.io/biocontainers/bowtie2:2.4.4--py38h72fc82f_0"
    memory:       "8 GB"
    cpu:          4
    disks:        "local-disk 100 SSD"
    preemptible:  1
  }
}

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
