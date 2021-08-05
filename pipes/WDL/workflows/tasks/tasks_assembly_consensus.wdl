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

task consensus {

  input {
    String    sra_id
    File      sorted_bam
    File      reference_seq
    Boolean     count_orphans = true
    Int         max_depth = "600000"
    Boolean     disable_baq = true
    Int         min_bq = "0"
    Int         min_qual = "20"
    Float       min_freq = "0.6"
    Int         min_depth = "10"
  }


  command {
    date | tee DATE
    samtools --version | head -n1 | tee VERSION

    bcftools mpileup -Ou -f ${reference_seq} -b ${bamfile} | bcftools call -mv Ob -o | vcfutils.pl vcf2fq > ${sra_id}_cns.fastq

    samtools mpileup -d 60000 -Q 0 --reference ${reference_seq} ${sorted_bam} | ivar consensus -p SRR12100584.consensus -q 20 -t 0.6 -m 10 -n N

  }

  output {
    File    samfile="${sra_id}.sam"
    String     date          = read_string("DATE")
    String     version       = read_string("VERSION")
  }

  runtime {
    docker:       "staphb/bcftools:1.12"
    memory:       "8 GB"
    cpu:          4
    disks:        "local-disk 100 SSD"
    preemptible:  1
  }
  }
