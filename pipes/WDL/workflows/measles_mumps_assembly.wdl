version 1.0

import "tasks/tasks_trim.wdl" as mm_trim
import "tasks/tasks_fastqc.wdl" as mm_fastqc
#import "tasks/tasks_bowtie2.wdl" as mm_bowtie2

workflow mm_trim_and_assemble {

  input {
    String    sra_id
    File      read1
    File      reference_seq
  }

  call mm_trim.trim {
    input:
      sra_id=sra_id,
      read1=read1
  }

    call mm_fastqc.fastqc {
      input:
        sra_id=sra_id,
        read1_trim=trim.read1_trim
  }

  call bowtie2_se {
  input:
    sra_id=sra_id,
    read1_trim=trim.read1_trim,
    reference_seq=reference_seq
  }

  output {
    File    read1_trim=trim.read1_trim
    File    fastqc_html=fastqc.fastqc_html
    File    sam_file=bowtie2_se.samfile
   #Float   gc_content=fastqc.gc_content

  }
}

task bowtie2_se {
  input {
    String  sra_id
    File    read1_trim
    File    reference_seq
  }

  command {
  bowtie2-build ${reference_seq} mumps_ref
  bowtie2 -x mumps_ref -U ${read1_trim} -S ${sra_id}.sam --local
  }

  output {
    File    samfile="${sra_id}.sam"
  }

  runtime {
    docker:       "quay.io/biocontainers/bowtie2:2.4.4--py38h72fc82f_0"
    memory:       "8 GB"
    cpu:          4
    disks:        "local-disk 100 SSD"
    preemptible:  1
  }
}
