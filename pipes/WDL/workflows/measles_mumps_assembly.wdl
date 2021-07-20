version 1.0

import "tasks/tasks_trim.wdl" as mm_trim
import "tasks/tasks_fastqc.wdl" as mm_fastqc
import "tasks/tasks_bowtie2.wdl" as mm_bowtie2


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

  call mm_bowtie2.bowtie2_se {
  input:
    sra_id=sra_id,
    read1_trim=trim.read1_trim
    reference_seq=bowtie2_se.reference_seq
}

  output {
    File    read1_trim=trim.read1_trim
    File    fastqc_html=fastqc.fastqc_html
    File    sam_file=bowtie2_se.samfile
   #Float   gc_content=fastqc.gc_content

  }
}
