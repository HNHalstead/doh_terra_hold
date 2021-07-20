version 1.0

import "tasks/tasks_trim.wdl" as mm_trim
import "tasks/tasks_fastqc.wdl" as mm_fastqc


workflow mm_trim_and_assemble {

  input {
    String    sra_id
    File      read1
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

  output {
    File    read1_trim=trim.read1_trim
    File    fastqc_html=fastqc.fastqc_html
   #Float   gc_content=fastqc.gc_content

  }
}
