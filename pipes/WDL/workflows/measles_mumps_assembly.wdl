version 1.0

import "https://github.com/HNHalstead/doh_terra_hold/blob/dev/pipes/WDL/tasks/tasks_trim.wdl" as trim

workflow mm_trim_and_assemble {
  meta {
        description: "Assembles single end Illumina reads for Measles and Mumps."
        author: "Holly Halstead"
        email:  "holly.halstead@doh.wa.gov"
    }

  input {
    String    sra_id
    File      read1
  }

  call trim {
    input:
      sra_id=sra_id,
      read1=read1
  }

    call fastqc {
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


task fastqc {

  input {
    String    sra_id
    File      read1_trim
  }

  command {
	set -euo pipefail
	fastqc ${read1_trim} -o .
  }

  output {
    File    fastqc_html=glob("*fastqc.html")[0]
    File	fastqc_zip=glob("*fastqc.zip")[0]
  }

  runtime {
    docker:       "staphb/fastqc:0.11.8"
    memory:       "8 GB"
    cpu:          4
    disks:        "local-disk 100 SSD"
    preemptible:  1
  }
}
