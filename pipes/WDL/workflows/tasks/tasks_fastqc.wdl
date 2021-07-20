version 1.0

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
