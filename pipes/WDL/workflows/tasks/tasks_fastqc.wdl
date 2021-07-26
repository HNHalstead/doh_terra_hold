version 1.0

task fastqc {

  input {
    String    sra_id
    File      read1_trim
    File      read2_trim
  }

  command {
	set -euo pipefail
	fastqc ${read1_trim} -o .
  fastqc ${read2_trim} -o .
  }

  output {
    File  fastqc_html_r1=glob("*fastqc.html")[0]
    File	fastqc_zip_r1=glob("*fastqc.zip")[0]
    File  fastqc_html_r2=glob("*fastqc.html")[1]
    File	fastqc_zip_r2=glob("*fastqc.zip")[1]
  }

  runtime {
    docker:       "staphb/fastqc:0.11.8"
    memory:       "8 GB"
    cpu:          4
    disks:        "local-disk 100 SSD"
    preemptible:  1
  }
}
