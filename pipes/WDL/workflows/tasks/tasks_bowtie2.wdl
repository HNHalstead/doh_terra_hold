version 1.0



task bowtie2_se {
  input {
    String  sra_id
    File    read1_trim
    File    reference_seq
  }

  command {
	set -euo pipefail
	#fastqc ${read1_trim} -o .
  bowtie2-build ${reference_seq} mumps_ref
  bowtie2 -x mumps_ref -U ${read1_trim} -S ${sra_id}.sam --local
  ls
  ls>ls.txt
  ls /data
  ls /data> ls_data.txt
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
