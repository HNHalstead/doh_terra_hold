version 1.0

task fastqc {

  input {
    String    sra_id
    File      read1_trim
    File      read2_trim
  }

  command {
	set -euo pipefail
	fastqc ${read1_trim} -o . /fastqc ${read2_trim} -o .
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

task kraken2 {
  input {
  	File        read1
	  File? 		  read2
	  String      sra_id
	  String?     kraken2_db = "/kraken2-db"
    Int?        cpus=4
    String      virus_name="Mumps"
  }

  command{
    date | tee DATE
    kraken2 --version | head -n1 | tee VERSION
    num_reads=$(ls *fastq.gz 2> /dev/nul | wc -l)
    if ! [ -z ${read2} ]; then
      mode="--paired"
    fi
    echo $mode
    kraken2 $mode \
      --classified-out cseqs#.fq \
      --threads ${cpus} \
      --db ${kraken2_db} \
      ${read1} ${read2} \
      --report ${sra_id}_kraken2_report.txt

    percentage_human=$(grep "Homo sapiens" ${sra_id}_kraken2_report.txt | cut -f 1)
     # | tee PERCENT_HUMAN
    percentage_virus=$(grep ${virus_name} ${sra_id}_kraken2_report.txt | cut -f1 )
    echo $virus "grepped virus" $percentage_virus
     # | tee PERCENT_COV
    if [ -z "$percentage_human" ] ; then percentage_human="0.00" ; fi
    if [ -z "$percentage_virus" ] ; then percentage_virus="0.00" ; fi
  }

  output {
    String     date          = read_string("DATE")
    String     version       = read_string("VERSION")
    File 	     kraken_report = "${sra_id}_kraken2_report.txt"
    Float 	   percent_human = read_string("PERCENT_HUMAN")
    Float 	   percent_virus   = read_string("PERCENT_VIRUS")
  }

  runtime {
    docker:       "staphb/kraken2:2.0.9-beta"
    memory:       "16 GB"
    cpu:          4
    disks:        "local-disk 100 SSD"
    preemptible:  0
  }
}
