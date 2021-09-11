version 1.0

task fastqc_pe {
  input {
    String    sra_id
    File      read1
    File      read2
    Int?      cpus = 2
    Int?      mem_gb = 8
  }

  command {
    set -euo pipefail
    fastqc ${read1} -o $PWD --threads ${cpus}/fastqc ${read2} -o $PWD --threads ${cpus}
    date | tee DATE
    fastqc --version | hgrep FastQC | tee VERSION
    ls
    ls>ls.txt
  }

  output {
    File  fastqc_html_r1=glob("*fastqc.html")[0]
    File	fastqc_zip_r1=glob("*fastqc.zip")[0]
    File  fastqc_html_r2=glob("*fastqc.html")[1]
    File	fastqc_zip_r2=glob("*fastqc.zip")[1]
    String	date=read_string("DATE")
    String	version=read_string("VERSION")
  }

  runtime {
    docker:       "staphb/fastqc:0.11.8"
    memory:       mem_gb + "GB"
    cpu:          4
    disks:        "local-disk 100 SSD"
    preemptible:  1
    continueOnReturnCode: "True"
  }
}

task fastqc_se {
  input {
    String    sra_id
    File      read
    Int?      cpus = 4
    String      memory = "16 GB"
    String stripped = basename(read, ".fastq.gz")
    String  docker_image="staphb/fastqc:0.11.8"
  }

  command {
    set -euo pipefail
    fastqc ${read} -o $PWD --threads ${cpus}
    date | tee DATE
    fastqc --version | tr -d 'FastQC' | tee VERSION
    fastqc_v=$(cat VERSION)

    unzip -p ${stripped}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 | tee READ_SEQS

    unzip -p ${stripped}_fastqc.zip */fastqc_data.txt | grep "%GC" | cut -f 2 | tee PERCENT_GC

    cat DATE>fastqc_se_software.txt
    echo -e "docker image:\t${docker_image}">>fastqc_se_software.txt
    echo -e "docker image platform:">>fastqc_se_software.txt
    uname -a>>fastqc_se_software.txt
    echo -e "main tool used:">>fastqc_se_software.txt
    echo -e "\tFastQC\t$fastqc_v\t\ta quality control tool for high throughput sequence data">>fastqc_se_software.txt
    echo -e "licenses available at:">>fastqc_se_software.txt
    echo -e "\thttps://github.com/s-andrews/FastQC/blob/master/LICENSE">>fastqc_se_software.txt
    printf '%100s\n' | tr ' ' ->>fastqc_se_software.txt
    dpkg -l>>fastqc_se_software.txt


  }

  output {
    File  fastqc_html=glob("*fastqc.html")[0]
    File	fastqc_zip=glob("*fastqc.zip")[0]
    File	image_software="fastqc_se_software.txt"
    String	date=read_string("DATE")
    String	version=read_string("VERSION")
    String	total_sequences=read_string("READ_SEQS")
    String  percent_gc=read_string("PERCENT_GC")
  }

  runtime {
    docker:       "${docker_image}"
    memory:       "${memory}"
    cpu:          cpus
    disks:        "local-disk 100 SSD"
    preemptible:  0
    continueOnReturnCode: "True"
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
    Int?        cpus = 4
    String      memory = "16 GB"
    String      docker_image="staphb/kraken2:2.0.9-beta"
  }

  command{
    date | tee DATE
    kraken2 --version | head -n1 | tee VERSION
    kraken2_v=$(cat VERSION)
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
    echo $percentage_human | tee PERCENT_HUMAN
    echo $percentage_virus | tee PERCENT_VIRUS

    cat DATE>kraken2_software.txt
    echo -e "docker image:\t${docker_image}">>kraken2_software.txt
    echo -e "docker image platform:">>kraken2_software.txt
    uname -a>>kraken2_software.txt
    echo -e "main tool used:">>kraken2_software.txt
    echo -e "\tkraken2\t$kraken2_v\t\ta taxonomic classification tool using exact k-mer matches to achieve high accuracy and fast classification speeds">>kraken2_software.txt
    echo -e "licenses available at:">>kraken2_software.txt
    echo -e "\thttps://github.com/DerrickWood/kraken2/blob/master/LICENSE">>kraken2_software.txt
    printf '%100s\n' | tr ' ' ->>kraken2_software.txt
    dpkg -l>>kraken2_software.txt
  }

  output {
    String     date          = read_string("DATE")
    String     version       = read_string("VERSION")
    File 	     kraken_report = "${sra_id}_kraken2_report.txt"
    Float 	   percent_human = read_string("PERCENT_HUMAN")
    Float 	   percent_virus   = read_string("PERCENT_VIRUS")
    File	     image_software="kraken2_software.txt"
  }

  runtime {
    docker:       "${docker_image}"
    memory:       "${memory}"
    cpu:          cpus
    disks:        "local-disk 100 SSD"
    preemptible:  0
  }
}

task stats_n_coverage {

  input {
    File        bamfile
    String      sra_id
    String      docker_image="staphb/samtools:1.10"
    Int?        cpus = 4
    String?     memory = "16 GB"
  }

  command{
    date | tee DATE
    samtools --version | head -n1 | tee VERSION
    samtools_v=$(cat VERSION)

    samtools stats ${bamfile} > ${sra_id}.stats.txt

    samtools coverage ${bamfile} -m -o ${sra_id}.cov.hist
    samtools coverage ${bamfile} -o ${sra_id}.cov.txt
    samtools flagstat ${bamfile} > ${sra_id}.flagstat.txt

    coverage=$(cut -f 6 ${sra_id}.cov.txt | tail -n 1)
    depth=$(cut -f 7 ${sra_id}.cov.txt | tail -n 1)
    meanbaseq=$(cut -f 8 ${sra_id}.cov.txt | tail -n 1)
    meanmapq=$(cut -f 9 ${sra_id}.cov.txt | tail -n 1)

    if [ -z "$coverage" ] ; then coverage="0" ; fi
    if [ -z "$depth" ] ; then depth="0" ; fi
    if [ -z "$meanbaseq" ] ; then meanbaseq="0" ; fi
    if [ -z "$meanmapq" ] ; then meanmapq="0" ; fi

    echo $coverage | tee COVERAGE
    echo $depth | tee DEPTH
    echo $meanbaseq | tee MEANBASEQ
    echo $meanmapq | tee MEANMAPQ

    cat DATE>stats_n_coverage_software.txt
    echo -e "docker image:\t${docker_image}">>stats_n_coverage_software.txt
    echo -e "docker image platform:">>stats_n_coverage_software.txt
    uname -a>>stats_n_coverage_software.txt
    echo -e "main tool used:">>stats_n_coverage_software.txt
    echo -e "samtools\t$samtools_v\t\tset of utilities for interacting with and post-processing short DNA sequence read alignments">>stats_n_coverage_software.txt
    echo -e "licenses available at:">>stats_n_coverage_software.txt
    echo -e "\thttps://github.com/samtools/samtools/blob/develop/LICENSE">>stats_n_coverage_software.txt
    printf '%100s\n' | tr ' ' ->>stats_n_coverage_software.txt
    dpkg -l>>stats_n_coverage_software.txt
  }

  output {
    String     date = read_string("DATE")
    String     samtools_version = read_string("VERSION")
    File       stats = "${sra_id}.stats.txt"
    File       assembly_cov_hist = "${sra_id}.cov.hist"
    File       assembly_cov_stats = "${sra_id}.cov.txt"
    File       flagstat = "${sra_id}.flagstat.txt"
    Float      coverage = read_string("COVERAGE")
    Float      depth = read_string("DEPTH")
    Float      meanbaseq = read_string("MEANBASEQ")
    Float      meanmapq = read_string("MEANMAPQ")
    File	     image_software="stats_n_coverage_software.txt"
  }

  runtime {
    docker:       "~{docker_image}"
    memory:       "~{memory}"
    cpu:          cpus
    disks:        "local-disk 100 SSD"
    preemptible:  0
  }
  meta {
        author: "Edited from StaPH-B Titan pipeline for SARS-CoV-2 analysis"
    }
}
