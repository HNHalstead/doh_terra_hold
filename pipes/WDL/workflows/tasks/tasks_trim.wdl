version 1.0

task trim {

  input {
    String    sra_id
    File      read1
    File      read2
    Int?      cpus = 4
    String      memory = "8 GB"
    String  docker_image="staphb/trimmomatic:0.39"
  }

  command {
    date | tee DATE
    trimmomatic --version | head -n1 | tee VERSION
    trimmomatic_v=$(cat VERSION)

    trimmomatic SE -phred33 ${read1} ${sra_id}_trimmed_1.fastq.gz \
    ILLUMINACLIP:../Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    trimmomatic SE -phred33 ${read2} ${sra_id}_trimmed_2.fastq.gz \
    ILLUMINACLIP:../Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    cat DATE>trim_software.txt
    echo -e "docker image\t${docker_image}">>trim_software.txt
    echo -e "docker image platform:">>trim_software.txt
    uname -a>>trim_software.txt
    echo -e "main tool used:">>trim_software.txt
    echo -e "strimmomatic\t$trimmomatic_v\t\ta flexible read trimming tool for Illumina NGS data">>trim_software.txt
    echo -e "licenses available at:">>trim_software.txt
    echo -e "\thttps://academic.oup.com/bioinformatics/article/30/15/2114/2390096">>trim_software.txt
    printf '%100s\n' | tr ' ' ->>trim_software.txt
    dpkg -l>>trim_software.txt

  }

  output {
    File    read1_trim="${sra_id}_trimmed_1.fastq.gz"
    File    read2_trim="${sra_id}_trimmed_2.fastq.gz"
    File	  image_software="trim_software.txt"
    String     date          = read_string("DATE")
    String     version       = read_string("VERSION")
  }

  runtime {
    docker:       "${docker_image}"
    memory:       "${memory}"
    cpu:          cpus
    disks:        "local-disk 100 SSD"
    preemptible:  1
  }
}
