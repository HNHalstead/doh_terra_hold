version 1.0

task trim {

  input {
    String    sra_id
    File      read1
    File      read2
  }

  command {
    trimmomatic SE -phred33 ${read1} ${sra_id}_trimmed_1.fastq.gz \
    ILLUMINACLIP:../Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    trimmomatic SE -phred33 ${read2} ${sra_id}_trimmed_2.fastq.gz \
    ILLUMINACLIP:../Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

  }

  output {
    File    read1_trim="${sra_id}_trimmed_1.fastq.gz"
    File    read2_trim="${sra_id}_trimmed_2.fastq.gz"
  }

  runtime {
    docker:       "staphb/trimmomatic:0.39"
    memory:       "8 GB"
    cpu:          4
    disks:        "local-disk 100 SSD"
    preemptible:  1
  }
}
