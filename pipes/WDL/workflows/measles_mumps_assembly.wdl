version 1.0

import "tasks/tasks_trim.wdl" as mm_trim
import "tasks/tasks_qc.wdl" as mm_qc
import "tasks/tasks_assembly_consensus.wdl" as mm_assembly_consensus

workflow mm_trim_and_assemble {

  input {
    String    sra_id
    File      read1
    File      read2
    File      reference_seq
    String      virus_name="Mumps"
  }

  call mm_qc.fastqc_se as fastqc_raw_r1 {
    input:
      sra_id=sra_id,
      read=read1,
  }

  call mm_qc.fastqc_se as fastqc_raw_r2 {
    input:
      sra_id=sra_id,
      read=read2
  }

  call mm_trim.trim {
    input:
      sra_id=sra_id,
      read1=read1,
      read2=read2
  }

  call mm_qc.fastqc_se as fastqc_trim_r1 {
    input:
      sra_id=sra_id,
      read=trim.read1_trim
  }

  call mm_qc.fastqc_se as fastqc_trim_r2 {
    input:
      sra_id=sra_id,
      read=trim.read2_trim
  }

  call mm_qc.kraken2 {
    input:
      sra_id=sra_id,
      read1=read1,
      read2=read2,
      virus_name=virus_name
  }

  call mm_assembly_consensus.bowtie2_se {
    input:
      sra_id=sra_id,
      read1_trim=trim.read1_trim,
      read2_trim=trim.read2_trim,
      reference_seq=reference_seq
  }

  call mm_assembly_consensus.sam_to_bam {
    input:
      sra_id=sra_id,
      samfile=bowtie2_se.samfile
  }

  output {
    File    r1_fastqc_html_raw=fastqc_raw_r1.fastqc_html
    File    r2_fastqc_html_raw=fastqc_raw_r2.fastqc_html
    String  r1_total_sequences_raw=fastqc_raw_r1.total_sequences
    String  r2_total_sequences_raw=fastqc_raw_r2.total_sequences

    File    read1_trim=trim.read1_trim
    File    read2_trim=trim.read2_trim
    File    fastqc_html_r1=fastqc_trim_r1.fastqc_html
    File    fastqc_html_r2=fastqc_trim_r2.fastqc_html
    String  total_sequences_trimmed_r1=fastqc_trim_r1.total_sequences
    String  total_sequences_trimmed_r2=fastqc_trim_r2.total_sequences
    #File    sam_file=bowtie2_se.samfile
    File    bamfile=sam_to_bam.bamfile
    File    sorted_bam=sam_to_bam.sorted_bam
    File    indexed_bam=sam_to_bam.indexed_bam
    File    kraken2_report=kraken2.kraken_report
    Float   percent_human=kraken2.percent_human
    Float   percent_virus=kraken2.percent_virus
   #Float   gc_content=fastqc.gc_content

  }
}
