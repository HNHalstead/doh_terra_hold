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

  call mm_qc.stats_n_coverage {
    input:
      sra_id=sra_id,
      bamfile=sam_to_bam.bamfile
  }

  call mm_assembly_consensus.consensus {
    input:
      sra_id=sra_id,
      sorted_bam=sam_to_bam.sorted_bam,
      reference_seq=reference_seq
  }

  output {
    File    r1_fastqc_html_raw=fastqc_raw_r1.fastqc_html
    File    r2_fastqc_html_raw=fastqc_raw_r2.fastqc_html
    String  r1_total_sequences_raw=fastqc_raw_r1.total_sequences
    String  r2_total_sequences_raw=fastqc_raw_r2.total_sequences
    File    read1_trim=trim.read1_trim
    File    read2_trim=trim.read2_trim
    File    r1_fastqc_html_trimmed=fastqc_trim_r1.fastqc_html
    File    r2_fastqc_html_trimmed=fastqc_trim_r2.fastqc_html
    String  r1_total_sequences_trimmed=fastqc_trim_r1.total_sequences
    String  r2_total_sequences_trimmed=fastqc_trim_r2.total_sequences
    String  r1_percent_gc_trimmed=fastqc_trim_r1.percent_gc
    String  r2_percent_gc_trimmed=fastqc_trim_r2.percent_gc
    #File    sam_file=bowtie2_se.samfile
    File    bamfile=sam_to_bam.bamfile
    File    sorted_bam=sam_to_bam.sorted_bam
    File    indexed_bam=sam_to_bam.indexed_bam
    File    kraken2_report=kraken2.kraken_report
    Float   percent_human=kraken2.percent_human

    Int     consensus_number_N=consensus.number_N
    Int     consensus_number_ATCG=consensus.number_ATCG
    Int     consensus_number_Degenerate=consensus.number_Degenerate
    Int     consensus_number_total=consensus.number_Total
    File    consensus_seq=consensus.consensus_seq
    File    consensus_qual_File=consensus.consensus_qual

    File       assembly_stats=stats_n_coverage.stats
    File       assembly_cov_hist=stats_n_coverage.assembly_cov_hist
    File       assembly_cov_stats=stats_n_coverage.assembly_cov_stats
    File       assembly_flagstat=stats_n_coverage.flagstat
    Float      assembly_coverage=stats_n_coverage.coverage
    Float      assembly_depth=stats_n_coverage.depth
    Float      assembly_meanbaseq=stats_n_coverage.meanbaseq
    Float      assemlbymeanmapq=stats_n_coverage.meanmapq

    #File    sample_variants =consensus.sample_variants
    #String  variant_num=consensus.variant_num

  }
}
