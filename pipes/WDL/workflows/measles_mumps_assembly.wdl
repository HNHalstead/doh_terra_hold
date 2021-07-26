version 1.0

import "tasks/tasks_trim.wdl" as mm_trim
import "tasks/tasks_fastqc.wdl" as mm_fastqc
import "tasks/tasks_bowtie2.wdl" as mm_bowtie2
import "tasks/tasks_samtools.wdl" as mm_samtools

workflow mm_trim_and_assemble {

  input {
    String    sra_id
    File      read1
    File      read2
    File      reference_seq
  }

  call mm_trim.trim {
    input:
      sra_id=sra_id,
      read1=read1,
      read2=read2
  }

  call mm_fastqc.fastqc {
    input:
      sra_id=sra_id,
      read1_trim=trim.read1_trim,
      read2_trim=trim.read2_trim
  }

  call mm_bowtie2.bowtie2_se {
    input:
      sra_id=sra_id,
      read1_trim=trim.read1_trim,
      read2_trim=trim.read2_trim,
      reference_seq=reference_seq
  }

  call mm_samtools.sam_to_bam {
    input:
      sra_id=sra_id,
      samfile=bowtie2_se.samfile
  }

  output {
    File    read1_trim=trim.read1_trim
    File    read2_trim=trim.read2_trim
    File    fastqc_html=fastqc.fastqc_html
    #File    sam_file=bowtie2_se.samfile
    File    bamfile=sam_to_bam.bamfile
    File	sorted_bam=sam_to_bam.sorted_bam
    File	indexed_bam=sam_to_bam.indexed_bam
   #Float   gc_content=fastqc.gc_content

  }
}
