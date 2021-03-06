version 1.0

task bowtie2_se {
  input {
    String  sra_id
    File    read1_trim
    File    read2_trim
    File    reference_seq
  }

  command {
    date | tee DATE
    bowtie2 --version | head -n1 | tee VERSION

    bowtie2-build ${reference_seq} mumps_ref
    bowtie2 -x mumps_ref -U ${read1_trim},${read2_trim} -S ${sra_id}.sam --local
  }

  output {
    File    samfile="${sra_id}.sam"
    String     date          = read_string("DATE")
    String     version       = read_string("VERSION")
  }

  runtime {
    docker:       "quay.io/biocontainers/bowtie2:2.4.4--py38h72fc82f_0"
    memory:       "8 GB"
    cpu:          4
    disks:        "local-disk 100 SSD"
    preemptible:  1
  }
}

task sam_to_bam {

  input {
    String    sra_id
    File      samfile
  }

  command {
    date | tee DATE
    samtools --version | head -n1 | tee VERSION


    samtools view -S -b ${samfile}>${sra_id}.bam
    samtools sort ${sra_id}.bam -o ${sra_id}.sorted.bam
    samtools index ${sra_id}.sorted.bam
  }

  output {
    File    bamfile="${sra_id}.bam"
    File	sorted_bam="${sra_id}.sorted.bam"
    File	indexed_bam="${sra_id}.sorted.bam.bai"
    String     date          = read_string("DATE")
    String     version       = read_string("VERSION")
  }

  runtime {
    docker:       "staphb/samtools:1.12"
    memory:       "8 GB"
    cpu:          4
    disks:        "local-disk 100 SSD"
    preemptible:  1
  }
}

task consensus {

  input {
    String    sra_id
    File      sorted_bam
    File      reference_seq
    Int         max_depth = "60000"
    Boolean     disable_baq = true
    Int         min_bq = "0"
    Int         min_qual = "20"
    Float       min_freq = "0.6"
    Int         min_depth = "10"
  }


  command {
    date | tee DATE
    samtools --version | head -n1 | tee VERSION
    ivar version | head -n1 | tee IVAR_VERSION

    samtools mpileup -d ${max_depth} -Q ${min_bq} --reference ${reference_seq} ${sorted_bam} | ivar consensus -p ${sra_id}.consensus -q ${min_qual} -t ${min_freq} -m ${min_depth} -n N

    ls
    ls>ls.txt

   num_N=$( grep -v ">" ${sra_id}.consensus.fa | grep -o 'N' | wc -l )
   if [ -z "$num_N" ] ; then num_N="0" ; fi
   echo $num_N | tee NUM_N

   num_ACTG=$( grep -v ">" ${sra_id}.consensus.fa | grep -o -E "C|A|T|G" | wc -l )
   if [ -z "$num_ACTG" ] ; then num_ACTG="0" ; fi
   echo $num_ACTG | tee NUM_ACTG

   num_degenerate=$( grep -v ">" ${sra_id}.consensus.fa | grep -o -E "B|D|E|F|H|I|J|K|L|M|O|P|Q|R|S|U|V|W|X|Y|Z" | wc -l )
   if [ -z "$num_degenerate" ] ; then num_degenerate="0" ; fi
   echo $num_degenerate | tee NUM_DEGENERATE

   num_total=$( grep -v ">" ${sra_id}.consensus.fa | grep -o -E '[A-Z]' | wc -l )
   if [ -z "$num_total" ] ; then num_total="0" ; fi
   echo $num_total | tee NUM_TOTAL

    #variants_num=$(grep "TRUE" ${sra_id}.variants.tsv | wc -l)
    #if [ -z "$variants_num" ] ; then variants_num="0" ; fi
    #echo $variants_num | tee VARIANT_NUM

  }

  output {
    Int       number_N = read_string("NUM_N")
    Int       number_ATCG = read_string("NUM_ACTG")
    Int       number_Degenerate = read_string("NUM_DEGENERATE")
    Int       number_Total = read_string("NUM_TOTAL")
    #File      sample_variants = "${sra_id}.variants.tsv"
    File      consensus_seq = "${sra_id}.consensus.fa"
    File      consensus_qual = "${sra_id}.consensus.qual.txt"
    String     date = read_string("DATE")
    String     samtools_version = read_string("VERSION")
    String     ivar_version = read_string("IVAR_VERSION")
    #String     variant_num       = read_string("VARIANT_NUM")
  }

  runtime {
    docker:       "staphb/ivar:1.3"
    memory:       "8 GB"
    cpu:          4
    disks:        "local-disk 100 SSD"
    preemptible:  1
  }
  }
