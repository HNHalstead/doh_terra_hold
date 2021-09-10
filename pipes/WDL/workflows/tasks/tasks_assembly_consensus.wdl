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

    dpkg -l>software.txt
  }

  output {
    File    samfile="${sra_id}.sam"
    File    image_software="software.txt"
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

task bowtie2_se_to_bam {
  input {
    String  sra_id
    File    read1_trim
    File    read2_trim
    File    reference_seq
    String  docker_image="hnh0303/bowtie2:0.1.0"
    String  alias_name="bowtie2_se_bam"
    Int     cpus=4
    String  memory="8 GB"
  }

  command {
    date | tee BOWTIE2DATE
    bowtie2 --version | head -n1 | tee BOWTIE2VERSION
    bowtie2_v=$(cat BOWTIE2VERSION)

    date | tee SAMTOOLSDATE
    samtools --version | head -n1 | tee SAMTOOLSVERSION
    samtools_v=$(cat SAMTOOLSVERSION)

    bowtie2-build ${reference_seq} mumps_ref
    bowtie2 -x mumps_ref -U ${read1_trim},${read2_trim} --local | samtools view -S -b >${sra_id}.bam
    samtools sort ${sra_id}.bam -o ${sra_id}.sorted.bam
    samtools index ${sra_id}.sorted.bam

    cat BOWTIE2DATE>bowtie2_se_to_bam_software.txt
    echo -e "docker image\t${docker_image}">>bowtie2_se_to_bam_software.txt
    echo -e "docker image platform:">>bowtie2_se_to_bam_software.txt
    uname -a>>bowtie2_se_to_bam_software.txt
    echo -e "main tool(s) used:">>bowtie2_se_to_bam_software.txt
    echo -e "\tbowtie2\t$bowtie2_v\t\tsequence alignment and sequence analysis">>bowtie2_se_to_bam_software.txt
    echo -e "\tsamtools\t$samtools_v\t\tset of utilities for interacting with and post-processing short DNA sequence read alignments">>bowtie2_se_to_bam_software.txt
    echo -e "licenses available at:">>bowtie2_se_to_bam_software.txt
    echo -e "\thttps://github.com/BenLangmead/bowtie2/blob/master/LICENSE">>bowtie2_se_to_bam_software.txt
    echo -e "\thttps://github.com/samtools/samtools/blob/develop/LICENSE">>bowtie2_se_to_bam_software.txt
    printf '%100s\n' | tr ' ' ->>bowtie2_se_to_bam_software.txt
    dpkg -l>>bowtie2_se_to_bam_software.txt

  }

  output {
    File    bamfile="${sra_id}.bam"
    File	sorted_bam="${sra_id}.sorted.bam"
    File	indexed_bam="${sra_id}.sorted.bam.bai"
    File    image_software="bowtie2_se_to_bam_software.txt"
    String     assembly_bowtie2_date=read_string("BOWTIE2DATE")
    String     assembly_bowtie2_version=read_string("BOWTIE2VERSION")
    String     assembly_samtools_date=read_string("SAMTOOLSDATE")
    String     assembly_samtools_version=read_string("SAMTOOLSVERSION")
  }

  runtime {
    docker:       "${docker_image}"
    memory:       "${memory}"
    cpu:          cpus
    disks:        "local-disk 100 SSD"
    preemptible:  1
  }
}

task sam_to_bam {

  input {
    String    sra_id
    File      samfile
    String    docker_image="staphb/samtools:1.12"
    String    memory="8 GB"
    Int       cpus=4
  }

  command {
    date | tee DATE
    samtools --version | head -n1 | tee VERSION
    samtools_v=$(cat VERSION)


    samtools view -S -b ${samfile}>${sra_id}.bam
    samtools sort ${sra_id}.bam -o ${sra_id}.sorted.bam
    samtools index ${sra_id}.sorted.bam

    cat DATE>sam_to_bam_software.txt
    echo -e "docker image\t${docker_image}">>sam_to_bam_software.txt
    echo -e "docker image platform:">>sam_to_bam_software.txt
    uname -a>>sam_to_bam_software.txt
    echo -e "main tool used:">>sam_to_bam_software.txt
    echo -e "samtools\t$samtools_v\t\tset of utilities for interacting with and post-processing short DNA sequence read alignments">>sam_to_bam_software.txt
    echo -e "licenses available at:">>sam_to_bam_software.txt
    echo -e "\thttps://github.com/samtools/samtools/blob/develop/LICENSE">>sam_to_bam_software.txt
    printf '%100s\n' | tr ' ' ->>sam_to_bam_software.txt
    dpkg -l>>sam_to_bam_software.txt
  }

  output {
    File    bamfile="${sra_id}.bam"
    File	sorted_bam="${sra_id}.sorted.bam"
    File	indexed_bam="${sra_id}.sorted.bam.bai"
    String     date          = read_string("DATE")
    String     version       = read_string("VERSION")
    File    image_software="sam_to_bam_software.txt"
  }

  runtime {
    docker:       "${docker_image}"
    memory:       "${memory}"
    cpu:          cpus
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
    String  docker_image="staphb/ivar:1.3"
    String  alias_name="bowtie2_se_bam"
    Int     cpus=4
    String  memory="8 GB"
  }


  command {
    date | tee DATE
    samtools --version | head -n1 | tee VERSION
    samtools_v=$(cat VERSION)
    ivar version | head -n1 | tee IVAR_VERSION
    ivar_v=$(cat IVAR_VERSION)

    samtools mpileup -d ${max_depth} -Q ${min_bq} --reference ${reference_seq} ${sorted_bam} | ivar consensus -p ${sra_id}_consensus -q ${min_qual} -t ${min_freq} -m ${min_depth} -n N

   num_N=$( grep -v ">" ${sra_id}_consensus.fa | grep -o 'N' | wc -l )
   if [ -z "$num_N" ] ; then num_N="0" ; fi
   echo $num_N | tee NUM_N

   num_ACTG=$( grep -v ">" ${sra_id}_consensus.fa | grep -o -E "C|A|T|G" | wc -l )
   if [ -z "$num_ACTG" ] ; then num_ACTG="0" ; fi
   echo $num_ACTG | tee NUM_ACTG

   num_degenerate=$( grep -v ">" ${sra_id}_consensus.fa | grep -o -E "B|D|E|F|H|I|J|K|L|M|O|P|Q|R|S|U|V|W|X|Y|Z" | wc -l )
   if [ -z "$num_degenerate" ] ; then num_degenerate="0" ; fi
   echo $num_degenerate | tee NUM_DEGENERATE

   num_total=$( grep -v ">" ${sra_id}_consensus.fa | grep -o -E '[A-Z]' | wc -l )
   if [ -z "$num_total" ] ; then num_total="0" ; fi
   echo $num_total | tee NUM_TOTAL

   cat DATE>consensus_software.txt
   echo -e "docker image\t${docker_image}">>consensus_software.txt
   echo -e "docker image platform:">>consensus_software.txt
   uname -a>>consensus_software.txt
   echo -e "main tool(s) used:">>consensus_software.txt
   echo -e "ivar\t$ivar_v\t\tpackage that contains functions broadly useful for viral amplicon-based sequencing">>consensus_software.txt
   echo -e "samtools\t$samtools_v\t\tset of utilities for interacting with and post-processing short DNA sequence read alignments">>consensus_software.txt
   echo -e "licenses available at:">>consensus_software.txt
   echo -e "\thttps://github.com/andersen-lab/ivar/blob/master/LICENSE">>consensus_software.txt
   echo -e "\thttps://github.com/samtools/samtools/blob/develop/LICENSE">>consensus_software.txt
   printf '%100s\n' | tr ' ' ->>consensus_software.txt
   dpkg -l>>consensus_software.txt


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
    File      consensus_seq = "${sra_id}_consensus.fa"
    File      consensus_qual = "${sra_id}_consensus.qual.txt"
    String     date = read_string("DATE")
    String     samtools_version = read_string("VERSION")
    String     ivar_version = read_string("IVAR_VERSION")
    File       image_software="consensus_software.txt"
    #String     variant_num       = read_string("VARIANT_NUM")
  }

  runtime {
    docker:       "${docker_image}"
    memory:       "${memory}"
    cpu:          cpus
    disks:        "local-disk 100 SSD"
    preemptible:  1
  }
  }
