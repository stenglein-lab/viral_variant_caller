#!/usr/bin/env nextflow

/*
    A pipeline to determine SARS-CoV-2 consensus sequences and prepare
    information for submissiont to GISAID and other public databases

    March 5, 2022

    Mark Stenglein
*/

/*
  Check input parameters
*/
def check_params_and_input () {

  // TODO:
  // check_metadata()
  // check indexes exist

  // This list includes a list of files or paths that are required
  // to exist.  Check that they exist and fail if not.
  checkPathParamList = [
    params.input_dir,
    // params.fastq_dir,
    params.refseq_dir,
    params.script_dir,
    params.refseq_fasta,
    params.refseq_genbank,
    params.primer_fasta_5p,
    params.primer_fasta_3p
  ]
  log.info("Checking for required input paths and files...")
  for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

}

/* 
  Check parameters and input
*/
check_params_and_input()

// possible multiple directories
// if fastq are in multiple directories the directories
// should be provided as a comma-separated list (no spaces)
def fastq_dirs = params.fastq_dir.tokenize(',')

// construct list of directories in which to find fastq
fastq_dir_list = []
for (dir in fastq_dirs){
   def file_pattern = "${dir}/${params.fastq_pattern}"
   fastq_dir_list.add(file_pattern)
}


/*
 These fastq files represent the main input to this workflow
*/
Channel
  // .fromFilePairs("${params.fastq_dir}/${params.fastq_pattern}", size: -1, checkIfExists: true, maxDepth: 1)
  .fromFilePairs(fastq_dir_list, size: -1, checkIfExists: true, maxDepth: 1)

  // this map gets rid of any _S\d+ at the end of sample IDs but leaves fastq
  // names alone.  E.g. strip _S1 from the end of a sample ID..  
  // This is typically sample #s from Illumina basecalling.
  // could cause an issue if sequenced the same sample with 
  // multiple barcodes so was repeated on a sample sheet. 
  .map{ [ it[0].replaceFirst( /_S\d+$/ ,"") , it[1] ] }                           
  .into {samples_ch_qc; samples_ch_trim; samples_ch_count; sample_ids_ch}


/*
  Setup some initial indexes and dictionaries needed by downstream processes.
  Only do this once at beginning.

  All this setup used to be in a single process but I separated it out because
  each process needs its own singularity image
  
  Chain the processes together so at the end of the chain all indexes will be
  complete.

*/
process setup_lofreq_index {

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py39h43839c7_5"
  } else {                                                                      
      container "quay.io/biocontainers/lofreq:2.1.5--py39h43839c7_5"
  }      

  output:
  val("lofreq_index_complete") into post_lofreq_index_ch

  script:
  """
  # ------------------
  # lofreq fasta index
  # ------------------
  lofreq faidx ${params.refseq_fasta}
  """
}


process setup_snpeff_index {

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/snpeff:5.0--hdfd78af_1"
  } else {                                                                      
      container "quay.io/biocontainers/snpeff:5.0--hdfd78af_1"
  }      

  input:
  val("lofreq_index_complete") from post_lofreq_index_ch

  output:
  val("snpeff_index_complete") into post_snpeff_index_ch

  script:
  """
  # ----------------
  # setup snpEFF db
  # ----------------

  # make a minimal snpEff.config 
  rm -f ${params.snpeff_cfg}
  echo "${params.refseq_name}.genome: ${params.refseq_name}" > ${params.snpeff_cfg}

  # make a directory for the snp eff db
  mkdir -p ${params.snpeff_data}/${params.refseq_name}

  # cp fasta and genbank format data for virus refseq to the directory location snpeff is expecting
  cp ${params.refseq_fasta} ${params.snpeff_data}/${params.refseq_name}/sequences.fa
  cp ${params.refseq_genbank} ${params.snpeff_data}/${params.refseq_name}/genes.gbk

  # could make it from genbank format file
  snpEff build -c ${params.snpeff_cfg} -nodownload -v -genbank -dataDir ${params.snpeff_data} ${params.refseq_name} > ${params.snpeff_data}/${params.refseq_name}.build
  """
}

process setup_bwa_index {

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/bwa:0.7.8--h7132678_7"
  } else {                                                                      
      container "quay.io/biocontainers/bwa:0.7.8--h7132678_7"
  }      

  input:
  val("snpeff_index_complete") from post_snpeff_index_ch

  output:
  val("bwa_index_complete") into post_bwa_index_ch

  script:
  """
  # -----------------------
  # bwa index viral refseq
  # -----------------------
  bwa index ${params.refseq_fasta}
  """
}

process setup_gatk_index {

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/gatk4:4.2.6.0--hdfd78af_0"
  } else {                                                                      
      container "quay.io/biocontainers/gatk4:4.2.6.0--hdfd78af_0"
  }      

  input:
  val("bwa_index_complete") from post_bwa_index_ch

  output:
  val("all_indexes_complete") into post_index_setup_ch

  script:
  """
  # -----------------
  # GATK index setup
  # -----------------
  # setup gatk indexes for BSQR

  # first, we need to make, and then index a dummy ignore_regions.bed file because gatk requires this fileb
  # this file will consist of the refseq name plus the coordinates 1 1, which is just the first base of the 
  # reference sequence.  GATK will ignore this base for basecall quality score recalibration 
  printf "%s\t1\t1\n" ${params.refseq_name} > ${params.ignore_regions}

  gatk IndexFeatureFile --input ${params.ignore_regions} 

  rm -f "${params.refseq_dir}/${params.refseq_name}.dict"
  gatk CreateSequenceDictionary -R ${params.refseq_fasta}

  """
}


/*
 Run fastqc on input fastq 
*/
process initial_qc {
  label 'lowmem_non_threaded'                                                                

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0"  
  } else {                                                                      
      container "quay.io/biocontainers/fastqc:0.11.9--0"                        
  }     

  input:
  tuple val(sample_id), path(initial_fastq) from samples_ch_qc

  output:
  val(sample_id) into post_initial_qc_ch

  script:
  """
  mkdir -p  ${params.initial_fastqc_dir} 
  fastqc -o ${params.initial_fastqc_dir} $initial_fastq 
  """
}

/*
 Count initial fastq
*/
process initial_fastq_count {
  label 'lowmem_non_threaded'                                                                
  publishDir "${params.counts_out_dir}", mode:'link'

  input:
  tuple val(sample_id), path(initial_fastq) from samples_ch_count

  output:
  path("${sample_id}_initial_count.txt") into post_count_initial_ch

  shell:
  // only count the first file because for paired-read data both files
  // will have the same # of reads.

  // for an explanation of the xargs command used for arithmetic in a pipe, see: 
  // https://askubuntu.com/questions/1203063/how-can-i-pipe-the-result-of-the-wc-command-into-an-arithmetic-expansion
  '''
  zcat -f !{initial_fastq[0]} | wc -l | xargs bash -c 'echo $(($0 / 4))' | awk '{print "!{sample_id}" "\tinitial\t" $1}' > "!{sample_id}_initial_count.txt"
  '''
}

/*
 Use multiqc to merge initial fastqc reports
*/
process initial_multiqc {
  publishDir "${params.outdir}", mode:'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/multiqc:1.11--pyhdfd78af_0"
  } else {                                                                      
      container "quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"              
  }    

  input:
  val(all_sample_ids) from post_initial_qc_ch.collect()

  output: 
  path("initial_qc_report.html")

  script:
  """
  multiqc -n "initial_qc_report.html" -m fastqc ${params.initial_fastqc_dir}
  """
}

/*
 Use cutadapt to trim off Illumina adapters and low quality bases
*/

process trim_illumina_adapters_and_low_quality {
  publishDir "${params.fastq_out_dir}", mode:'link'
  label 'lowmem_threaded'                                                                

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/cutadapt:3.5--py39h38f01e4_0"
  } else {                                                                      
    container "quay.io/biocontainers/cutadapt:3.5--py39h38f01e4_0"              
  }     

  input:
  tuple val(sample_id), path(initial_fastq) from samples_ch_trim

  output:
  tuple val(sample_id), path("*_f1.fastq") into samples_ch_second_trim

  script:

  // this handles paired-end data, in which case must specify a paired output file
  def paired_output   = initial_fastq[1] ? "-p ${sample_id}_R2_f1.fastq" : ""

  // fixed # of bases to trim
  def paired_trimming = initial_fastq[1] ? "-U $params.always_trim_5p_bases -U -${params.always_trim_3p_bases}" : ""

  // trim TruSeq-style cutadapt
  // see: https://cutadapt.readthedocs.io/en/stable/guide.html#illumina-truseq
  // def truseq_cutadapt = "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
  // def truseq_cutadapt = '-a "AGATCGGAAGAGC;min_overlap=8" -A "AGATCGGAAGAGC;min_overlap=8"'
  // def truseq_cutadapt = '-a "AGATCGGAAGAGC" -A "AGATCGGAAGAGC"'
  def truseq_cutadapt = '-a "AGATCGGAAGAGC;min_overlap=1" -A "AGATCGGAAGAGC;min_overlap=1"'

  """
  cutadapt \
   $truseq_cutadapt \
   -q 30,30 \
   --minimum-length ${params.post_trim_min_length} \
   -u ${params.always_trim_5p_bases} \
   -u -${params.always_trim_3p_bases} \
   $paired_trimming \
   -o ${sample_id}_R1_f1.fastq \
   --cores=${task.cpus} \
   $paired_output \
   $initial_fastq 
  """

}

/*
 Use cutadapt to trim off PCR primers
*/
process trim_PCR_primers {
  publishDir "${params.fastq_out_dir}", mode:'link'
  label 'lowmem_threaded'                                                                
 
  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/cutadapt:3.5--py39h38f01e4_0"
  } else {                                                                      
    container "quay.io/biocontainers/cutadapt:3.5--py39h38f01e4_0"              
  }     

  input:
  tuple val(sample_id), path(initial_fastq) from samples_ch_second_trim

  output:
  tuple val(sample_id), path("*_f.fastq") optional true into post_trim_qc_ch
  tuple val(sample_id), path("*_f.fastq") optional true into post_trim_ch
  tuple val(sample_id), path("*_f.fastq") optional true into post_trim_count_ch
  // path ("*cutadapt.json") 

  script:

  // this handles paired-end data, in which case must specify a paired output file
  def paired_output   = initial_fastq[1] ? "-p ${sample_id}_R2_f.fastq" : ""

  """
  cutadapt \
   -O 10 \
   -q 30,30 \
   --minimum-length ${params.post_trim_min_length} \
   -o ${sample_id}_R1_f.fastq \
   --json=${sample_id}.cutadapt.json \
   --cores=${task.cpus} \
   $paired_output \
   $initial_fastq 
  """
}

/*
 Count post-trimming fastq
*/
process trimmed_fastq_count {
  label 'lowmem_non_threaded'                                                                
  publishDir "${params.counts_out_dir}", mode:'link'

  input:
  tuple val(sample_id), path(trimmed_fastq) from post_trim_count_ch

  output:
  path("${sample_id}_trimmed_count.txt") into post_count_trim_ch

  shell:
  // only count the first file because for paired-read data both files
  // will have the same # of reads.
  '''
  zcat -f !{trimmed_fastq[0]} | wc -l | xargs bash -c 'echo $(($0 / 4))' |  awk '{print "!{sample_id}" "\tpost_trimming\t" $1}' > "!{sample_id}_trimmed_count.txt"
  '''
}

/*
 Use fastqc to do QC on post-trimmed fastq
*/
process post_trim_qc {
  label 'lowmem_non_threaded'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0"  
  } else {                                                                      
      container "quay.io/biocontainers/fastqc:0.11.9--0"                        
  }     

  input:
  tuple val(sample_id), path(input_fastq) from post_trim_qc_ch

  output:
  val(sample_id) into post_trim_multiqc_ch

  script:

  """
  mkdir -p  ${params.post_trim_fastqc_dir} 
  fastqc -o ${params.post_trim_fastqc_dir} $input_fastq
  """
}

/*
 Use multiqc to merge post-trimming fastq reports
*/
process post_trim_multiqc {
  publishDir "${params.outdir}", mode:'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/multiqc:1.11--pyhdfd78af_0"
  } else {                                                                      
      container "quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"              
  }    

  input:
  val(all_sample_ids) from post_trim_multiqc_ch.collect()

  output:
  path("post_trim_qc_report.html")

  """
  multiqc -n "post_trim_qc_report.html" -m fastqc -m cutadapt ${params.post_trim_fastqc_dir}
  """
}

/*
  Use bowtie2 to remove host-derived reads
*/
// TODO: make host filtering optional
// TODO: switch to bwa for host filtering too?  The reason to use bowtie2 in
//       instead of bwa is that it has convenient built-in options for 
//       outputting unmapped reads (the --un or --un-conc) options, whereas
//       for bwa you have to go through additional steps with samtools / bedtools 
process host_filtering {
  // publishDir "${params.outdir}", pattern: "*_R1_fh.fastq"
  label 'lowmem_threaded'                                                                

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/bowtie2:2.4.5--py39ha4319a6_1"
  } else {                                                                      
      container "quay.io/biocontainers/bowtie2:2.4.5--py39ha4319a6_1"           
  }        

  input:
  tuple val(sample_id), path(input_fastq) from post_trim_ch
  val("indexes_complete") from post_index_setup_ch

  output:
  tuple val(sample_id), path("*_fh.fastq") optional true into post_host_ch_variants
  tuple val(sample_id), path("*_fh.fastq") optional true into post_host_ch_dvg
  tuple val(sample_id), path("*_fh.fastq") optional true into post_host_ch_count


  // TODO: multiqc analysis of bowtie output (host filtering) (?)

  script:

  // handle single-end or paired-end inputs
  def r1 = input_fastq[0] 
  def r2 = input_fastq[1] ? input_fastq[1] : ""
  def bowtie_file_input  = input_fastq[1] ? "-1 $r1 -2 $r2" : "-U $r1"
  def bowtie_file_output = input_fastq[1] ? "--un-conc ${sample_id}_R%_fh.fastq" : "--un ${sample_id}_R1_fh.fastq"

  """
  bowtie2 \
  -x "${params.host_bt_index}" \
  --local \
  -q \
  $bowtie_file_input \
  --sensitive \
  --score-min "C,${params.host_bt_min_score},0" \
  -p ${task.cpus} \
  $bowtie_file_output > /dev/null 
  """
}


/*
 Count post-host-filtering fastq
*/
process host_filtered_fastq_count {
  label 'lowmem_non_threaded'                                                                
  publishDir "${params.counts_out_dir}", mode:'link'

  input:
  tuple val(sample_id), path(filtered_fastq) from post_host_ch_count

  output:
  path("${sample_id}_host_filtered_count.txt") into post_count_host_ch

  shell:
  // only count the first file because for paired-read data both files
  // will have the same # of reads.
  '''
  zcat -f !{filtered_fastq[0]} | wc -l | xargs bash -c 'echo $(($0 / 4))' |  awk '{print "!{sample_id}" "\tpost_host_filtered\t" $1}' > "!{sample_id}_host_filtered_count.txt"
  '''
}


/*
 Use bwa to align host-filtered reads to the viral reference sequence
 output to bam
*/
process bwa_align_to_refseq {
  label 'lowmem_threaded'                                                                

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/bwa:0.7.8--h7132678_7"
  } else {                                                                      
      container "quay.io/biocontainers/bwa:0.7.8--h7132678_7"
  }      

  input:
  // the filter{size()} functionality here checks if fastq file is empty,
  // which causes bwa to crash
  // see: https://stackoverflow.com/questions/47401518/nextflow-is-an-input-file-empty
  // this means that fastq that are empty at this stage will just stop progress through pipeline
  tuple val(sample_id), path(input_fastq) from post_host_ch_variants.filter{ it[1].getAt(0).size() > 0 }

  output:
  tuple val(sample_id), path("${sample_id}.unsorted.bam") into post_bwa_align_sort_ch


  // in the following 'shell' code block, 
  // !{} will be expanded with the values of variables in the nextflow context, and
  // ${} will be left alone to be used as bash variables (nextflow doesn't expand: leaves for bash)
  //
  // see: https://www.nextflow.io/docs/latest/process.html#shell
  shell:
  '''
  #
  # this complicated bit constructs a RG header needed by GATK in downstream bam-processing steps
  #
  # see: https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
  #
  # this solution from: https://www.biostars.org/p/280837/
  #
  # the first fastq header from the file
  header=$(head -n 1 !{input_fastq[0]}) 
  id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')
  sm=$(echo $header | head -n 1 | grep -Eo "[ATGCN]+$")
  rg=$(echo "@RG\\tID:$id\\tSM:$id"_"$sm\\tLB:$id"_"$sm\\tPL:ILLUMINA")
  
  bwa mem \
  -L !{params.bwa_clipping_penalty} \
  -t !{params.refseq_bwa_threads} \
  -R $rg \
  !{params.refseq_fasta} !{input_fastq}  > !{sample_id}.unsorted.bam
  '''
}

/*
 Sort bwa-aligned bam
*/
process sort_refseq_aligned_bam {
  label 'lowmem_threaded'                                                                

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0"
  } else {                                                                      
      container "quay.io/biocontainers/samtools:1.14--hb421002_0"               
  }     

  input:
  tuple val(sample_id), path(unsorted_bam) from post_bwa_align_sort_ch

  output:
  tuple val(sample_id), path("${sample_id}.bam") into post_bwa_align_ch

  script:
  """
  samtools sort -@${task.cpus} -o ${sample_id}.bam $unsorted_bam
  """
}

/*
 Use iVar to soft-trim primer sites
*/
process ivar_trim {
  label 'lowmem_non_threaded'                                                                
  publishDir "${params.bam_out_dir}", mode:'link', pattern: "*.bam"             

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/ivar:1.3.1--hecb563c_3"
  } else {                                                                      
      container "quay.io/biocontainers/ivar:1.3.1--hecb563c_3"
  }      

  input:
  tuple val(sample_id), path(bam) from post_bwa_align_ch

  output:
  tuple val(sample_id), path("${sample_id}.ivar_trim_unsorted.bam") into post_ivar_trim_ch

  script:

  // should we do ivar-trimming or not?
  if (params.ivar_trim) {

  // do ivar trimming - just pass bam through
  // Usage: ivar trim -i <input.bam> -b <primers.bed> -p <prefix> [-m <min-length>] [-q <min-quality>] [-s <sliding-window-width>] 
  """
  ivar trim \
    -i $bam \
    -x $params.ivar_trim_offset \
    -e \
    -b $params.primer_bed \
    -p "${sample_id}.ivar_trim_unsorted" \
  """

  } else {

  // don't do ivar trimming - just pass bam through

  """
  ln $bam ${sample_id}.ivar_trim_unsorted.bam
  """
  }
}

/*
 Sort ivar-trimmed bam (ivar doesn't output sorted bam)
*/
process sort_ivar_trimmed_bam {
  label 'lowmem_threaded'                                                                

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0"
  } else {                                                                      
      container "quay.io/biocontainers/samtools:1.14--hb421002_0"               
  }     

  input:
  tuple val(sample_id), path(unsorted_bam) from post_ivar_trim_ch

  output:
  tuple val(sample_id), path("${sample_id}.ivar_trim.bam") into post_ivar_sort_ch

  script:
  """
  samtools sort -@${task.cpus} -o ${sample_id}.ivar_trim.bam $unsorted_bam
  """
}


/* 
  Count the # of mapped reads in the bam file and output to channel
  This is to avoid gatk crashing in case of empty bam files
*/
process count_bam_records {
  label 'lowmem_non_threaded'                                                                

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0"
  } else {                                                                      
      container "quay.io/biocontainers/samtools:1.14--hb421002_0"               
  }     

  input:
  tuple val(sample_id), path(input_bam) from post_ivar_sort_ch

  output:
  tuple env(bam_records), val(sample_id), path(input_bam) into post_bwa_count_ch

  """
  # count # of aligned reads: this assumes bam only contains mapped reads...
  bam_records=`samtools view -S -F 0x4 $input_bam | grep -v ^@ | wc -l`
  """

}

/*
 This process performs gatk "Base Quality Score Recalibration" (BQSR).
 BQSR is "a data pre-processing step that detects systematic errors made by 
 the sequencing machine when it estimates the accuracy of each base call."
 see: 
 https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-

 The use of this tool is recommended for downstream variant calling using lofreq
 see: https://csb5.github.io/lofreq/commands/
*/
process apply_bsqr {
  label 'lowmem_non_threaded'                                                                
  publishDir "${params.bam_out_dir}", mode:'link', pattern: "*.bam"             

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/gatk4:4.2.6.0--hdfd78af_0"
  } else {                                                                      
      container "quay.io/biocontainers/gatk4:4.2.6.0--hdfd78af_0"
  }      

  input:
  // exclude bam with no mapped reads
  tuple env(bam_records), val(sample_id), path(input_bam) from post_bwa_count_ch.filter { Integer.parseInt(it[0]) > 0 }

  output:
  tuple val(sample_id), path("${sample_id}.${params.refseq_name}.bam") into post_bsqr_ivar_ch

  output:
  tuple val(sample_id), path("${sample_id}.${params.refseq_name}.bam") into post_bsqr_snv_ch
  tuple val(sample_id), path("${sample_id}.${params.refseq_name}.bam") into post_bsqr_depth_ch
  tuple val(sample_id), path("${sample_id}.${params.refseq_name}.bam") into post_bsqr_indel_ch
  tuple val(sample_id), path("${sample_id}.${params.refseq_name}.bam") into post_bsqr_count_ch
  tuple val(sample_id), path("${sample_id}.${params.refseq_name}.bam") into post_bsqr_consensus_ch
  tuple val(sample_id), path("${sample_id}.${params.refseq_name}.bam") into post_bsqr_stats_ch

  """
  gatk BaseRecalibrator \
     -I ${input_bam} \
     -R ${params.refseq_fasta} \
     --known-sites ${params.ignore_regions} \
     -O recalibration.table

  gatk ApplyBQSR \
     -R ${params.refseq_fasta} \
     -I ${input_bam} \
     --bqsr-recal-file recalibration.table \
     -O ${sample_id}.${params.refseq_name}.bam
  """
}


/*
 Count # of reads aligned to refseq
*/
process refseq_aligned_read_count {
  label 'lowmem_non_threaded'                                                                
  publishDir "${params.counts_out_dir}", mode:'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0"
  } else {                                                                      
      container "quay.io/biocontainers/samtools:1.14--hb421002_0"               
  }     

  input:
  tuple val(sample_id), path(bam) from post_bsqr_count_ch

  output:
  path("${sample_id}_refseq_aligned_count.txt") into post_count_refseq_aligned_ch

  shell:
  // See: http://www.metagenomics.wiki/tools/samtools/number-of-reads-in-bam-file
  // for an explanation of what -F 2436 means in the samtools command.  
  //
  // 2436 is the sum of the following bitwise flags:
  // any alignment with any of these will be excluded:
  // 
  //  4     unmapped read
  //  128   second in pair 
  //  256   not primary alignment
  //  2048  supplementary alignment
  //  ----
  //  2436  sum of above
  //
  // Basically, we are counting reads that are mapped to the viral refseq
  // We are *not* counting paired reads (aligned R2) because previous fastq counts 
  // only count R1 reads , so for all paired-read data, counts are reported 
  // consistetnely in terms of read pairs and not total reads.  
  //
  '''
  samtools view -S !{bam} -F 2436 | wc -l | awk '{print "!{sample_id}" "\trefseq_aligned\t" $1}' > "!{sample_id}_refseq_aligned_count.txt"
  '''
}

/*
 Calculate some mapping statistics from the bam file
*/
process tabulate_mapping_stats_one {
  label 'lowmem_non_threaded'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0"
  } else {                                                                      
      container "quay.io/biocontainers/samtools:1.14--hb421002_0"               
  }     

  when:
  params.plot_mapping_stats

  input:
  tuple val(sample_id), path(input_bam) from post_bsqr_stats_ch

  output:
  path("${input_bam}.mapping_stats") optional true into post_stats_ch

  shell:
  '''
  # process output using samtools
  samtools stats -i 2000 !{input_bam} | grep ^IS | cut -f 2- | awk '{ print "!{sample_id}" "\t" $0; }' > !{input_bam}.mapping_stats
  '''
}

/*
 This process concatenates all the mapping stats files 
*/
process tabulate_stats {
  publishDir "${params.outdir}", mode:'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity') {                              
      container "library://stenglein-lab/r_variant_tools/r_variant_tools:1.0.0"
  } 

  input: 
  path(mapping_stats_files) from post_stats_ch.collect()

  output: 
  path("mapping_stats_plot.pdf") 
  path("all.mapping_stats") 

  script:
  """
  cat $mapping_stats_files > all.mapping_stats
  Rscript ${params.script_dir}/plot_mapping_stats.R ${params.script_dir} all.mapping_stats
  """
}



/*
 Tabulate depth of coverage over refseq

 This prepends coverage depth info for all samples with the 
 sample ID as a first column so it'll be tidy format for import into 
 R and processing with tidyverse packages
 
*/
process tabulate_depth_one {
  label 'lowmem_non_threaded'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0"
  } else {                                                                      
      container "quay.io/biocontainers/samtools:1.14--hb421002_0"               
  }     

  input:
  tuple val(sample_id), path(input_bam) from post_bsqr_depth_ch

  output:
  path("${input_bam}.depth") optional true into post_depth_ch

  shell:
  '''
  # process output using samtools
  samtools depth -a -d 0 !{input_bam} | awk '{ print "!{sample_id}" "\t" $0; }'  > !{input_bam}.depth
  '''
}

/*
 This process concatenates all the depth files 
 into a single file using the collectFile operator
*/
process tabulate_depth {
  publishDir "${params.outdir}", mode:'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity') {                              
      container "library://stenglein-lab/r_variant_tools/r_variant_tools:1.0.0"
  } 

  input: 
  path(depth_files) from post_depth_ch.collect()

  output: 
  path("all.depth") into tabulate_dvg_depth_ch
  path("all.depth") into analyze_variants_depth_ch
  path("coverage_plot.pdf") 
  path("Average_depths.xlsx") into average_depth_ch


  script:
  """
  cat $depth_files > all.depth
  Rscript ${params.script_dir}/plot_depth.R ${params.script_dir} all.depth
  """
}

/*
 Call DIs using DI-tector
*/
process call_dvgs {
  label 'lowmem_threaded'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/python:3.10"
  } else {                                                                      
      container "quay.io/biocontainers/python:3.10"
  }     

  input:
  tuple val(sample_id), path(input_fastq) from post_host_ch_dvg

  // this will skip execution of dvg calling with di-tector unless
  // this param is set to true
  // this will also cause downstream processes to be skipped
  when:
  params.call_variants && params.run_ditector

  output:
  // TODO: DI-tector fails when no reads: check for this?
  tuple val(sample_id), path("DI_counts.txt") optional true into post_dvg_call_ch

  script:
  // TODO: parameterize
  // -x threads
  // -p polarity (+/- sense: for calling something 5' or 3' SB
  // -n number of supporting reads necessary
  // -l minimum length of dvg to report

  // handle single end or paired data
  def r1 = input_fastq[0]
  def r2 = input_fastq[1] ? input_fastq[1] : ""

  """
  # combined R1 and R2 if necessary
  cat $r1 $r2 > ${sample_id}_R12_fh.fastq

  # run DI-tector
  python3 ${params.ditector_script} ${params.refseq_fasta} ${sample_id}_R12_fh.fastq -o "." -t DI -x ${task.cpus} -p 0 -n 4 -l 0
  """
}

process process_dvg_calls {
  label 'lowmem_non_threaded'
  publishDir "${params.ditector_out_dir}", mode:'link'

  input:
  tuple val(sample_id), path(di_counts) from post_dvg_call_ch

  output:
  path("${sample_id}_di_counts.txt") into post_dvg_process_ch

  // in the following 'shell' code block, 
  // !{} will be expanded with the values of variables in the nextflow context, and
  // ${} will be left alone to be used as bash variables (nextflow doesn't expand: leaves for bash)
  //
  // see: https://www.nextflow.io/docs/latest/process.html#shell
  shell:
  '''
  cat !{di_counts} | grep 'DVG' | grep -v -e "^=" -e "DVG's" | awk '{ print "!{sample_id}" "\t" $0; }'  > "!{sample_id}_di_counts.txt"
  '''
}

process tabulate_dvg_calls {
  label 'lowmem_non_threaded'
  publishDir "${params.outdir}", mode:'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity') {                              
      container "library://stenglein-lab/r_variant_tools/r_variant_tools:1.0.0"
  } 

  input:
  path(all_depth) from tabulate_dvg_depth_ch
  path(di_count_files) from post_dvg_process_ch.collect()

  output:
  path("dvg_summary.xlsx") optional true

  script:
  """
  Rscript ${params.script_dir}/process_ditector_output.R ${params.script_dir} \
    $all_depth \
    $di_count_files
  """
}

/*
 Call SNVs using lofreq
*/
process call_snvs {
  label 'lowmem_threaded'
  publishDir "${params.vcf_out_dir}", mode:'link'                               

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py39h43839c7_5"
  } else {                                                                      
      container "quay.io/biocontainers/lofreq:2.1.5--py39h43839c7_5"
  }      

  when:
  params.call_variants 

  input:
  tuple val(sample_id), path(input_bam) from post_bsqr_snv_ch

  output:
  tuple val(sample_id), path("${input_bam}.snv.vcf") into post_snv_call_ch

  script:
  """
  # ---------------
  # Variant calling
  # ---------------

  # use lofreq to call variants
  lofreq index ${input_bam}

  # lofreq call to quantify variant frequencies 
  lofreq call-parallel --pp-threads ${task.cpus} --no-default-filter -f ${params.refseq_fasta} -o ${input_bam}.pre_vcf ${input_bam}

  # call lofreq filter separately to avoid doing strand-bias filtering
  # -v N --> requires minimum Nx coverage (e.g. 40 = call variants at positions with > 40x coverage)
  # -V 0 --> no coverage maximum                                                
  # -a N --> call variants above fraction N (e.g. 0.01 = call variants with >1% allele freq)
  # -A 0 --> no maximum allele frequency 
  lofreq filter -v ${params.min_depth_for_variant_call} -V 0 -a ${params.min_allele_freq} -A 0 --no-defaults -i  ${input_bam}.pre_vcf -o ${input_bam}.snv.vcf
  """
}


/*
 Call Indels also using lofreq
*/
process call_indels {
  label 'lowmem_threaded'
  publishDir "${params.vcf_out_dir}", mode:'link'                               

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py39h43839c7_5"
  } else {                                                                      
      container "quay.io/biocontainers/lofreq:2.1.5--py39h43839c7_5"
  }      

  when:
  params.call_variants

  input:
  tuple val(sample_id), path(input_bam) from post_bsqr_indel_ch

  output:
  tuple val(sample_id), path("${input_bam}.indel.vcf") into post_indel_call_ch

  shell:
  '''
  lofreq indelqual --dindel -f !{params.refseq_fasta} -o !{input_bam}.indelqual.bam !{input_bam}

  lofreq index !{input_bam}.indelqual.bam
                                                                                
  lofreq call-parallel --pp-threads !{task.cpus} --no-default-filter --call-indels --only-indels -f !{params.refseq_fasta} !{input_bam}.indelqual.bam -o !{input_bam}.indel.pre_vcf
                                                                                
  lofreq filter -v !{params.min_depth_for_variant_call} -V 0 -a !{params.min_allele_freq} -A 0 --no-defaults -i  !{input_bam}.indel.pre_vcf -o !{input_bam}.indel.vcf
  '''
}


/* 
 Call consensus sequences for each dataset 

 Here we are using mpileup for variant calling instead of lofreq 
 because I couldn't get lofreq vcf to work as input for bcftools

*/
process call_dataset_consensus_samtools {
  // publishDir "${params.consensus_out_dir}", mode:'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0"
  } else {                                                                      
      container "quay.io/biocontainers/samtools:1.14--hb421002_0"               
  }     

  input:
  tuple val(sample_id), path(input_bam) from post_bsqr_consensus_ch

  output:
  tuple val(sample_id), path("${sample_id}.consensus.bam") into consensus_samtools_ch

  script:
  """
  # consensus calling according to the strategy outlined here: https://www.biostars.org/p/367626/
  samtools sort $input_bam -o sorted_input.bam
  samtools mpileup -uf ${params.refseq_fasta} sorted_input.bam > ${sample_id}.consensus.bam
  """
}

/* 
 Call consensus sequences step 2:  bcftools 
*/
process call_dataset_consensus_bcftools {

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/bcftools:1.15.1--h0ea216a_0"
  } else {                                                                      
      container "quay.io/biocontainers/bcftools:1.15.1--h0ea216a_0"
  }     

  input:
  tuple val(sample_id), path(consensus_bam) from consensus_samtools_ch

  output:
  tuple val(sample_id), path("${sample_id}_consensus.fastq") into consensus_fastq_ch

  script:
  """
  # consensus calling according to the strategy outlined here: https://www.biostars.org/p/367626/
  bcftools call -c $consensus_bam | vcfutils.pl vcf2fq > ${sample_id}_consensus.fastq
  """
}


/* 
 Call consensus sequences for each dataset step 3: seqtk.  

 Use seqtk to convert fastq generated by samtools/bcftools -> a consensus fasta
*/
process call_dataset_consensus {
  // publishDir "${params.consensus_out_dir}", mode:'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3"
  } else {                                                                      
      container "quay.io/biocontainers/seqtk:1.3--h5bf99c6_3"                   
  }       

  input:
  tuple val(sample_id), path(consensus_fastq) from consensus_fastq_ch

  output:
  tuple val(sample_id), path("${sample_id}_consensus.fasta") into consensus_completeness_fasta_ch
  tuple val(sample_id), path("${sample_id}_consensus.fasta") into assign_pangolin_ch

  script:
  """
  # Convert .fastq to .fasta and set bases of quality lower than 20 to N
  # the sed here adds in sample ID to the fasta header line, otherwise they all just have the same name 
  seqtk seq -aQ33 -q20 -n N $consensus_fastq | sed s/">.*"/">${sample_id}"/ > ${sample_id}_consensus.fasta
  """
}


/*
   Calculate the fraction of non-N bases (fraction completeness) of consensus sequence(s) for a dataset
*/
process calculate_consensus_completeness {

  input:
  tuple val(sample_id), path(consensus_fasta) from consensus_completeness_fasta_ch

  output:
  path("*consensus_completeness.txt") into individual_consensus_completeness_ch
  tuple path("*fraction_complete.txt"), path(consensus_fasta) into triage_consensus_ch
  tuple path("*fraction_complete.txt"), path(consensus_fasta) into triage_consensus_fail_ch

  shell:
  '''
  !{params.script_dir}/determine_consensus_completeness.pl !{consensus_fasta} > !{sample_id}_fraction_complete.txt
  cat !{sample_id}_fraction_complete.txt | awk '{ print "!{sample_id}" "\t" $0; }'  > "!{sample_id}_consensus_completeness.txt"
  '''
}


/*
   Triage consensus sequences based on whether they satisfy a particular pass/fail criteria
   For now, must be > 95% complete (95% of bases have to be not Ns)
   Set this cutoff using params.minimum_fraction_called
*/
triage_consensus_ch
  .filter{ Float.parseFloat(it[0].text) >= params.minimum_fraction_called }
  .map { it[1] }
  .into { consensus_fasta_individual_ch ; consensus_fasta_sufficient_ch }

triage_consensus_fail_ch
  .filter{ Float.parseFloat(it[0].text) < params.minimum_fraction_called }
  .map { it[1] }
  .set { consensus_fasta_insufficient_ch }


/*
  Output fasta with sufficient coverage

  This simply uses publishDir to put a copy of the fasta in a designated output directory
*/
process output_fasta_with_sufficient_coverage {
  publishDir "${params.consensus_out_dir}/${params.consensus_pass_dir}", mode:'link'

  input:
  path (consensus_fasta) from consensus_fasta_sufficient_ch

  output:
  path (consensus_fasta) into sufficient_cov_ch

  script:
  """
  touch $consensus_fasta 
  """
}

/*
  Output fasta with insufficient coverage

  This simply uses publishDir to put a copy of the fasta in a designated output directory
*/
process output_fasta_with_insufficient_coverage {
  publishDir "${params.consensus_out_dir}/${params.consensus_fail_dir}", mode:'link'

  input:
  path (consensus_fasta) from consensus_fasta_insufficient_ch

  output:
  path (consensus_fasta) into insufficient_cov_ch

  script:
  """
  touch $consensus_fasta 
  """
}

/*
  This process optionally obfuscate sample IDs for submission to public databases
*/
process obfuscate_sample_IDs {
  publishDir "${params.outdir}", mode:'link'                               

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity') {                              
      container "library://stenglein-lab/r_variant_tools/r_variant_tools:1.0.0"
  } 

  input:
  // the map{it[0]} here pulls out the first element (the sample IDs)
  // from the tuple.  I.e. it ignores fastq files (2nd element of the tuple).
  // collect merges all these samples_ids into a single list
  val (sample_ids) from sample_ids_ch.map{it[0]}.collect()

  output:
  path("obfuscated_sample_ids.txt") into obfuscated_ids_report_ch
  path("obfuscated_sample_ids.txt") into obfuscated_ids_gisaid_ch

  script:
  def sample_ids_text = sample_ids.join(" ")
  """
  Rscript ${params.script_dir}/obfuscate_sample_ids.R ${params.key_file} $sample_ids_text 
  """
}


/*
   Concatenate the consensus completeness info into a single file to be fed to a reporting script
*/
process tabulate_consensus_completeness {
  publishDir "${params.outdir}", mode:'link'                               

  input:
  path(files) from individual_consensus_completeness_ch.collect()

  output:
  path("all_consensus_completeness.txt") into consensus_completeness_ch

  script:
  """
  # mash the output together
  cat $files > all_consensus_completeness.txt
  """
}

/*
   Update pangolin info for up-to-date lineage assignment
*/

process update_pangolin_info {
  
  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/pangolin:4.0.4--pyhdfd78af_0"
  } else {                                                                      
      container "quay.io/biocontainers/pangolin:4.0.4--pyhdfd78af_0"
  }       

  output: 
  val ("pangolin_update_complete") into pangolin_update_ch

  script:
  """
  mkdir -p $params.pangolin_datadir
  pangolin --update-data --datadir $params.pangolin_datadir
  """
}

/*
   Assign consensus sequence to a PANGO lineage using Pangolin

   Might be better to do this on a single concatenated file
   containing all the sequences...
*/

process assign_to_pango_lineage {
  publishDir "${params.pangolin_out_dir}", mode:'link'
  
  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/pangolin:4.0.4--pyhdfd78af_0"
  } else {                                                                      
      container "quay.io/biocontainers/pangolin:4.0.4--pyhdfd78af_0"
  }       

  input: 
  tuple val(sample_id), path(consensus_fasta) from assign_pangolin_ch
  val (pangolin_updated_flag) from pangolin_update_ch

  output: 
  path ("*lineage_report.csv") into collect_pangolin_ch

  script:
  """
  pangolin --datadir $params.pangolin_datadir --outfile ${sample_id}_lineage_report.csv ${consensus_fasta}
  """
}

/*
   Concatenate the pangolin assignment info into a single file to be fed to a reporting script
*/
process tabulate_pangolin_info {
  publishDir "${params.outdir}", mode:'link'                               

  input:
  path(files) from collect_pangolin_ch.collect()

  output:
  path("pangolin_lineage_report.csv") into pangolin_report_ch

  script:
  """
  # pull out header line from first file (assumes word "taxon" in header)
  head -1 $files | grep taxon | head -1 > pangolin_lineage_report.csv

  # pull out the data lines from rest of files 
  # (assumes word "taxon" in header but not data lines)
  cat $files | grep -v taxon >> pangolin_lineage_report.csv
  """
}


/*
  Report on datasets
*/
process report_on_datasets {
  publishDir "${params.outdir}", mode:'link'                               

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity') {                              
      container "library://stenglein-lab/r_variant_tools/r_variant_tools:1.0.0"
  } 

  input:
  path(consensus_completeness) from consensus_completeness_ch
  path(average_depth_xls) from average_depth_ch
  path(pangolin_lineage_report) from pangolin_report_ch
  path(obfuscated_ids) from obfuscated_ids_report_ch

  output:
  path("*dataset_summary.xlsx") into dataset_summary_ch
  path("*.xlsx") 
  path("*.pdf") 

  script:
  """
  Rscript ${params.script_dir}/report_on_datasets.R ${params.script_dir} $average_depth_xls $consensus_completeness $pangolin_lineage_report ${params.minimum_fraction_called} ${params.output_prefix} $obfuscated_ids
  """
}

/* 
 use snpEff to annotate vcfs 
 */
process annotate_variants {
  publishDir "${params.vcf_out_dir}", mode:'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/snpeff:5.1--hdfd78af_1"
  } else {                                                                      
      container "quay.io/biocontainers/snpeff:5.1--hdfd78af_1"
  }      

  input:
  tuple val(sample_id), path(vcf) from post_indel_call_ch.concat(post_snv_call_ch)

  output:
  tuple val(sample_id), path("${vcf}.snp_eff") into post_variant_annotate_ch

  script:

  // setup optional custom annotation using the -interval command line option
  // see: https://pcingola.github.io/SnpEff/se_commandline/
  def custom_annotation  = params.custom_annotations ? "-interval $params.custom_annotations_bed" : ""

  """
  snpEff ann \
    -c ${params.snpeff_cfg} \
    -dataDir ${params.snpeff_data} \
    $custom_annotation \
    -no-downstream \
    -no-upstream \
    -no-intergenic \
    -no-intron \
    ${params.refseq_name} $vcf > ${vcf}.snp_eff
  """
}

/*
 use SnpSift to extract SnpEff annotations
 */
process extract_annotated_variant_fields {
  publishDir "${params.vcf_out_dir}", mode:'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/snpsift:5.1--hdfd78af_0"
  } else {                                                                      
      container "quay.io/biocontainers/snpsift:5.1--hdfd78af_0"
  }      

  input:
  tuple val(sample_id), path(snp_eff) from post_variant_annotate_ch

  output:
  tuple val(sample_id), path("${snp_eff}.snp_sift") into post_snp_sift_ch

  script:
  """
  SnpSift extractFields -e "." -s "," ${snp_eff} CHROM POS REF ALT AF DP SB INDEL ANN[*].EFFECT ANN[*].IMPACT ANN[*].GENE ANN[*].AA_POS ANN[*].HGVS_P ANN[*].FEATUREID > ${snp_eff}.snp_sift
  """
}

/*
 This process prepends SnpSift info for all samples with the 
 sample ID as a first column so it'll be tidy format for import into 
 R and processing with tidyverse packages
*/
process prepend_snp_sift_output {
  label 'lowmem_non_threaded'
  publishDir "${params.vcf_out_dir}", mode:'link'

  input:
  tuple val(sample_id), path(snp_sifts) from post_snp_sift_ch
  .groupTuple() // group tuple here will merge the snv and indel variant files into a single one (using cat below)

  output:
  path("${sample_id}.variants.tsv") optional true into post_prepend_snp_sift_ch

  shell:
  '''
  # grep: remove header lines (contain the text CHROM at beginning)
  # awk: add a new first column w/ sample id
  cat !{snp_sifts} | grep -v '^CHROM' | awk '{ print "!{sample_id}" "\t" $0; }' > "!{sample_id}.variants.tsv"
  '''
}


/*
 tabulate snpeff snv annotations for all datasets using snpsift output

 this will tabulate all SNV and indel variants together
*/
process tabulate_snpeff_variants {
  publishDir "${params.outdir}", mode:'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity') {                              
      container "library://stenglein-lab/r_variant_tools/r_variant_tools:1.0.0"
  } 

  input:
  path(all_depth) from analyze_variants_depth_ch
  path(snp_sifts) from post_prepend_snp_sift_ch.collect()

  output:
  path("variant_summary.xlsx") 
  path("sample_correlation_heatmap.pdf") optional true

  script:
  """
  Rscript ${params.script_dir}/analyze_snpeff_variants.R ${params.script_dir} \
    ${params.min_allele_freq} \
    $all_depth \
    ${params.min_depth_for_variant_call} \
    $snp_sifts
  """
}

process tabulate_fastq_counts {
  publishDir "${params.outdir}", mode: 'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity') {                              
      container "library://stenglein-lab/r_variant_tools/r_variant_tools:1.0.0"
  } 

  input:
  path(all_count_files) from post_count_initial_ch.concat(post_count_trim_ch, post_count_host_ch, post_count_refseq_aligned_ch).collect()

  output:
  path ("*all_read_counts.txt") 
  path ("*summarized_read_counts.txt") 
  path ("*filtering_plots.pdf") 

  script:
  """
  Rscript ${params.script_dir}/process_fastq_counts.R ${params.script_dir} ${params.output_prefix} ${all_count_files} 
  """
}

/*
  Replace the fasta sample IDs with sample IDs in submission metadata
  so that fasta sequence names match names in metadata
*/
process collect_consensus_fasta {
  publishDir "${params.outdir}", mode:'link'                               

  input:
  path(all_fasta) from consensus_fasta_individual_ch.collect()
  
  output:
  path("*consensus_sequences_original_ids.fasta") into consensus_fasta_ch

  script:
  """
  # concatenate all fasta into one file with original sample IDs
  cat $all_fasta > "${params.output_prefix}consensus_sequences_original_ids.fasta"
  """

}

/*
  This process prepares files suitable for uploading to GISAID
*/
process prepare_gisaid_submission_files {
  publishDir "${params.outdir}", mode:'link'                               
  
  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity') {                              
      container "library://stenglein-lab/r_variant_tools/r_variant_tools:1.0.0"
  } 

  when:
  params.prepare_gisaid

  input:
  path(all_fasta) from consensus_fasta_ch
  path(dataset_summary) from dataset_summary_ch
  path(obfuscated_ids) from obfuscated_ids_gisaid_ch

  output:
  path ("*.fasta") into gisaid_fasta_ch
  path ("*.xlsx") into gisaid_xls_ch

  script:
  """
  Rscript ${params.script_dir}/prepare_gisaid_submission_files.R ${params.script_dir} $all_fasta $dataset_summary $obfuscated_ids ${params.minimum_fraction_called} ${params.gisaid_metadata_file} ${params.output_prefix}
  """
}
