#!/usr/bin/env nextflow

/*
    Stenglein lab variant analysis nextflow pipeline 

    October 22, 2020 

    Mark Stenglein
*/


// TODO: command line options and parameter checking


params.fastq_dir = "$baseDir/fastq/"
params.initial_fastqc_dir = "$baseDir/initial_fastqc/" 
params.post_trim_fastqc_dir = "$baseDir/post_trim_fastqc/" 
params.outdir = "results"                                                       


params.maxForks_lowmem = "16"
params.maxForks_highmem = "4"

// ------------------
// Trimming settings
// ------------------
// SARS-CoV-2 ARTIC v3 primers are 22-30 bp long, so always trim 30
// TODO: don't do this for non-amplicon datasets! 
params.always_trim_5p_bases = "30" 
params.always_trim_3p_bases = "1" 
params.post_trim_min_length = "60" 

// --------------------
// Host cell filtering
// --------------------
// Vero cell: African green monkey genome for host filtering
params.host_bt_index = "/home/databases/primates/agm_genome"
params.host_bt_suffix = "agm_genome"
params.host_bt_min_score = "60"
params.host_bt_threads = "8"


// SARS-CoV-2 wa1 refseq
params.viral_refseq_name = "wa1"
params.viral_fasta = "${baseDir}/viral_refseq/${params.viral_refseq_name}.fasta"
params.viral_gb = "${baseDir}/viral_refseq/${params.viral_refseq_name}.gb"
params.viral_gff = "${baseDir}/viral_refseq/${params.viral_refseq_name}.gff"
params.viral_bt_index = "${baseDir}/viral_refseq/${params.viral_refseq_name}"
params.viral_bt_min_score = "120"

params.viral_bwa_threads = "8"


// where are R scripts found...
params.R_bindir="${baseDir}"

// conda for snpEFF
params.snpeff_jar = "/home/apps/snpEff/snpEff.jar" 
params.snpeff_cfg = "${baseDir}/snpEff.config" 
params.snpeff_data = "${baseDir}/viral_refseq/snpeff_data/"

// conda or something for gatk?
params.gatk_exe="/home/apps/gatk-4.1.9.0/gatk"

// ignoring regions with 
params.ignore_regions="${baseDir}/viral_refseq/ignore_regions.bed"

// flag to optionally run cd-hit to collapse non-unique reads
// skip this step by default
params.skip_collapse_to_unique = true
// cd-hit-dup cutoff for collapsing reads with >= this much fractional similarity
params.duplicate_cutoff = "0.98"


params.min_depth_for_variant_call="50"
params.min_allele_freq="0.03"


// TODO: command line arg processing and validating 

// TODO: check that appropriate refseq files exist (fasta, gb, etc.)

// TODO: handle single end or paired end, possibly automatically

// TODO: better control of concurrency in terms of processes

// TODO: move scripts to a bin subdir

// TODO: put all up on github

// TODO: conda

// TODO: BSQR fails in case of no mapping reads... deal with this possibility 

/*
 These fastq files represent the main input to this workflow
*/
// TODO: autodetect single end vs paired end
Channel
    .fromFilePairs("${params.fastq_dir}/*_R{1,2}*.fastq", size: -1, checkIfExists: true, maxDepth: 1)
    .into {samples_ch_qc; samples_ch_trim}



/*
   Setup some initial indexes and dictionaries needed by downstream processes.
   Only do this once at beginning.
   TODO: remove initial .dict file
*/
process setup_indexes {

  output:
  val("indexes_complete") into post_index_setup_ch

  script:
  """
  # ------------------
  # lofreq fasta index
  # ------------------
  lofreq faidx ${params.viral_fasta}

  # ----------------
  # setup snpEFF db
  # ----------------

  # make a minimal snpEff.config 
  rm -f snpEff.config
  echo "${params.viral_refseq_name}.genome: ${params.viral_refseq_name}" > ${params.snpeff_cfg}

  # make a directory for the snp eff db
  mkdir -p ${params.snpeff_data}/${params.viral_refseq_name}

  # cp fasta and gb for virus refseq to the directory location snpeff is expecting
  # these must exist! 
  cp ${params.viral_fasta} ${params.snpeff_data}/${params.viral_refseq_name}/sequences.fa
  cp ${params.viral_gb} ${params.snpeff_data}/${params.viral_refseq_name}/genes.gbk

  # build the snpEff db
  java -jar ${params.snpeff_jar} build -c ${params.snpeff_cfg} -nodownload -v -genbank -dataDir ${params.snpeff_data} ${params.viral_refseq_name}

  # -----------------------
  # bwa index viral refseq
  # -----------------------
  bwa index ${params.viral_fasta}

  # -----------------
  # GATK index setup
  # -----------------
  # setup gatk indexes for BSQR
  ${params.gatk_exe} IndexFeatureFile -I ${params.ignore_regions} 

  rm -f "${baseDir}/viral_refseq/${params.viral_refseq_name}.dict"
  ${params.gatk_exe} CreateSequenceDictionary -R ${params.viral_fasta}
  """
}

/*
 Run fastqc on input fastq 
*/
process initial_qc {
  maxForks params.maxForks_lowmem

  input:
  tuple val(sample_id), path(initial_fastq) from samples_ch_qc

  output:
  val(sample_id) into post_initial_qc_ch
  val(sample_id) into write_datasets_ch
  // TODO: count

  script:
  """
  mkdir -p  ${params.initial_fastqc_dir} 
  fastqc -o ${params.initial_fastqc_dir} $initial_fastq 
  """
}

/* 
Collect and compress all raw fastq files --> deliverables
*/
// TODO: 

/*
 Use multiqc to merge initial fastqc reports
*/
process initial_multiqc {
  publishDir "${params.outdir}", mode:'link'

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
 Write all the sample IDs to a file
*/
/*
// TODO: fix
process write_ids {

  input:
  val(sample_ids) from write_datasets_ch.collect()

  exec:
  filename="${params.outdir}/dataset_ids.txt"
  // delete if exists
  new File(filename).delete()  
  File file = new File(filename)
  sample_ids.each {
    file.append("${it}\n")
  }
  println file.text
}
*/

/*
 Use cutadapt to trim off adapters and low quality bases
*/
process trim_adapters_and_low_quality {
  // publishDir "${params.outdir}"
  maxForks params.maxForks_lowmem

  input:
  tuple val(sample_id), path(initial_fastq) from samples_ch_trim

  output:
  // TODO: count
  // TODO: multiqc trimming report
  tuple val(sample_id), path("*_f.fastq") optional true into post_trim_qc_ch
  tuple val(sample_id), path("*_f.fastq") optional true into post_trim_ch

  // TODO: put adapters as a param
  script:

  // this handles paired-end data, in which case must specify a paired output file
  def paired_output   = initial_fastq[1] ? "-p ${sample_id}_R2_f.fastq" : ""
  def paired_adapters = initial_fastq[1] ? "-A AGATCGGAAGAGC -G GCTCTTCCGATCT -A AGATGTGTATAAGAGACAG -G CTGTCTCTTATACACATCT" : ""
  // TODO: don't trim this much for non-amplicon data!
  def paired_trimming = initial_fastq[1] ? "-U $params.always_trim_5p_bases -U -${params.always_trim_3p_bases}" : ""

  """
  cutadapt \
   -a AGATCGGAAGAGC -g GCTCTTCCGATCT -a AGATGTGTATAAGAGACAG -g CTGTCTCTTATACACATCT \
   $paired_adapters \
   -q 30,30 \
   --minimum-length ${params.post_trim_min_length} \
   -u ${params.always_trim_5p_bases} \
   -u -${params.always_trim_3p_bases} \
   $paired_trimming \
   -o ${sample_id}_R1_f.fastq \
   $paired_output \
   $initial_fastq 
  """

}

//
// The code below enables optional collapsing of non-unique reads using cd-hit
// 
// For a discussion of conditional execution in nextflow, see
// https://nextflow-io.github.io/patterns/index.html#_conditional_process_executions
//
// (This pattern is from there)
//
// TODO: this is not working...
//
/*
(post_trim_ch, post_collapse_ch) = ( params.skip_collapse_to_unique
                                     ? [Channel.empty(), post_trim_optional_ch]
                                     : [post_trim_optional_ch, Channel.empty()] )
*/

/*
(post_trim_qc_ch, post_collapse_qc_ch) = ( params.skip_collapse_to_unique
                                         ? [Channel.empty(), post_trim_optional_qc_ch]
                                         : [post_trim_optional_qc_ch, Channel.empty()] )
*/

/*
 Use cd-hit to collapse duplicate reads
*/
/*
process collapse_to_unique {
  // publishDir "${params.outdir}"
  maxForks params.maxForks_lowmem

  input:
  // post_trim_ch will be empty if params.skip_collapse_to_unique is true,   
  // in which case, this process will not run...
  tuple val(sample_id), path(r1_fastq), path(r2_fastq) from post_trim_ch

  output:
  // TODO: count
  tuple val(sample_id), path("${sample_id}_R1_fu.fastq"), path("${sample_id}_R2_fu.fastq") into post_collapse_qc_ch
  tuple val(sample_id), path("${sample_id}_R1_fu.fastq"), path("${sample_id}_R2_fu.fastq") into post_collapse_ch

  script:
  """
  cd-hit-dup -i ${r1_fastq} -o2 ${r2_fastq} -o ${sample_id}_R1_fu.fastq -o2 ${sample_id}_R2_fu.fastq -e ${params.duplicate_cutoff} 
  """
}
*/

/*
 Use fastqc to do QC on post-trimmed fastq
*/
process post_trim_qc {
  maxForks params.maxForks_lowmem

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
  publishDir "${params.outdir}"

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
// TODO: switch to bwa for host filtering too?
process host_filtering {
  // publishDir "${params.outdir}", pattern: "*_R1_fh.fastq"
  maxForks params.maxForks_highmem

  input:
  tuple val(sample_id), path(input_fastq) from post_trim_ch

  output:
  tuple val(sample_id), path("*_fh.fastq") optional true into post_host_ch_variants
  // TODO: count
  // TODO: multiqc analysis of bowtie output (host filtering)

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
  -p "${params.host_bt_threads}" \
  $bowtie_file_output 2> "${sample_id}.host_filtering_bt.log" > /dev/null 
  """
}

/* 
Collect and compress all host-filtered fastq files --> deliverables
*/
// TODO: 



/*
 Use bwa to align host-filtered reads to the viral reference sequence
 output to bam
*/
process bwa_align_to_refseq {
  maxForks params.maxForks_lowmem

  input:
  tuple val(sample_id), path(input_fastq) from post_host_ch_variants

  output:
  tuple val(sample_id), path("${sample_id}.bam") into post_bwa_align_ch
  tuple val(sample_id), path("${sample_id}.bam") into post_bwa_align_depth_ch


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
  -t !{params.viral_bwa_threads} \
  -R $rg \
  !{params.viral_fasta} !{input_fastq} | samtools sort -@12 -o !{sample_id}.bam  -
  '''
}


/*
 This process uses gatk "Base Quality Score Recalibration" (BQSR)
 BQSR is "a data pre-processing step that detects systematic errors made by the sequencing machine when it estimates the accuracy of each base call."
 see: https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-

 The use of this tool is recommended for downstream variant calling using lofreq
 see: https://csb5.github.io/lofreq/commands/
*/
process apply_bsqr {
  maxForks params.maxForks_lowmem
  publishDir "${params.outdir}", pattern: "*.bam"

  input:
  tuple val(sample_id), path(input_bam) from post_bwa_align_ch
  // this input value will prevent this process from running before the 
  // appropriate indexes are setup by the setup_indexes process 
  val("indexes_complete") from post_index_setup_ch

  output:
  tuple val(sample_id), path("${sample_id}.${params.viral_refseq_name}.bam") into post_bsqr_snv_ch
  tuple val(sample_id), path("${sample_id}.${params.viral_refseq_name}.bam") into post_bsqr_depth_ch
  tuple val(sample_id), path("${sample_id}.${params.viral_refseq_name}.bam") into post_bsqr_indel_ch

  """
  ${params.gatk_exe} BaseRecalibrator \
     -I ${input_bam} \
     -R ${params.viral_fasta} \
     --known-sites ${params.ignore_regions} \
     -O recalibration.table

  ${params.gatk_exe} ApplyBQSR \
     -R ${params.viral_fasta} \
     -I ${input_bam} \
     --bqsr-recal-file recalibration.table \
     -O ${sample_id}.${params.viral_refseq_name}.bam
  """
}

/*
 Tabulate depth of coverage over refseq
*/
process tabulate_depth_one {
  maxForks params.maxForks_lowmem

  input:
  tuple val(sample_id), path(input_bam) from post_bsqr_depth_ch

  output:
  tuple val(sample_id), path("${input_bam}.depth") optional true into post_depth_ch

  script:
  """
  # process output using samtools
  samtools depth -d 0 ${input_bam} > ${input_bam}.depth
  """
}


/*
 Call SNVs using lofreq
*/
process call_snvs {
  maxForks params.maxForks_highmem
  publishDir "${params.outdir}", mode:'link'

  input:
  tuple val(sample_id), path(input_bam) from post_bsqr_snv_ch

  output:
  // path("${input_bam}.snv.vcf") into post_snv_call_ch
  path("${sample_id}.variant_alleles.txt") into post_variant_call_ch
  tuple val(sample_id), path(depth) optional true into post_variant_call_depth_ch

  script:
  """
  # ---------------
  # Variant calling
  # ---------------

  # use lofreq to call variants
  
  # lofreq call to quantify variant frequencies 
  # lofreq call-parallel --no-default-filter --pp-threads ${params.host_bt_threads} -f ${params.viral_fasta} -o pre_vcf ${input_bam}
  lofreq call --no-default-filter -f ${params.viral_fasta} -o pre_vcf ${input_bam}

  # call lofreq filter separately to avoid doing strand-bias filtering
  # -v 40 --> requires minimum 40x coverage
  # -V 0 --> no coverage maximum
  # -a 0.01 --> call variants above 1% (0.01 = min allele frequency)
  # -A 0 --> no maximum allele frequency
  lofreq filter -v ${params.min_depth_for_variant_call} -V 0 -a ${params.min_allele_freq} -A 0 --no-defaults -i  pre_vcf -o vcf
  
  # ----------------
  # variant analysis
  # ----------------
  
  # analyze variants will determine whether these are non synonymous or synonymous variants
  # and identify the impacted CDS, etc.
  ${baseDir}/analyze_variants ${params.viral_gff} vcf >  "${sample_id}.variant_alleles.txt"
  
  """
}

/*
 tabulate snv calls for all datasets
*/
/*
process tabulate_snvs {
  publishDir "${params.outdir}", mode:'link'

  input:
  path(vcf) from post_snv_call_ch.collect()

  output:
  path("Single_nucleotide_variant_summary.xlsx") into post_snv_variant_tabulate_ch

  script:
  """
  Rscript ${baseDir}/analyze_snv_vcf.R $vcf
  """
}
*/


/*
 Call Indels also using lofreq
*/
process call_indels {
  maxForks params.maxForks_highmem
  publishDir "${params.outdir}", mode:'link'

  input:
  tuple val(sample_id), path(input_bam) from post_bsqr_indel_ch

  output:
  path("${input_bam}.indel.vcf") into post_indel_call_ch

  shell:
  '''
  lofreq indelqual --dindel -f !{params.viral_fasta} -o !{input_bam}.indelqual.bam !{input_bam}

  lofreq call --call-indels --only-indels -f !{params.viral_fasta} !{input_bam}.indelqual.bam | \
  awk '{ print $0 ";" "!{sample_id}"; }' > !{input_bam}.indel.vcf
  '''
}

/*
 tabulate indel calls for all datasets
*/
// TODO: this failed in the case where there was a single variant in the vcf file...
process tabulate_indel_variants {
  publishDir "${params.outdir}", mode:'link'

  input:
  path(vcf) from post_indel_call_ch.collect()

  output:
  path("Structural_variant_summary.xlsx") into post_indel_variant_tabulate_ch

  script:
  """
  Rscript ${baseDir}/analyze_indel_vcf.R ${params.R_bindir} $vcf
  """
}

/*
 tabulate variants for all datasets
*/
process tabulate_variants {
  publishDir "${params.outdir}", mode:'link'

  input:
  path(variant_alleles) from post_variant_call_ch.collect()

  output:
  path("variant_table.txt") into post_variant_tabulate_ch

  script:
  """
  ${baseDir}/tabulate_variants $variant_alleles > variant_table.txt
  """
}

/*
 This process prepends coverage depth info for all samples with the 
 sample ID as a first colum so it'll be tidy format for import into 
 R and processing with tidyverse packages
*/
process prepend_depth {
  maxForks params.maxForks_lowmem

  input:
  tuple val(sample_id), path(depth) from post_depth_ch

  output:
  path("${sample_id}_prepended_depth") optional true into post_prepend_depth_ch

  shell:
  '''
  cat !{depth} | awk '{ print "!{sample_id}" "\t" $0; }' > "!{sample_id}_prepended_depth"
  '''
}

/*
 This process concatenates all the depth files 
 into a single file using the collectFile operator
*/
process tabulate_depth {
  publishDir "${params.outdir}", mode:'link'

  input: 
  path(depth_files) from post_prepend_depth_ch.collect()

  output: 
  path("all.depth") 
  path("coverage_plot.pdf") 


  script:
  """
  cat $depth_files > all.depth
  Rscript ${baseDir}/plot_depth.R ${params.R_bindir} all.depth
  """
}





/*
process count_fastq {
  publishDir "${params.outdir}", mode: 'link'

  input:
  path(di_counts) from post_di_tabulate_ch

  output:
  path ("fastq_counts.txt") 

  script:
  """
  $basedDir/count_all_fastq > fastq_counts.txt
  """
}
*/
