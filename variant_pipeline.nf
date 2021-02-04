#!/usr/bin/env nextflow

/*
    Stenglein lab variant analysis nextflow pipeline 

    October 22, 2020 

    Mark Stenglein
*/


// TODO: command line options and parameter checking


params.fastq_dir = "$baseDir/fastq/"
params.outdir = "$baseDir/results"                                                       
params.initial_fastqc_dir = "${params.outdir}/initial_fastqc/" 
params.post_trim_fastqc_dir = "${params.outdir}/post_trim_fastqc/" 
params.counts_out_dir = "${params.outdir}/fastq_counts/"                        
params.fastq_out_dir = "${params.outdir}/trimmed_fastq/"                        
params.bam_out_dir = "${params.outdir}/bam/"                                    
params.vcf_out_dir = "${params.outdir}/vcf/"  
params.consensus_out_dir = "${params.outdir}/consensus_sequences/"  
params.ditector_out_dir = "${params.outdir}/ditector/"  

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


// SARS-CoV-2 wa1 refseq, or a reference seq of your choosing
params.refseq_dir = "${baseDir}/refseq/"
params.refseq_name = "NC_045512"
// params.refseq_name = "MN985325"
params.refseq_fasta = "${params.refseq_dir}/${params.refseq_name}.fasta"
params.refseq_genbank = "${params.refseq_dir}/${params.refseq_name}.gb"
params.refseq_bt_index = "${params.refseq_dir}/${params.refseq_name}"
params.refseq_bt_min_score = "120"

params.refseq_bwa_threads = "8"


// where are R scripts found...
params.R_bindir="${baseDir}/scripts"
params.scripts_bindir="${baseDir}/scripts"

// DI-tector info
params.ditector_script="${params.scripts_bindir}/DI-tector_06.py"
params.run_ditector = false

// conda for snpEFF
params.snpeff_cfg = "${params.refseq_dir}/snpEff.config" 
params.snpeff_data = "${params.refseq_dir}/snpeff_data/"
params.snpeff_threads = "8"

// regions to omit from BSQR step
params.ignore_regions="${params.refseq_dir}/ignore_regions.bed"

// optional custom variant annotation
params.custom_annotations = true
params.custom_annotations_bed = "${params.refseq_dir}/voc_positions.bed"

// flag to optionally run cd-hit to collapse non-unique reads
// skip this step by default
params.skip_collapse_to_unique = true
// cd-hit-dup cutoff for collapsing reads with >= this much fractional similarity
params.duplicate_cutoff = "0.98"

// min depth and allele freq for calling variants and indels                    
params.min_depth_for_variant_call="40"
params.min_allele_freq="0.03"



// TODO: command line arg processing and validating 

// TODO: check that appropriate refseq files exist (fasta, gb, etc.)

// TODO: BSQR fails in case of no mapping reads... deal with this possibility 

// TODO: handle fastq.gz compressed files 

/*
 These fastq files represent the main input to this workflow
*/
Channel
    .fromFilePairs("${params.fastq_dir}/*_R{1,2}*.fastq*", size: -1, checkIfExists: true, maxDepth: 1)
    .into {samples_ch_qc; samples_ch_trim; samples_ch_count}


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
  lofreq faidx ${params.refseq_fasta}

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

  # -----------------------
  # bwa index viral refseq
  # -----------------------
  bwa index ${params.refseq_fasta}

  # -----------------
  # GATK index setup
  # -----------------
  # setup gatk indexes for BSQR
  gatk IndexFeatureFile --feature-file ${params.ignore_regions} 

  rm -f "${params.refseq_dir}/${params.refseq_name}.dict"
  gatk CreateSequenceDictionary -R ${params.refseq_fasta}
  """
}

/*
 Run fastqc on input fastq 
*/
process initial_qc {
  label 'lowmem_non_threaded'                                                                

  input:
  tuple val(sample_id), path(initial_fastq) from samples_ch_qc

  output:
  val(sample_id) into post_initial_qc_ch
  // TODO: count

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
  cat !{initial_fastq[0]} | wc -l | xargs bash -c 'echo $(($0 / 4))' | awk '{print "!{sample_id}" "\tinitial\t" $1}' > "!{sample_id}_initial_count.txt"
  '''
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
 Use cutadapt to trim off adapters and low quality bases
*/
process trim_adapters_and_low_quality {
  // publishDir "${params.outdir}"
  label 'lowmem_non_threaded'                                                                

  input:
  tuple val(sample_id), path(initial_fastq) from samples_ch_trim

  output:
  tuple val(sample_id), path("*_f.fastq") optional true into post_trim_qc_ch
  tuple val(sample_id), path("*_f.fastq") optional true into post_trim_ch
  tuple val(sample_id), path("*_f.fastq") optional true into post_trim_count_ch

  // TODO: put adapters as a param (?)
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
  cat !{trimmed_fastq[0]} | wc -l | xargs bash -c 'echo $(($0 / 4))' |  awk '{print "!{sample_id}" "\tpost_trimming\t" $1}' > "!{sample_id}_trimmed_count.txt"
  '''
}

/*
 Use fastqc to do QC on post-trimmed fastq
*/
process post_trim_qc {
  label 'lowmem_non_threaded'

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
  $bowtie_file_output 2> "${sample_id}.host_filtering_bt.log" > /dev/null 
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
  cat !{filtered_fastq[0]} | wc -l | xargs bash -c 'echo $(($0 / 4))' |  awk '{print "!{sample_id}" "\tpost_host_filtered\t" $1}' > "!{sample_id}_host_filtered_count.txt"
  '''
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
  label 'lowmem_threaded'                                                                

  input:
  tuple val(sample_id), path(input_fastq) from post_host_ch_variants

  output:
  tuple val(sample_id), path("${sample_id}.bam") into post_bwa_align_ch


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
  -t !{params.refseq_bwa_threads} \
  -R $rg \
  !{params.refseq_fasta} !{input_fastq} | samtools sort -@!{task.cpus} -o !{sample_id}.bam  -
  '''
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


  input:
  tuple val(sample_id), path(input_bam) from post_bwa_align_ch
  // this input value will prevent this process from running before the 
  // appropriate indexes are setup by the setup_indexes process 

  output:
  tuple val(sample_id), path("${sample_id}.${params.refseq_name}.bam") into post_bsqr_snv_ch
  tuple val(sample_id), path("${sample_id}.${params.refseq_name}.bam") into post_bsqr_depth_ch
  tuple val(sample_id), path("${sample_id}.${params.refseq_name}.bam") into post_bsqr_indel_ch
  tuple val(sample_id), path("${sample_id}.${params.refseq_name}.bam") into post_bsqr_count_ch
  tuple val(sample_id), path("${sample_id}.${params.refseq_name}.bam") into post_bsqr_consensus_ch

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
 Tabulate depth of coverage over refseq
*/
process tabulate_depth_one {
  label 'lowmem_non_threaded'

  input:
  tuple val(sample_id), path(input_bam) from post_bsqr_depth_ch

  output:
  tuple val(sample_id), path("${input_bam}.depth") optional true into post_depth_ch

  script:
  """
  # process output using samtools
  samtools depth -a -d 0 ${input_bam} > ${input_bam}.depth
  """
}

/*
 This process prepends coverage depth info for all samples with the 
 sample ID as a first colum so it'll be tidy format for import into 
 R and processing with tidyverse packages
*/
process prepend_depth {
  label 'lowmem_non_threaded'

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
  path("all.depth") into tabulate_dvg_depth_ch
  path("all.depth") into analyze_variants_depth_ch
  path("coverage_plot.pdf") 


  script:
  """
  cat $depth_files > all.depth
  Rscript ${params.R_bindir}/plot_depth.R ${params.R_bindir} all.depth
  """
}

/*
 Call DIs using DI-tector
*/
process call_dvgs {
  label 'lowmem_threaded'

  input:
  tuple val(sample_id), path(input_fastq) from post_host_ch_dvg

  // this will skip execution of dvg calling with di-tector unless
  // this param is set to true
  // this will also cause downstream processes to be skipped
  when:
  params.run_ditector

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

  input:
  path(all_depth) from tabulate_dvg_depth_ch
  path(di_count_files) from post_dvg_process_ch.collect()

  output:
  path("dvg_summary.xlsx") optional true

  script:
  """
  Rscript ${params.R_bindir}/process_ditector_output.R ${params.R_bindir} \
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


  input:
  tuple val(sample_id), path(input_bam) from post_bsqr_indel_ch

  output:
  tuple val(sample_id), path("${input_bam}.indel.vcf") into post_indel_call_ch

  shell:
  '''
  lofreq indelqual --dindel -f !{params.refseq_fasta} -o !{input_bam}.indelqual.bam !{input_bam}

  lofreq index !{input_bam}.indelqual.bam
                                                                                
  lofreq call-parallel --pp-threads !{task.cpus} --no-default-filter --call-indels --only-indels -f !{params.refseq_fasta} !{input_bam}.indelqual.bam > !{input_bam}.indel.pre_vcf
                                                                                
  lofreq filter -v !{params.min_depth_for_variant_call} -V 0 -a !{params.min_allele_freq} -A 0 --no-defaults -i  !{input_bam}.indel.pre_vcf -o !{input_bam}.indel.vcf
  '''
}



/* 
 Call consensus sequences for each dataset 

 Here we are using mpileup for variant calling instead of lofreq 
 because I couldn't get lofreq vcf to work as input for bcftools

*/
process call_dataset_consensus {
  publishDir "${params.consensus_out_dir}", mode:'link'

  input:
  tuple val(sample_id), path(input_bam) from post_bsqr_consensus_ch

  output:
  tuple val(sample_id), path("${sample_id}_consensus.fasta") into post_merge_vcf_ch

  script:
  """
  # consensus calling according to the strategy outlined here: https://www.biostars.org/p/367626/
  samtools sort $input_bam -o sorted_input.bam
  samtools mpileup -uf ${params.refseq_fasta} sorted_input.bam | bcftools call -c | vcfutils.pl vcf2fq > consensus.fastq
  # Convert .fastq to .fasta and set bases of quality lower than 20 to N
  # the sed here adds in sample ID to the fasta header line, otherwise they all just have the same name 
  seqtk seq -aQ64 -q20 -n N consensus.fastq | sed s/">"/">${sample_id}_"/ > ${sample_id}_consensus.fasta
  """
}



/* 
  Consensus calling using same lofreq vcfs that defined called variants

  ** couldn't get lofreq vcf to work as input for bcftools ** 

*/
/*
process call_dataset_consensus {
  publishDir "${params.outdir}", mode:'link'

  input:
  tuple val(sample_id), path(vcfs) from post_indel_call_consensus_ch
    .concat(post_snv_call_consensus_ch)
    .groupTuple()

  output:
  tuple val(sample_id), path("${sample_id}_consensus.fasta") into post_merge_vcf_ch

  script:
  """
  # convert to compressed/indexed vcf to make bcftools happy
  bcftools view ${vcfs[0]} -O z > ${vcfs[0]}.gz
  bcftools index ${vcfs[0]}.gz
  bcftools view ${vcfs[1]} -O z > ${vcfs[1]}.gz
  bcftools index ${vcfs[1]}.gz

  bcftools concat ${vcfs[0]}.gz ${vcfs[1]}.gz | bcftools call -c | vcfutils.pl vcf2fq > consensus.fastq
  seqtk seq -aQ64 -q20 -n N consensus.fastq > ${sample_id}_consensus.fasta
  """
}
*/


/* 
 use snpEff to annotate vcfs 
 */
process annotate_variants {
  publishDir "${params.vcf_out_dir}", mode:'link'

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

  input:
  path(all_depth) from analyze_variants_depth_ch
  path(snp_sifts) from post_prepend_snp_sift_ch.collect()

  output:
  path("variant_summary.xlsx") 
  path("sample_correlation_heatmap.pdf") optional true

  script:
  """
  Rscript ${params.R_bindir}/analyze_snpeff_variants.R ${params.R_bindir} \
    ${params.min_allele_freq} \
    $all_depth \
    ${params.min_depth_for_variant_call} \
    $snp_sifts
  """
}


process tabulate_fastq_counts {
  publishDir "${params.outdir}", mode: 'link'

  input:
  path(all_count_files) from post_count_initial_ch.concat(post_count_trim_ch, post_count_host_ch, post_count_refseq_aligned_ch).collect()

  output:
  path ("all_read_counts.txt") 
  path ("filtering_plots.pdf") 

  script:
  """
  Rscript ${params.R_bindir}/process_fastq_counts.R ${params.R_bindir} ${all_count_files}
  """
}
