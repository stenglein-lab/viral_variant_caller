# viral_variant_caller

This repository contains a pipeline for calling single nucleotide and structural variants in viral populations. 


### Installation

To download this pipeline, run:

```
git clone https://github.com/stenglein-lab/viral_variant_caller.git
```

### Usage

Before you run this pipeline, you'll have to make sure that dependencies (software on which this pipeline depends) are in place.  The best way to do this is to use the included [conda environment](./environment_setup/variant_conda_environment.yaml).  Please see the [instructions to install and setup this environment](./environment_setup/README.md) for instructions on how to install and initialize this environment.

To run the pipeline, you'll need to move your datasets (fastq files) into a new directory named `fastq` in the `viral_variant_caller` directory with the datasets (fastq) files that you wish to analyze.  Then run:

```
nextflow run variant_pipeline.nf 
```

To resume the pipeline if it stopped for some reason

```
nextflow run variant_pipeline.nf -resume
```

[See here](https://github.com/stenglein-lab/taxonomy_pipeline/blob/master/docs/tutorial.md#section_screen) for an explanation of using the screen utility to avoid dropped connections and premature termination of a pipeline.


### Viral reference sequence

This pipeline is by default set up to call variants relative to a SARS-CoV-2 reference genome, but it could be reconfigured pretty easily to call variants relative to a different reference sequence.  This can be done by changing these parameters in `variant_pipeline.nf`:

```
// SARS-CoV-2 USA/WA1 (Genbank accession MN985325), Wuhan-1 (NC_045512) reference sequence, 
// or a reference seq of your choosing                                          
params.refseq_dir = "${baseDir}/refseq/"                                        
params.refseq_name = "MN985325"                                                 
// params.refseq_name = "NC_045512"                                             
params.refseq_fasta = "${params.refseq_dir}/${params.refseq_name}.fasta"        
params.refseq_genbank = "${params.refseq_dir}/${params.refseq_name}.gb" 
```

The pipeline is expecting your reference sequence to exist in the refseq directory in fasta and genbank format.  These can both be downloaded from NCBI, or exported from software like geneious.  The refseq directory should be populated with appropriate files.  So, in the above example, the pipeline is expecting the refseq directory to contain files named MN985325.fasta and MN985325.gb.  This is the USA/WA1 SARS-CoV-2 sequence.  

### Host filtering

A step in this pipline removes host-derived reads.  To modify the host genome used, modify these parameters in the main nextflow pipeline file:
```
params.host_bt_index = "/home/databases/primates/agm_genome"                    
params.host_bt_suffix = "agm_genome"                                            
```

(agm here stands for African green monkey, because when this pipeline was developed, we were analyzing virus grown in Vero cells, which are from that species of monkey.

You will need an [existing bowtie2 index](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer) to exist for the pipeline to do host filtering.  


### Dependencies

This pipeline depends on the scripts in this repository as well as the following tools:

- [nextflow](https://www.nextflow.io/)
- [lofreq](https://csb5.github.io/lofreq/)
- [bwa](https://github.com/lh3/bwa)
- [gatk](https://gatk.broadinstitute.org/hc/en-us)
- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [multiqc](https://multiqc.info/)
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [samtools](http://samtools.github.io/)
- [R](https://www.r-project.org/)

These dependencies are handled via the conda environment defined in this repository.  [See here for instructions](./environment_setup/README.md) on creating this environment. 


### Results

The results of this pipeline will be output to a results subdirectory. The main pipeline outputs are:

- **variant_summary.xlsx:** this spreadsheet contains a matrix of variant frequencies of all detected variants in all datasets. 
- **consensus_sequences:** a directory containing a consensus sequence for each dataset (variants >50% will change consensus relative to the original refseq)
- **sample_correlation_heatmap.pdf:** A heatmap of variants and clustering of variants and datasets
- **initial_qc_report.html, post_trim_qc_report.html:** QC reports of raw data before and after adapter/quality trimming
- **coverage_plots.pdf:** coverage plots over virus reference sequence 
- **Median_dephts.xlsx:** median coverage values over virus reference sequence in each dataset
- **filtering_plots.pdf:** plots of #/fraction of reads remaining after various filtering steps

[See here](https://github.com/stenglein-lab/taxonomy_pipeline/blob/master/docs/tutorial.md#-transferring-files-to-and-from-servers) for an explanation of how to transfer files from a server using a utility like sftp or cyberduck

#### A note on the snpEff database

I was getting 'ERROR_CHROMOSOME_NOT_FOUND' errors in the snpEff vcf output, and to fix it had to modify the genbank file I downloaded from genbank to change the version field from: 
```
VERSION     MN985325.1
```

to:

```
VERSION     MN985325
```

This could possibly be avoided by using GTF annotation instead of GB?
