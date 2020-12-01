# viral_variant_caller

This repository contains a pipeline for calling single nucleotide and structural variants in viral populations. 


### Installation

To download this pipeline, run:

```
git clone https://github.com/stenglein-lab/viral_variant_caller.git
```

### Usage

To run this pipeline, create and populate a directory named `fastq` in the `viral_variant_caller` directory with the datasets (fastq) files that you wish to analyze.  Then run:

```
nextflow run variant_pipeline.nf 
```


#### Viral reference sequence

This pipeline is by default set up to call variants relative to the SARS-CoV-2 WA1 genome.  This can be changed by changing these parameters in `variant_pipeline.nf`:

```
// SARS-CoV-2 wa1 refseq                                                        
params.viral_refseq_name = "wa1"                                                
params.viral_fasta = "${baseDir}/viral_refseq/${params.viral_refseq_name}.fasta"
params.viral_gb = "${baseDir}/viral_refseq/${params.viral_refseq_name}.gb"      
params.viral_gff = "${baseDir}/viral_refseq/${params.viral_refseq_name}.gff"  
```

(The viral_refseq directory should be populated with appropriate sequence files in .fasta, .gb, and .gff format.)

#### Host filtering

A step in this pipline removes host-derived reads.  To modify the host genome used, modify these parameters in the main nextflow pipeline file:
```
params.host_bt_index = "/home/databases/primates/agm_genome"                    
params.host_bt_suffix = "agm_genome"                                            
```


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

We are working to integrate a conda environment to facilitate use of these tools for this pipeline but for the moment, these tools must be installed and in the user's PATH.
