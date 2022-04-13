# viral_variant_caller / cdphe_sequencing branch

This repository branch contains a pipeline to generate SARS-CoV-2 sequences for a project to do SARS-CoV-2 genome sequencing from Colorado.

### To analyze a dataset

To download this pipeline, run:

```
git clone -b cdphe_sequencing https://github.com/stenglein-lab/viral_variant_caller.git
```


### Usage

This pipeline has 2 ways to handle dependencies: (1) using singularity images, or (2) using conda environments.  Singularity is the preferred method.

#### Singularity containers
  
The pipeline can use singularity containers to run programs like bowtie2 and bwa.  To use these containers you must be running the pipeline on a computer that has [singularity](https://sylabs.io/singularity) [installed](https://sylabs.io/guides/latest/admin-guide/installation.html).  To test if singularity is installed, run:

```
singularity --version
```

To run with singularity containers include the option `-profile singularity` in the nextflow command line, for instance:

```
nextflow run main.nf -resume -profile local,singularity --fastq_dir ../2022_3_1_run_4_fastq
```
Singularity containers will be automatically downloaded and stored in a directory named `singularity_cacheDir` in your home directory.  They will only be downloaded once.

This invocation:
- takes advantage of nextflow's built-in ability to resume pipelines (-resume flag)
- specifies the local and singularity profiles, both specified in [nextflow.config](./nextflow.config)
- points to the directory containing the fastq files that will be analyzed (--fastq_dir parameter)

#### Conda environment

The pipeline can also use an all-in-one conda environment.  This requires conda to be installed on your computer.  To test if conda is installed, run:

```
conda -V
```

The conda environments needed are defined in [yaml files in this directory](./environment_setup/) and will be automatically created if you run the pipeline using the conda profile.  To run the pipeline with conda, include `-profile conda` in the nextflow command line, for instance:

```
nextflow run main.nf -resume -profile local,conda --fastq_dir ../2022_3_1_run_4_fastq
```

The conda environment will be created in a directory in your home directory named `conda_cacheDir` and will only be created once.

You should specify either `-profile conda` or `-profile singularity` or the pipeline will output an error message and halt.  

#### Running with screen

[See here](https://github.com/stenglein-lab/taxonomy_pipeline/blob/master/docs/tutorial.md#section_screen) for an explanation of using the screen utility to avoid dropped connections and premature termination of a pipeline.


### Viral reference sequence

This pipeline is by default set up to call variants relative to a SARS-CoV-2 reference genome, but it could be reconfigured pretty easily to call variants relative to a different reference sequence.  This can be done by changing these parameters in [nextflow.config](./nextflow.config) (or by [overriding them](https://www.nextflow.io/docs/latest/cli.html#pipeline-parameters) on the command line).  

```
  // SARS-CoV-2 Wuhan-1 (NC_045512) reference sequence,
  // or a reference seq of your choosing
  refseq_name          = "NC_045512"
  refseq_fasta         = "${refseq_dir}/${refseq_name}.fasta"
  refseq_genbank       = "${refseq_dir}/${refseq_name}.gb"
```

The pipeline is expecting your reference sequence to exist in the refseq directory in fasta and genbank format.  These can both be downloaded from NCBI, or exported from software like geneious.  The refseq directory should be populated with appropriate files.  So, in the above example, the pipeline is expecting the refseq directory to contain files named NC_045512.fasta and NC_045512.gb.  This is the RefSeq SARS-CoV-2 sequence.  

### Host filtering

A step in this pipline removes host-derived reads.  By default it maps against the human genome.  To modify which host genome is used, modify these parameters in `nextflow.config` or by overriding these values as [command line arguments to nextflow](https://www.nextflow.io/docs/latest/cli.html#pipeline-parameters).

```
  // Human samples: use human genome for host filtering
  host_bt_index      = "/home/databases/human/GCRh38"
```

You will need an [existing bowtie2 index](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer) to exist for the pipeline to do host filtering.  In the above example, a bowtie2 index named GCRh38 exists in the directory `/home/databases/human/`.

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
- [Pangolin](https://cov-lineages.org/resources/pangolin.html)

These dependencies are handled via the conda environment [defined](./environment_setup/variant_conda_environment.yaml) in this repository.  To use this conda environment, you must run the pipeline with `-profile conda` specified.   The first time you run with this profile, the pipeline will create the environment and cache it in a directory named conda_cacheDir in your home directory.  You can override this cache location by changing this parameter defined in `nextflow.config`:

```
    conda.cacheDir         = "$HOME/conda_cacheDir"
```

### Results

The results of this pipeline will be output to a results subdirectory.  Results include:

- **consensus_sequences.fasta:** consenensus sequences ready to be deposited to GISAID or another database
- **gisaid_submission.xlsx:** an excel metadata file for GISAID submission.  Note this has to be manually converted to an .xls format file by opening and saving from Excel (GISAID can't handle .xlsx files and the openxlsx R package can't write .xls files).
- **dataset_summary.xlsx:** an excel file containing information about all the samples, even those that are not complete enough to submit to GISAID.
- **initial_qc_report.html, post_trim_qc_report.html:** QC reports of raw data before and after adapter/quality trimming
- **coverage_plots.pdf:** coverage plots over virus reference sequence 
- **Median_dephts.xlsx:** median coverage values over virus reference sequence in each dataset
- **filtering_plots.pdf:** plots of #/fraction of reads remaining after various filtering steps

