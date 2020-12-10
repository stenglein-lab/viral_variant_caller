library(tidyverse)
library(vcfR)
library(ape)
library(openxlsx)

# TODO: 
# - parse out metadata from sample IDs and join it (transpose) to the output data frame

# source("parse_lofreq_vcf.R")

# parses one vcf file created by lofreq
# and returns a dataframe with information about the variants
parse_vcf <- function(vcf_file_name) {

  # read the VCF file using the read.vcfR function from vcfR package
  vcf <- read.vcfR( vcf_file_name,
                    verbose = T,
                    cols = 8)

  # this is metadata from vcf file
  # don't handle for now
  # vcf@meta
  # queryMETA(vcf)
  # queryMETA(vcf, element = 'DP4')
  
  # fixed data
  # the columns with refseq, position, etc
  fix_df <- getFIX(vcf)
  fix_df
  
  # info: the DP, AF, etc fields in the 8th column
  info <- getINFO(vcf)
  
  # merge columns
  vcf_df <- as.data.frame(cbind(fix_df, info))
  
  # parse out total depth and allele frequency from info column
  # TODO: case for more than one allele? 
  dp_af <- str_match(vcf_df$info, "DP=(\\d+);AF=([+-]?([0-9]*[.])?[0-9]+)" )
  vcf_df$depth=dp_af[,2]
  vcf_df$allele_freq=dp_af[,3]
  
  # parse out DP4 info: coverage for reference and alternative alleles on forward and reverse mapping reads
  dp4 <- str_match(vcf_df$info, "DP4=(\\d+),(\\d+),(\\d+),(\\d+);")
  vcf_df$dp4_ref_f = dp4[,2]
  vcf_df$dp4_ref_r = dp4[,3]
  vcf_df$dp4_alt_f = dp4[,4]
  vcf_df$dp4_alt_r = dp4[,5]
  
  
  # add the sample name as a column
  vcf_df <- vcf_df %>% mutate(vcf_name = vcf_file_name)
  
  return(vcf_df)

}



# vcfs = list.files(path = "./results", pattern = "*.indel.vcf$", full.names=T)

args = commandArgs(trailingOnly=TRUE)

# TODO: get from a command line argument?
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

# vcf file names passed as command line arguments
vcfs = args

# vcfs <- rep("N19_1_20_INF_pass4_dpi3_flask1_tubeid6_R1_R1_fh.fastq.indel.vcf", 2)

# this runs parse_vcf on all the vcf names and returns a list of dataframes
# see: https://www.rdocumentation.org/packages/base/versions/3.3.1/topics/lapply?
df <- lapply(vcfs, "parse_vcf")

# this rbinds all the dfs together to make one big tidy df
# see: https://stat.ethz.ch/pipermail/r-help/2007-April/130594.html
df <- do.call("rbind", df)

# are these indels insertions or a deletion?
# if the reference bases are longer than the alt bases, this is a deletion
# since these are indels the other option is that it's an insertion
df <- df %>% mutate(vcf_type = 
                      if_else(str_length(as.character(REF)) > str_length(as.character(ALT)), 
                              "deletion", 
                              "insertion"))

# get DNA sequence and also gff file
# dna <- ape::read.dna("viral_refseq/wa1.fasta", format = "fasta")
# gff <- read.table("viral_refseq/wa1.gff.nofasta", sep="\t", quote="")

# output a summary table

df_wide <- df %>% pivot_wider(id_cols = c(CHROM, POS, REF, ALT), names_from=vcf_name, values_from=c(allele_freq, depth))

wb <- createWorkbook("Structural_variant_summary.xlsx")
addWorksheet(wb, "structural_variants")
writeData(wb, "structural_variants", df_wide)
saveWorkbook(wb, "Structural_variant_summary.xlsx", overwrite = TRUE)


