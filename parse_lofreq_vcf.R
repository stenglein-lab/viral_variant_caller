library(tidyverse)
library(vcfR)

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

