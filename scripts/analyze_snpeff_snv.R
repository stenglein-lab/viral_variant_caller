library(tidyverse)
library(ape)
library(openxlsx)
library(readxl)

# if being run from the command line
if (!interactive()) {
  
  args = commandArgs(trailingOnly=TRUE)
  
  if (length(args)==0) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
  } 
  
  # snp sift file names passed as command line arguments
  r_bindir=args[1]
  # -c(1) --> all but the first element of the list
  snp_sifts = args[-c(1)]
}

# takes one tsv file created by snpsift
# and returns a dataframe with information about the variants
variant_df <- function(snp_sift_file_name) {
  
  # read the file in
  df <- read.delim(snp_sift_file_name, header = T)
  
  as.data.frame(df)
  names(df) <- c("Chromosome", "Position", "Reference", "Variant", "Frequency", "Depth", "Strand_bias", "Indel", "Effect", "Impact", "Gene")
  
  # add the sample name as a column
  df <- df %>% mutate(sample_ID = snp_sift_file_name)
  
  return(df)
  
}

# this runs variant_df on all the file names and returns a list of dataframes
# see: https://www.rdocumentation.org/packages/base/versions/3.3.1/topics/lapply?
df <- lapply(snp_sifts, "variant_df")

# this rbinds all the dfs together to make one big tidy df
# see: https://stat.ethz.ch/pipermail/r-help/2007-April/130594.html
df <- do.call("rbind", df)

# extract dataset ID from file name 
df <- df %>% mutate(sample_ID = str_replace(sample_ID, "results/", "") )
df <- df %>% mutate(sample_ID = str_replace(sample_ID, ".wa1.bam.snv.vcf.snp_eff.snp_sift", ""))
df <- df %>% select(sample_ID, everything())
names(df)[1] <- "Dataset_ID"

#?as.numeric
# make sure numeric info is numeric
df$Position <- as.numeric(as.character(df$Position))
df$Frequency <- as.numeric(as.character(df$Frequency))
df$Depth <- as.numeric(as.character(df$Depth))

# omit all variants with a frequency < 3%
min_allele_freq = 0.03
df <- df %>% filter(Frequency >= min_allele_freq)

# output a summary table

wb <- createWorkbook("Single_nucleotide_variant_snpeff_summary.xlsx")
addWorksheet(wb, "single_nucleotide_variants")
writeData(wb, "single_nucleotide_variants", df,borders="all")
?writeData
saveWorkbook(wb, "Single_nucleotide_variant_snpeff_summary.xlsx", overwrite = TRUE)