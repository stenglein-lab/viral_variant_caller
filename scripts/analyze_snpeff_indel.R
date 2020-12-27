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
} else {
  # for troubleshooting via RStudio interface
  r_bindir <- "."
  snp_sifts = list.files(path = "results", pattern = "*.snp_sift", full.names = T)
}

# takes one tsv file created by snpsift
# and returns a dataframe with information about the variants
variant_df <- function(snp_sift_file_name) {
  
  # read the file in
  df <- read.delim(snp_sift_file_name, header = T)
  
  as.data.frame(df)
  names(df) <- c("Chromosome", "Position", "Reference", "Variant", "Frequency", "Depth", "Strand_bias", "Indel", "Effect", "Impact", "Gene")
  
  # add the sample name as a column
  df <- df %>% mutate(sample_name = snp_sift_file_name)
  
  return(df)
  
}

# this runs variant_df on all the file names and returns a list of dataframes
# see: https://www.rdocumentation.org/packages/base/versions/3.3.1/topics/lapply?
df <- lapply(snp_sifts, "variant_df")

# this rbinds all the dfs together to make one big tidy df
# see: https://stat.ethz.ch/pipermail/r-help/2007-April/130594.html
df <- do.call("rbind", df)

# extract dataset ID from file name 
df <- df %>% mutate(Dataset_ID = str_replace(sample_name, "results/", "") )
df <- df %>% mutate(Dataset_ID = str_replace(Dataset_ID, ".wa1.bam.indel.vcf.snp_eff.snp_sift", ""))

#?as.numeric
# make sure numeric info is numeric
df$Position <- as.numeric(as.character(df$Position))
df$Allele_frequency <- as.numeric(as.character(df$Frequency))
df$Depth <- as.numeric(as.character(df$Depth))

# output a summary table

df_wide <- df %>% pivot_wider(id_cols = c(Chromosome, Position, Reference, Variant), 
                              names_from=Dataset_ID, 
                              #values_from=c(allele_freq, depth),
                              values_from=c(Frequency, Depth),
                              names_sort = T)

#?pivot_wider
df_wide <- df_wide %>% arrange(Position)

wb <- createWorkbook("Structural_variant_snpeff_summary.xlsx")
addWorksheet(wb, "structural_variants")
writeData(wb, "structural_variants", df_wide,borders="all")
?writeData
saveWorkbook(wb, "Structural_variant_snpeff_summary.xlsx", overwrite = TRUE)
