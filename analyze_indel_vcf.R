library(tidyverse)
library(vcfR)
library(ape)
library(openxlsx)
library(readxl)

# TODO: 
# - figure out how to source a file from nextflow
# - parse out metadata from sample IDs and join it (transpose) to the output data frame
# - ignore variants < 3%

# if being run from the command line
if (!interactive()) {

  args = commandArgs(trailingOnly=TRUE)
  
  if (length(args)==0) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
  } 
  
  # vcf file names passed as command line arguments
  r_bindir=args[1]
  # -c(1) --> all but the first element of the list
  vcfs = args[-c(1)]
} else {
  # for troubleshooting via RStudio interface
  r_bindir <- "."
  vcfs = list.files(path = "results", pattern = "*.vcf", full.names = T)
}


# get the parse_vcf function
source(paste0(r_bindir, "/parse_lofreq_vcf.R"))

# this runs parse_vcf on all the vcf names and returns a list of dataframes
# see: https://www.rdocumentation.org/packages/base/versions/3.3.1/topics/lapply?
df <- lapply(vcfs, "parse_vcf")

# this rbinds all the dfs together to make one big tidy df
# see: https://stat.ethz.ch/pipermail/r-help/2007-April/130594.html
df <- do.call("rbind", df)

# extract dataset ID from vcf name 
df <- df %>% mutate(dataset_id = str_replace(vcf_name, "results/", "") )
df <- df %>% mutate(dataset_id = str_replace(dataset_id, ".wa1.bam.indel.vcf", ""))

# read in metadata 
# metadata <- read_excel("Metadata.xlsx")
# df <- left_join(df, metadata, by="dataset_id")
# ?left_join


# are these indels insertions or a deletion?
# if the reference bases are longer than the alt bases, this is a deletion
# since these are indels the other option is that it's an insertion
df <- df %>% mutate(vcf_type = 
                      if_else(str_length(as.character(REF)) > str_length(as.character(ALT)), 
                              "deletion", 
                              "insertion"))

# omit all indels with a frequency < 3%
min_allele_freq = 0.03
df <- df %>% filter(allele_freq >= min_allele_freq)

?as.numeric
# make sure numeric info is numeric
df$POS <- as.numeric(as.character(df$POS))
df$allele_freq <- as.numeric(as.character(df$allele_freq))
df$depth <- as.numeric(as.character(df$depth))

# get DNA sequence and also gff file
# dna <- ape::read.dna("viral_refseq/wa1.fasta", format = "fasta")
# gff <- read.table("viral_refseq/wa1.gff.nofasta", sep="\t", quote="")

# output a summary table

# df_biomarch <- df %>% filter()


df_wide <- df %>% pivot_wider(id_cols = c(CHROM, POS, REF, ALT), 
                              names_from=dataset_id, 
                              #values_from=c(allele_freq, depth),
                              values_from=c(allele_freq, depth),
                              names_sort = T)

# remove allele_freq_ from column names
colnames(df_wide) <- str_replace(colnames(df_wide), "allele_freq_", "")
colnames(df_wide) <- str_replace(colnames(df_wide), "CHROM", "Reference_sequence")
colnames(df_wide) <- str_replace(colnames(df_wide), "POS", "Position")
colnames(df_wide) <- str_replace(colnames(df_wide), "REF", "Reference_base(s)")
colnames(df_wide) <- str_replace(colnames(df_wide), "ALT", "Variant_base(s)")


?pivot_wider
df_wide <- df_wide %>% arrange(Position)

wb <- createWorkbook("Structural_variant_summary.xlsx")
addWorksheet(wb, "structural_variants")
writeData(wb, "structural_variants", df_wide,borders="all")
?writeData
saveWorkbook(wb, "Structural_variant_summary.xlsx", overwrite = TRUE)



# how do replicates look?
# df <- df %>% mutate(replicate = str_match(dataset_id,"(*)_Rep([1-9])")[,2])
      # df %>% mutate(replicate = str_match(dataset_id,"(\\*)_Rep([1-9])")[,2])
      # df %>% mutate(dataset)

# df %>% group_by(CHROM, POS, REF, ALT, )

ggplot(df) + geom