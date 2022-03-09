library(tidyverse)
library(readxl)
library(openxlsx)

# This script prepares output files that can be submitted to GISAID for sequence deposition
# 
# Mark Stenglein 3/7/2022

# deal with expected input
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  # TODO: check CLAs
  consensus_fasta = args[1]
  dataset_summary = args[2]
  obfuscated_ids = args[3]
  min_fraction_called = args[4]
  metadata_file = args[5]
  output_dir="./"
} else {
  # if running via RStudio (for developing/troubleshooting)
  # put in some sensible defaults
  consensus_fasta = "../results/consensus_sequences_original_ids.fasta"
  dataset_summary = "../results/dataset_summary.xlsx"
  obfuscated_ids = "../results/obfuscated_sample_ids.txt"
  min_fraction_called = 0.95
  metadata_file = "../input/gisaid_metadata.tsv"
  output_dir="../results/"
}

# read in input files
dataset_summary <- read_excel(dataset_summary)
obfuscated_ids <- read.delim(obfuscated_ids, sep="\t", header=T)
metadata_in <- read.delim(metadata_file, sep="\t", header=F, comment.char = "#", strip.white = T , stringsAsFactors = F)

# handle metadata file:
# a 2-column tab-delimited file that maps gisaid fields to their values
# convert to a named list for easy access to values
# see: https://stackoverflow.com/questions/33418288/how-to-convert-a-matrix-to-dictionary-like-a-list
colnames(metadata_in) <- c("key", "value")
metadata <- as.list(t(metadata_in[, "value"]))
names(metadata) <- t(metadata_in[, "key"])

# filter out rows below cutoff and any NTC data
dataset_pass <- dataset_summary %>% 
  filter(fraction_called >= min_fraction_called & !str_detect(dataset, "NTC"))

# join in obfuscated IDS
dataset_pass <- left_join(dataset_pass, obfuscated_ids, by=c("dataset" = "original_id"))

# prepare data for upload to gisaid
# these names taken from GISAID excel template 
gisaid_columns <- c("submitter", "fn","covv_virus_name","covv_type",
                    "covv_passage", "covv_collection_date","covv_location","covv_add_location",
                    "covv_host","covv_add_host_info","covv_sampling_strategy","covv_gender",
                    "covv_patient_age","covv_patient_status","covv_specimen","covv_outbreak",
                    "covv_last_vaccinated","covv_treatment","covv_seq_technology","covv_assembly_method",
                    "covv_coverage","covv_orig_lab","covv_orig_lab_addr","covv_provider_sample_id",
                    "covv_subm_lab","covv_subm_lab_addr","covv_subm_sample_id","covv_authors")

# more descriptive versions of column names, from GISAID
gisaid_long_names <- c("Submitter", "FASTA filename","Virus name","Type",
                       "Passage details/history","Collection date","Location","Additional location information",
                       "Host","Additional host information","Sampling Strategy","Gender",
                       "Patient age","Patient status","Specimen source","Outbreak",
                       "Last vaccinated","Treatment","Sequencing technology","Assembly method",
                       "Coverage","Originating lab","Address","Sample ID given by originating laboratory",
                       "Submitting lab","Address","Sample ID given by the submitting laboratory","Authors")

# how many rows do we need?
# the # of sequences above min. completeness cutoff plus one for the long names row
num_rows_needed <- nrow(dataset_pass) 

# make a data frame with the appropriate columns that GISAID wants
gisaid_df <- data.frame(matrix(ncol = length(gisaid_columns), nrow = num_rows_needed))
colnames(gisaid_df) <- gisaid_columns

# prepare a 2nd row containing more descriptive names 
# gisaid wants this second row in the submitted metadata file
long_names_row <- data.frame(matrix(ncol = length(gisaid_columns), nrow = 0))
colnames(long_names_row) <- gisaid_columns
long_names_row[nrow(long_names_row) + 1,] <- gisaid_long_names


# this function pulls out date and returns it in a GISAID-compatible format
prepare_gisaid_date <- function(dataset_id){
  
  # pull out date in GISAID format
  # parse out date metadata from sample IDs
  # expected to be in the sample ID (dataset name)
  # in the format: YYYYMMDD or YYYY
  date_regex <- "_(20[0-9]{2})([0-9]{2})([0-9]{2})"
  date_fields <- str_match(dataset_id, date_regex)
  
  whole_match = date_fields[,1]
  year        = date_fields[,2]
  month       = date_fields[,3]
  day         = date_fields[,4]
  if (is.na(year)) {
    
    message(paste0("WARNING: could not parse year out of sample ID: ", dataset_id))
    message(       "         this will result in invalid GISAID submission data.")
    gisaid_date <- NA_character_
    
  } else if (!is.na(year) & !is.na(month) & !is.na(day)) {
    
    # the field covv_collection_date should be in the format:
    # 2021-03-31
    gisaid_date <- paste0(year, "-", month, "-", day)
    
  }
  else {
    
    # or if no month/day info, just the year
    # 2021
    gisaid_date <- year
    
  }
  
  gisaid_date
}

# parse date from dataset ID and return in gisaid-expected format
dataset_pass$covv_collection_date <- lapply(dataset_pass$dataset, prepare_gisaid_date)

# pull year out of gisaid date
dataset_pass$year <- str_match(dataset_pass$covv_collection_date, "(^[:digit:]{4})")[,2]

# per GISAID conventions, name should be in the format:
# hCoV-19/Country/Identifier/Year
# e.g. 
# hCoV-19/USA/Colorado_CSU_XYZ/2021
dataset_pass <- dataset_pass %>% 
  mutate(covv_virus_name = paste0("hCoV-19/USA/Colorado_CSU_", obfuscated_id, "/", year))

# covv_coverage 
dataset_pass <- dataset_pass %>% 
  mutate(covv_coverage <- mean_depth)


# ---------------
# Consensus FASTA 
# ---------------
# read in the existing fasta file into one string
# this fasta doesn't have GISAID-formatted sequence names
fasta_string <- readChar(consensus_fasta, file.info(consensus_fasta)$size)

  
# for each dataset ID
# replace the name in the fasta file with the GISAID-formatted ID

replaced_fasta_string <- fasta_string

for (row in 1:nrow(dataset_pass)) {
  original_id <- filter(dataset_pass, row_number() == row) %>% pull(dataset)
  replacement_id <- filter(dataset_pass, row_number() == row) %>% pull(covv_virus_name)
  
  # replace one fasta header
  # this assumes that sequences have unique IDs
  replaced_fasta_string <- str_replace(replaced_fasta_string, original_id, replacement_id)
}

# output fasta with replaced names
cat(replaced_fasta_string, file=paste0(output_dir, metadata$fn), sep="")

# fill in ouput data table

# fill in fixed data from the metadata input
for (name in names(metadata)){
  gisaid_df[name] <- metadata[name]
}

# fill in dynamic data from dataset info
for (name in names(dataset_pass)){
  if (name %in% names(gisaid_df)) {
    gisaid_df[name] <- dataset_pass[name]
  }
}

# add in the 2nd row of descriptive names
gisaid_df <- rbind(long_names_row, gisaid_df)

# write out the excel file 
wb <- createWorkbook(paste0(output_dir, "gisaid_submission.xls"))
modifyBaseFont(wb, fontSize = 11, fontColour = "black", fontName = "Helvetica")
addWorksheet(wb, "Submissions")
writeData(wb, "Submissions", gisaid_df, na.string=NA_character_)
saveWorkbook(wb, paste0(output_dir, "gisaid_submission.xlsx"), overwrite = TRUE)




