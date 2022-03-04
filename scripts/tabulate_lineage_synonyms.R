library(tidyverse)
library(jsonlite)

# This script will output a tsv file that maps
# pango lineages to WHO synonyms.  
#
# E.g.  B.1.1.529 (Pango) == Omicron (WHO)
# 
# The input to this script is a directory of json files from the 
# cov-lineages constellations repository:
#
# See: https://github.com/cov-lineages/constellations 
#
# These json files contain the mapping information, which this script 
# extracts and tabulates
#
# This script outputs to standard output
# 
# Mark Stenglein, 3/1/2022
#

# deal with expected input
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  json_dir = args[1]
  output_dir="./"
} else {
  # if running via RStudio
  json_dir <- "./constellations/constellations/definitions/" 
  output_dir="../results/"
}

# list all the json files in the json-containing directory 
files <- list.files(json_dir, pattern="*.json")

# if there are any json files, output header
if (length(files) > 0) {
  header_text <- sprintf("pango_lineage\twho_label\n")
  cat(header_text)
}

# for each json
for (file in files) {
  
  filename <- paste0(json_dir, file)
  
  json_df <- read_json(filename)

  # if there is a WHO label for this constellation
  if ("WHO_label" %in% names(json_df$variant)) {
    
    who_label = json_df$variant$WHO_label
    
    pango_lineages = json_df$variant$Pango_lineages

    # a constellation may contain multiple pango lineages
    # output a line for each
    for (lin in pango_lineages) {
      output_text <- sprintf("%s\t%s\n", lin, who_label)
      cat(output_text)
    }
  }
}


