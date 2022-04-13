library(tidyverse)
library(DescTools)

# this script optionally obfuscates sample IDs 
# using the Vigenere cypher
# 
# Although the sample IDs are not identifying, they do contain info
# like study subject # 
#
# So obfucscating them adds an additional layer of separation 
# in public database IDs
#
# If you wish to not obfuscate IDs, just leave key_file as NA 
#
# Mark Stenglein March 6, 2022


# deal with expected input
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  # pop the first argument: this will be the path to a file containing a key
  # this key file will not be distributed with the repository
  key_file <- args[1]
  args <- args[-(1)]
  sample_ids = args
  output_dir = "./"
} else {
  # if running via RStudio
  key_file <- NA
  sample_ids =  c("AA_B123_20210132", "BB_C234_20220111", "CC_D456_20210304", "SomeOtherFormat")
  output_dir =  "../results/"
}

# -----------------------
# optional ID obfuscation
# -----------------------
# obfuscation key
key <- NA
if (file.exists(key_file)) {
  key <- readChar(key_file, file.info(key_file)$size)
  key <- str_replace_all(key, "[\r\n]" , "")
} else {
  message("INFO: sample ID obfuscation key file not found.  Will not obfuscate sample IDs.")
}

# generate a sample ID for one sample
generate_sample_id <- function (sample_id, obfuscation_key){
  
  # expect samples to be in the format
  # XX_XXXXX_20210531
  # (date is last part)
  
  id_regex <- "([:alpha:]+)_([:alnum:]+)_([:digit:]+)"
  matches <- str_match(sample_id, id_regex)
  
  # if we were able to parse the expected format
  if (!is.na(matches[,2]) & !is.na(matches[,3]) & !is.na(matches[,4])) {
    id_1 <- matches[,2]
    id_2 <- matches[,3]
    date <- matches[,4]
    
    # use the Vigenere function to "encrypt" the ID parts
    # since this is not strong encryption, call it obfuscation
    # if no key, will leave IDs alone
    if (!is.na(obfuscation_key)) {
      id_1 = Vigenere(id_1, key=obfuscation_key)
      id_2 = Vigenere(id_2, key=obfuscation_key)
    }
    
    new_id <- paste0(id_1, id_2, "_", date)
  }
  else {
    # ID wasn't in expected format: leave as is
    new_id <- sample_id
  }
} 

# obfuscate all IDs
obfuscated_ids <- lapply(sample_ids, generate_sample_id, key)

# prepare output
output_df <- cbind(sample_ids, obfuscated_ids)
colnames(output_df) <- c("original_id", "obfuscated_id")

# output tab-delimited table of old->new IDs
write.table(output_df, "obfuscated_sample_ids.txt", quote=F, row.names=F, col.names=T, sep="\t")

