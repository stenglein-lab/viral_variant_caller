library(DescTools)

# this script encrypts sample IDs using the Vigenere algorithm
# as implemented in the DescTools package
# 
# Although the sample IDs are not identifying, they do contain info
# like study subject # 
# this just adds an additional layer of separation between the people
# and the data in public databases
#
# Mark Stenglein Jan 26, 2022


# deal with expected input
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  # pop the first argument: this will be the path to a file containing an encryption key
  # this key file will not be distributed with the repository
  key_file <- args[1]
  args <- args[-(1)]
  sample_ids = args
  output_dir = "./"
} else {
  # if running via RStudio
  key_file <- NA
  sample_ids =  c("AAAA", "BBBB", "CCCC")
  output_dir =  "../results/"
}

# encryption key
key <- NA
if (file.exists(key_file)) {
  key <- readChar(key_file, file.info(key_file)$size)
} else {
  message("INFO: sample ID encryption key file not found.  Will not encrypt sample IDs.")
}

if (!is.na(key)) {
  encoded_ids = as.character(lapply(sample_ids, Vigenere, key=key))
  output_df <- cbind(sample_ids, encoded_ids)
} else {
  # don't encode the IDS: just output the same IDs twice
  output_df <- cbind(sample_ids, sample_ids)
}

colnames(output_df) <- c("sample_id", "encoded_id")

# setting filename to "" writes to stdout
write.table(output_df, "", quote=F, row.names=F, col.names=T, sep="\t")
?write.table



