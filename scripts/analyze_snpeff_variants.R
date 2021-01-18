library(tidyverse)
library(openxlsx)
library(readxl)

# This script reads in files output by snpsift that includes SNV and indel variants 
# and tabulates these variant calls in all datasets
#

# This code block allows separate argument handling when
# run from the command line 
# or run interactively (i.e. via RStudio).  

# if being run from the command line
if (!interactive()) {
  
  args = commandArgs(trailingOnly=TRUE)
  
  if (length(args)==0) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
  } 
  
  # snp sift file names passed as command line arguments
  r_bindir=args[1]
  min_allele_freq = args[2]
  depth_file = args[3]
  min_depth_to_call_variant = args[4]
  # -c(1:4) --> all but the first four element of the list
  snp_sifts = args[-c(1:4)]
  output_directory="./"
} else {
  # setup argument defaults for use in RStudio interactive environment
  r_bindir="."
  min_allele_freq = 0.03
  depth_file = "../results/all.depth"
  min_depth_to_call_variant = 40
  snp_sifts = list.files(path = "../results", pattern = "*.variants.tsv$", full.names = T)
  output_directory="../results/"
}
  
# read in coverage depth info for all datasets
depth_df <- read.delim(depth_file, sep="\t", header=FALSE)
colnames(depth_df) <- c("sample_id", "reference_sequence", "position", "depth")

# replace NA depth values with 0, since no depth called means 0 depth of coverage.
depth_df <- depth_df %>% mutate(depth = if_else(is.na(depth), 0L, depth))

# calculate median depths
# median_depths <- depth_df %>% group_by(sample_id, reference_sequence) %>% summarize(median_depth = median(depth), .groups="drop")

# this function  takes one tsv file created by snpsift
# and returns a dataframe with information about the variants
read_variant_file <- function(snp_sift_file_name) {
  
  # read the file 
  # wrapped with try here because empty files cause read.delim to fail/error
  df <- try(read.delim(snp_sift_file_name, header = F, stringsAsFactors = F))
  if (!inherits(df, 'try-error')) {
    
    names(df) <- c("sample_id", "reference_sequence", "position", "reference_base", 
                   "variant_base", "frequency", "depth", "strand_bias", "indel", 
                   "effect", "impact", "gene", "variant")
    
    return(df)
  } else {
    # a vcf file with no variants return no dataframe
    print(paste0("no variants for dataset corresponding to vcf file: ", snp_sift_file_name))
    return()
  }
}

# this runs variant_df on all the file names and returns a list of dataframes
# see: https://www.rdocumentation.org/packages/base/versions/3.3.1/topics/lapply?
df_list <- lapply(snp_sifts, "read_variant_file")

# this rbinds all the dfs together to make one big df
# see: https://stat.ethz.ch/pipermail/r-help/2007-April/130594.html
df <- do.call("rbind", df_list)

# TODO: report / include datasets with no variants in the final table 
# datasets_with_variants    <- df %>% group_by(sample_id) %>% summarize(.groups = "drop") 
# datasets_with_no_variants <- median_depths[!(median_depths$sample_id %in% datasets_with_variants$sample_id),]

# make sure numeric info is numeric
df$position <- as.numeric(as.character(df$position))
df$frequency <- as.numeric(as.character(df$frequency))
df$depth <- as.numeric(as.character(df$depth))


# make a single-column key that uniquely identifies each particular variant 
df <- df %>% 
  mutate(variant_key = paste0(reference_sequence, position, reference_base, variant_base))

# this lookup table will keep track of the separated variant_key values
variant_key_map <- df %>% 
  group_by(variant_key, reference_sequence, position, reference_base, variant_base) %>%
  summarize(.groups="drop")

# create a df that describes variants, with all frequency/depth information stripped out
# this df will be linked later to the frequency df by the variant_key
variants <- df %>% 
  group_by(variant_key, indel, effect, impact, gene, 
           variant) %>% 
  mutate(number_occurences = n()) %>%
  filter(row_number() == 1) %>%
  select(-c(sample_id, frequency, depth, strand_bias)) %>%
  ungroup()

# convert 3 letter AA code to 1 letter code
# see: https://stackoverflow.com/questions/12760271/how-do-i-convert-the-three-letter-amino-acid-codes-to-one-letter-code-with-pytho/12790512
AA_code_3_to_1 <- function(string){
  code3 <- c("Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", 
             "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", 
             "Tyr", "Val")
  code1 <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", 
             "M", "F", "P", "S", "T", "W", "Y", "V")
  
   for (i in 1:length(code3))
   {
     string <- gsub(code3[i],code1[i],string,ignore.case=TRUE)
   }
  return (string)
}

# convert 3 letter AA codes to 1-letter codes
variants$variant <- lapply(variants$variant, AA_code_3_to_1)


# 
# create clearer descriptions from SNPSift annotations
#

# get rid of "p." at the beginning 
variants$variant <- str_replace(variants$variant, "^p.", "")
# get rid of ",." at the end 
variants$variant <- str_replace(variants$variant, ",.$", "")


# get rid of the text ",intragenic_variant" at the end
variants$effect <- 
  str_replace(variants$effect, ",intragenic_variant$", "")

# how common are variants?
# ggplot(variants) + geom_histogram(aes(x=number_occurences), bins=60) + theme_bw() + scale_y_log10()


# create df that only contains info about the variant frequencies and their depths 
# we will ignore strand bias info for now since for amplicon data (which is an expected data type)
# strand bias is expected 
frequencies <- df %>% select(variant_key, sample_id, frequency, depth)

# ignore variants with no instances of frequency > min_allele_freq in any of the datasets
freq_over_cutoff <- frequencies %>% 
  group_by(variant_key) %>% 
  mutate(max_frequency = max(frequency)) %>% 
  filter(max_frequency > min_allele_freq) %>% 
  select(-max_frequency) %>% 
  ungroup()

# this complete() call fills out the dataset with all possible
# combinations of sample_id and variant_key 
# this will allow us to assign a variant frequency to each variant
# in each dataset
freq_complete <- freq_over_cutoff %>% complete(sample_id, variant_key)

# pull back in refseq and position data using info cached in variant_key
freq_complete <- left_join(freq_complete, variant_key_map, by="variant_key")

# bring in depth info, even for positions with no variant info 
freq_complete_with_depth <- left_join(freq_complete, depth_df, by=c("sample_id", "reference_sequence", "position"))

# after this left join, there are 2 depths:
# 1) depth.y, from the original depth_df, which is the depth from bwa/samtools tabulation of coverage depth
# 2) depty.x, the depth from lofreq (from the vcf/snpsift file).  This is the depth from variant calling
#
# These numbers are expected to be very similar but not necessarily identical 
#

# how similar?
freq_complete_with_depth %>% filter(!is.na(depth.x)) %>% ggplot() +
  geom_point(aes(x=depth.x, y=depth.y)) + theme_bw() +
  xlab("Depth from lofreq") +
  ylab("Depth from bwa") + 
  scale_y_log10() +
  scale_x_log10()

# plotting confirms they are similar enough to drop one: we'll drop the one from lofreq
# because we actually only have depth from lofreq for positions with called variants
# but the whole point is that we want info for all positions, whether or not a variant was called
# so that we can distinguish variants that were not called because there wasn't enough coverage to
# call a variant and variants that were not called because there was enough data but no evidence of 
# a variant
freq_complete_with_depth <- freq_complete_with_depth %>% 
  rename(depth = depth.y) %>% 
  select(-depth.x, -reference_sequence, -position, -reference_base, -variant_base)

# assign a frequency of 0 to variants that were not detected by lofreq, despite enough coverage
# and a frequency of NA to variants that were not detected but don't have enough coverage to confidently
# say that there is no variant
# leave the frequency as is if already called
freq_complete_with_depth <- freq_complete_with_depth %>% 
  mutate(frequency = case_when(
    !is.na(frequency) ~ frequency,
    (depth > min_depth_to_call_variant) ~ 0.0,
    TRUE ~ NA_real_
  ))

# calculate median depth 
freq_complete_with_depth %>% group_by(sample_id, variant_key) %>% mutate(median_depth = median(depth), .groups="drop")

# subset variants that are actually in the filtered freq table
variants_in_freq <- freq_complete_with_depth %>%
  group_by(variant_key) %>% 
  summarize(.groups="drop")

variants_to_report <- variants %>% filter(variant_key %in% variants_in_freq$variant_key)

# now join together the frequency table with the variant info table
# this will create a dataframe with all the info necessary for reporting
df_to_report <- left_join(freq_complete_with_depth, variants_to_report, 
                          by = "variant_key")


# TODO: flag variants from datasets with too low (below some limit?) median coverage

# a summary table to output
df_wide_enough_data <- 
  # filter(df_to_report, median_depth > min_depth_to_call_variant) %>% 
  df_to_report %>% 
  pivot_wider(id_cols = 
                c(reference_sequence, 
                  position, 
                  gene, 
                  indel,
                  variant, 
                  reference_base, 
                  variant_base, 
                  effect), 
              names_from=sample_id, 
              values_from=frequency,
              names_sort = T) %>%
  arrange(position)

# these values describe the # of header rows and columns in the big wide data table
num_header_row = 1 
num_header_col = 8 

# TODO: Merge in dataset metadata and add to top of table

# ----------------------------
# write data to an Excel file
# using openxlsx functions
# ----------------------------

# create the workbook
wb <- createWorkbook("variant_summary.xlsx")

# create a worksheet 
addWorksheet(wb, "variants")
writeData(wb, "variants", df_wide_enough_data, borders="all")


# add some style to this workbook

all_cell_style <- createStyle( 
  fontName = "Helvetica",
  fontSize = 11,
  numFmt = "PERCENTAGE",
  border = "TopBottomLeftRight",
  borderColour = getOption("openxlsx.borderColour", "grey"),
  borderStyle = getOption("openxlsx.borderStyle", "thin"),
  halign = "left"
)

addStyle(wb=wb, sheet="variants", 
        style=all_cell_style,
        cols = num_header_col+1:ncol(df_wide_enough_data),
        rows = 2:nrow(df_wide_enough_data),
        gridExpand = T
)

row_header_style <- createStyle( 
  fontName = "Helvetica",
  fontSize = 11,
  numFmt = "TEXT",
  textDecoration = "bold",
  wrapText = TRUE
)

col_header_style <- createStyle( 
  fontName = "Helvetica",
  fontSize = 11,
  numFmt = "TEXT"
)

# style header cols
addStyle(wb=wb, sheet="variants", 
        style=col_header_style,
        cols = 1:num_header_col,
        rows = 1:nrow(df_wide_enough_data),
        gridExpand = T
)

# style header row
addStyle(wb=wb, sheet="variants", 
        style=row_header_style,
        cols = 1:nrow(df_wide_enough_data),
        rows = 1:num_header_row,
        gridExpand = T
)

# freeze top and left panes
freezePane(wb, "variants", firstActiveRow = num_header_row+1, firstActiveCol = num_header_col+1)

# add conditional formatting
# ?conditionalFormatting
conditionalFormatting(wb=wb, sheet="variants", 
                      "colourScale",
                      cols = num_header_col+1:ncol(df_wide_enough_data), rows = num_header_row+1:nrow(df_wide_enough_data),
                      style = c("white", "green"),
                      rule = c(0, 1),
                      type = "colourScale"
)

# add an Excel filter
addFilter(wb=wb, sheet="variants", 
          row=1:num_header_row, 
          cols = 1:num_header_col)

# set auto column widths for header columns
setColWidths(
  wb, sheet="variants",
  cols = 1:num_header_col,
  widths = "auto"
)



# write out the spreadsheet
saveWorkbook(wb, "variant_summary.xlsx", overwrite = TRUE)



