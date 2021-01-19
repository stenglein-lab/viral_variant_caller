# This script reads in coverage and DVG data and does some pre-processing on it
#
# Mark Stenglein Jan 18, 2021 

# This code block allows separate argument handling when
# run from the command line 
# or run interactively (i.e. via RStudio).  

library (tidyverse)
library (openxlsx)

# if being run from the command line
if (!interactive()) {

  args = commandArgs(trailingOnly=TRUE)

  if (length(args)==0) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
  }

  # snp sift file names passed as command line arguments
  r_bindir=args[1]
  depth_file = args[2]
  # -c(1:2) --> all but the first 2 element of the list
  di_file_names = args[-c(1:2)]
  output_directory="./"
} else {
  # setup argument defaults for use in RStudio interactive environment
  r_bindir="."
  depth_file = "../results/all.depth"
  di_file_names = list.files(path = "../results", pattern = "*di_counts.txt$", full.names = T)
  output_directory="../results/"
}




# ------------------------
# import read depth info.   
#
# Depth of coverage for all viruses at all positions.
# ------------------------
depth_df <- read.delim(depth_file, sep="\t", header=FALSE)
colnames(depth_df) <- c("dataset", "reference_sequence", "position", "depth")

# calculate median depths
median_depths <- depth_df %>% group_by(dataset, reference_sequence) %>% summarize(median_depth = median(depth), .groups="drop")


# ------------------------------
# Import DI-tect quantified DVGs
# ------------------------------
# read DI-tect output of deletion breakpoint split-mapped reads.
#

read_di_file <- function(ditector_file_name) {

  # read the file 
  # wrapped with try here because empty files cause read.delim to fail/error
  df <- try(read.delim(ditector_file_name, header = F, stringsAsFactors = F))
  if (!inherits(df, 'try-error')) {

    colnames(df) <- c("dataset", "type", "length", "bp", "ri", "delta", "ref", "counts", "pct_to_virus") 

    return(df)
  } else {
    # a vcf file with no variants return no dataframe
    print(paste0("no variants for dataset corresponding to vcf file: ", ditector_file_name))
    return()
  }
}

# this runs variant_df on all the file names and returns a list of dataframes
# see: https://www.rdocumentation.org/packages/base/versions/3.3.1/topics/lapply?
df_list <- lapply(di_file_names, "read_di_file")

# this rbinds all the dfs together to make one big df
# see: https://stat.ethz.ch/pipermail/r-help/2007-April/130594.html
di_df <- do.call("rbind", df_list)

typeof(di_df)
di_row <- nrow(di_df)

if (is.null(nrow(di_df))){
  print("no DVGs found in any dataset")  
  quit()
} 
  
  

# parse out virus name from DI-tector output
di_df$reference_sequence <- str_extract(di_df$ref, "(\\S+)\\|") %>% str_replace("\\|","")
di_df$pct_to_virus <- str_extract(di_df$pct_to_virus, "(\\S+)\\|") %>% str_replace("\\|","") %>% str_replace("%", "") 


# filter to create dataframes with only particular types of DVGs 
del_df <- filter(di_df, str_detect(type, "Deletion DVG"))
sb3_df <- filter(di_df, str_detect(type, "3' cb/sb DVG"))
sb5_df <- filter(di_df, str_detect(type, "5' cb/sb DVG"))
ins_df <- filter(di_df, str_detect(type, "Insertion DVG"))

# this function filters unsupported DVGs and 
# normalizes their abundances
process_dvg_df <- function(dvg_df) {
  
  # filter only breakpoints with >N supporting reads
  # note this is a paramter option when running DI-tector from the command line too
  # min_reads_over_breakpoint <- 3
  # dvg_df <- filter(dvg_df, counts >= min_reads_over_breakpoint) 
  
  # merge in depth info
  dvg_df <- left_join(dvg_df, median_depths, by=c("dataset", "reference_sequence"))
  
  # abundance = counts over breakpoint / median depth
  # TODO: just use pct_to_virus instead?
  # TODO: normalize by local depth?
  dvg_df <- dvg_df %>% mutate(rel_abundance = counts / median_depth)  %>% select (-ref, -pct_to_virus)
  
  return(dvg_df)
}

del_df <- process_dvg_df(del_df)
sb3_df <- process_dvg_df(sb3_df)
sb5_df <- process_dvg_df(sb5_df)
ins_df <- process_dvg_df(ins_df)



# this could be used to filter out shorter deletions
# del_df %>% filter(abs(bp-ri) < 100)


make_dvg_table <- function(dvg_df){
  
  # make a table that summarizes DVGs
  dvg_df_summary_table <- dvg_df %>% select(-median_depth) %>% pivot_wider(names_from = dataset, values_from=c(counts, rel_abundance)) %>% arrange(reference_sequence, bp)

  dvg_df_summary_table[is.na(dvg_df_summary_table)] = 0
  
  dvg_df_summary_table <- rename(dvg_df_summary_table, breakpoint_position = bp, reinitiation_position = ri)
  
  return(dvg_df_summary_table)

}

del_summary <- make_dvg_table(del_df)
sb3_summary <- make_dvg_table(sb3_df)
sb5_summary <- make_dvg_table(sb5_df)
ins_summary <- make_dvg_table(ins_df)

# TODO: write to one workbook
# create the workbook
wb <- createWorkbook("dvg_summary.xlsx")

# create a worksheet 
addWorksheet(wb, "deletion")
addWorksheet(wb, "insertion")
addWorksheet(wb, "sb3")
addWorksheet(wb, "sb5")

writeData(wb, "deletion",  del_summary, borders="all")
writeData(wb, "insertion", ins_summary, borders="all")
writeData(wb, "sb3", sb3_summary, borders="all")
writeData(wb, "sb5", sb5_summary, borders="all")


# write out the spreadsheet
saveWorkbook(wb, paste0(output_directory, "dvg_summary.xlsx"), overwrite = TRUE)


