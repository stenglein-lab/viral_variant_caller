library(tidyverse)
library(readxl)
library(openxlsx)

# deal with expected input
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  # TODO: check CLAs
  script_dir = args [1]
  depth_xls = args[2]
  completeness_file = args[3]
  pangolin_file = args[4]
  minimum_fraction_called = args[5]
  output_prefix = args[6]
  obfuscated_ids_file = args[7]
  output_dir="./"
} else {
  # if running via RStudio
  script_dir = "./"
  depth_xls = "../results/Average_depths.xlsx"
  completeness_file = "../results/all_consensus_completeness.txt"
  pangolin_file = "../results/pangolin_lineage_report.csv"
  minimum_fraction_called = 0.95
  output_prefix = paste0(Sys.Date(), "_")
  obfuscated_ids_file = "../results/obfuscated_sample_ids.txt"
  output_dir="../results/"
}

# make sure a number
minimum_fraction_called <- as.numeric(minimum_fraction_called)


# --------------------
# Read in data files
# --------------------

# average depth for each dataset
depth_df <- read_excel(depth_xls)

# fraction refseq completeness for each dataset
completeness_df <- read.delim(completeness_file, sep="\t", header=F)
colnames(completeness_df) <- c("dataset", "fraction_called")

# make names nicer
completeness_df$dataset <- str_replace(completeness_df$dataset, "_consensus.fasta", "")

# pull date out of sample name
# this sourced file contains a function named extract_date_from_dataset_id 
# that pulls out dates from sample IDs
source (paste0(script_dir, "/extract_metadata_from_sample_id.R"))

# parse date from dataset ID 
completeness_df$collection_date <- sapply(completeness_df$dataset, extract_date_from_dataset_id)

# obfuscated sample IDs
obfuscated_ids <- read.delim(obfuscated_ids_file, sep = "\t", header = T)

# pangolin_lineage_report
pangolin_df <- read.delim(pangolin_file, sep=",", header=T)

# parse out N content from Pangolin
pangolin_df$pangolin_fraction_N <-  as.numeric(str_match(pangolin_df$note, "N_content:([01].\\d+)")[,2])

# merge tables
df <- left_join(completeness_df, depth_df)
df <- left_join(df, pangolin_df, by=c("dataset" = "taxon"))
df <- left_join(df, obfuscated_ids, by=c("dataset" = "original_id"))

# relationship between depth of coverage and fraction of genome called
ggplot(df) +
  geom_point(aes(x=mean_depth, y=fraction_called),
             shape=21, color="black", fill="cornflowerblue",
             size=2.5, stroke=0.25, alpha=0.75) +
  xlab("Mean depth of coverage") +
  ylab("Fraction of genome called (not N)") + 
  scale_x_log10() +
  ggtitle("Relationship between coverage and fraction of genome recovered", 
          subtitle="CDPHE SARS-CoV-2 sequencing, Run 1: 10/25/2021") +
  theme_classic(base_size=16)

ggsave(paste0(output_dir, output_prefix, "fraction_called_vs_coverage.pdf"), height=7.5, width=10, units="in")

# same but non-log scale
ggplot(df) +
  geom_point(aes(x=mean_depth, y=fraction_called),
             shape=21, color="black", fill="cornflowerblue",
             size=2.5, stroke=0.25, alpha=0.75) +
  xlab("Mean depth of coverage") +
  ylab("Fraction of genome called (not N)") + 
  xlim(c(0,500)) +
  ggtitle("Relationship between coverage and fraction of genome recovered", 
          subtitle="CDPHE SARS-CoV-2 sequencing, Run 1: 10/25/2021") +
  theme_classic(base_size=16)
  

ggsave(paste0(output_dir, output_prefix, "fraction_called_vs_coverage_non_log.pdf"), height=7.5, width=10, units="in")

# ggplot(df) +
  # geom_histogram(aes(x=fraction_called), bins=50,
             # color="black", fill="coral3",
             # size=0.25) +
  # xlab("Fraction of genome called (not N)") +
  # ylab("Number samples") +
  # ggtitle("Number of samples with particular levels of completeness",
          # subtitle="CDPHE SARS-CoV-2 sequencing, Run 1: 10/25/2021") +
  # theme_classic(base_size=16)

# ggsave("completeness_histogram.pdf", height=7.5, width=10, units="in")



ggplot(df) +
  stat_ecdf(aes(y=fraction_called), geom = "step", color="blue") +
  xlab("Fraction of  samples") +
  ylab("Fraction of genome called (not N)") +
  ylim(c(0,1)) + 
  annotate("segment", x = 0, xend = 1, y = minimum_fraction_called, yend = minimum_fraction_called,
           colour = "red", size=0.5, alpha=0.5, linetype=2) +
  ggtitle("Fraction of samples with particular levels of completeness",
          subtitle="CDPHE SARS-CoV-2 sequencing, Run 1: 10/25/2021") +
  theme_bw(base_size=16)

ggsave(paste0(output_dir, output_prefix, "completeness_cdf.pdf"), height=7.5, width=10, units="in")


# make an Excel report  for all samples
report_df <- df %>% 
  arrange(desc(fraction_called), desc(collection_date)) %>% 
  select(dataset, fraction_called, collection_date,
         reference_sequence, 
         median_depth, mean_depth, 
         lineage, ambiguity_score, 
         scorpio_call, scorpio_support, 
	 version, pangolin_version, 
         constellation_version, is_designated, qc_status, qc_notes, note, obfuscated_id)  %>%
  rename(gisaid_sequence_id = obfuscated_id)


# create the workbook
wb <- createWorkbook(paste0(output_prefix, "dataset_summary.xlsx"))

# create a worksheet 
addWorksheet(wb, "summary")
writeData(wb, "summary", report_df, borders="all")

# add some style to this workbook
all_cell_style <- createStyle( 
  fontName = "Helvetica",
  fontSize = 11,
  numFmt = "0.00",
  border = "TopBottomLeftRight",
  borderColour = getOption("openxlsx.borderColour", "grey"),
  borderStyle = getOption("openxlsx.borderStyle", "thin"),
  halign = "left"
)

addStyle(wb=wb, sheet="summary", 
        style=all_cell_style,
        cols = 1:ncol(report_df)+1,
        rows = 1:nrow(report_df)+1,
        gridExpand = T
)

# integer numbers for some columns
integer_num_style <- createStyle( 
  fontName = "Helvetica",
  fontSize = 11,
  numFmt = "0",
  border = "TopBottomLeftRight",
  borderColour = getOption("openxlsx.borderColour", "grey"),
  borderStyle = getOption("openxlsx.borderStyle", "thin"),
  halign = "left"
)

addStyle(wb=wb, sheet="summary", 
        style=integer_num_style,
        cols = 5:6,
        rows = 1:nrow(report_df)+1,
        gridExpand = T
)

# row headers
row_header_style <- createStyle( 
  fontName = "Helvetica",
  fontSize = 11,
  numFmt = "TEXT",
  textDecoration = "bold",
  wrapText = TRUE
)

# style header row
addStyle(wb=wb, sheet="summary", 
         style=row_header_style,
         cols = 1:ncol(report_df),
         rows = 1:1,
         gridExpand = T
)


# column headers
col_header_style <- createStyle( 
  fontName = "Helvetica",
  fontSize = 11,
  numFmt = "TEXT"
)

# style header cols
addStyle(wb=wb, sheet="summary", 
         style=col_header_style,
         cols = 1:1,
         rows = 1:nrow(report_df)+1,
         gridExpand = T
)

# custom formatting style
# green background color
custom_annotation_style <- createStyle( 
  bgFill = "#E5FFCC"
)

# add conditional formatting to highlight any rows with custom annotations
conditionalFormatting(wb=wb, sheet="summary", 
                      type = "expression",
                      cols = 2:2,
                      rows = 1:nrow(report_df)+1,
                      style = custom_annotation_style,
                      rule = ">=0.95"
)

# add an Excel filter
addFilter(wb=wb, sheet="summary", 
          row= 1:1, 
          cols = 1:ncol(report_df))

# set auto column widths for header columns
setColWidths(
  wb, sheet="summary",
  cols = 1:ncol(report_df),
  widths = "auto"
)

# write out the spreadsheet
saveWorkbook(wb, paste0(output_dir, output_prefix, "dataset_summary.xlsx"), overwrite = TRUE)





# plot concordance between Pangolin's fraction of N bases (N_content field) and ours
# (1:1 concordance)
ggplot(df) +
  geom_point(aes(x=fraction_called, y=(1-pangolin_fraction_N))) +
  theme_bw() +
  coord_fixed() +
  xlim(c(0,1)) +
  ylim(c(0,1))

