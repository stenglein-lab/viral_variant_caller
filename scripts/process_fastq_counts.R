library (tidyverse)
library (openxlsx)
library (egg)

# This script reads in a table of fastq counts and generates some plots and summaries
#
# Mark Stenglein Jan 19, 2021

# ------------------------
# import fastq / bam read count info.   
# ------------------------
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  # TODO: check CLAs
  r_bindir=args[1]
  # -c(1:1) --> all but the first 1 element of the list
  count_file_names = args[-c(1:1)]
  output_directory="./"
} else {
  # if running via RStudio
  r_bindir = "." 
  count_file_names = list.files(path = "../results/fastq_counts", pattern = "*count.txt$", full.names = T)
  output_directory="../results/"
}

read_count_file <- function(file_name) {
  # read the file 
  # wrapped with try here because empty files cause read.delim to fail/error
  df <- try(read.delim(file_name, header = F, stringsAsFactors = F))
  if (!inherits(df, 'try-error')) {

    colnames(df) <- c("sample_id", "count_type", "count")
    return(df)

  } else {
    # a vcf file with no variants return no dataframe
    print(paste0("no count data in file: ", file_name))
    return()
  }
}

# this runs variant_df on all the file names and returns a list of dataframes
# see: https://www.rdocumentation.org/packages/base/versions/3.3.1/topics/lapply?
df_list <- lapply(count_file_names, "read_count_file")

# this rbinds all the dfs together to make one big df
# see: https://stat.ethz.ch/pipermail/r-help/2007-April/130594.html
counts_df <- do.call("rbind", df_list)

counts_df$count_type <- fct_relevel(counts_df$count_type, c("initial",
                                    "post_trimming",
                                    "post_host_filtered",
                                    "refseq_aligned"))

write.table(counts_df, file=paste0(output_directory, "all_read_counts.txt"), sep="\t", row.names=F, col.names=T, quote=F)

normalized_counts_df <- counts_df  %>% 
  pivot_wider(names_from = "count_type", values_from="count")  %>%
  mutate(post_host_filtered = post_host_filtered / initial,
         refseq_aligned = refseq_aligned / initial,
         post_trimming = post_trimming / initial, 
         initial = initial / initial) %>%
  pivot_longer(cols=-sample_id, names_to = "count_type", values_to="count") 

normalized_counts_df$count_type <- fct_relevel(normalized_counts_df$count_type, c("initial",
                                    "post_trimming",
                                    "post_host_filtered",
                                    "refseq_aligned"))


raw_counts_p <- ggplot(counts_df) +
 geom_line(aes(x=count_type, y=count, group=sample_id), size=0.25, linetype="dotted") +
 geom_boxplot(aes(x=count_type, y=count), fill="white", alpha=0.9) +
 theme_bw(base_size = 11) +
 scale_y_log10() + 
 ylab("Reads remaining") + 
 xlab("")

raw_counts_p

norm_counts_p <- ggplot(normalized_counts_df) +
 geom_line(aes(x=count_type, y=count, group=sample_id), size=0.25, linetype="dotted") +
 geom_boxplot(aes(x=count_type, y=count), fill="white", alpha=0.9) +
 theme_bw(base_size = 11) +
 ylab("Fraction reads remaining") + 
 xlab("")

norm_counts_p

p <- ggarrange(raw_counts_p, norm_counts_p, ncol=1)

ggsave(paste0(output_directory, "filtering_plots.pdf"), p, units="in", width=7, height=9)

