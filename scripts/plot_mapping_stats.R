library (tidyverse)

# This script reads in a table of read-mapping statistics and generates plots
#
# Mark Stenglein Mar 23, 2021

# ---------------------------------------------------------------
# import mapping stats info.   
# Depth of coverage for all reference_sequence at all positions.
# ---------------------------------------------------------------
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  # TODO: check CLAs
  r_bindir=args[1]
  stats_file_name=args[2]
  output_dir = "./"
} else {
  # if running via RStudio
  r_bindir = "." 
  stats_file_name = "../results/all.mapping_stats"
  output_dir = "../results/"
}

stats_df <- read.delim(stats_file_name, sep="\t", header=FALSE)

# from samtools stats output
# Insert sizes. Use `grep ^IS | cut -f 2-` to extract this part. The columns are: insert size, pairs total, inward oriented pairs, outward oriented pairs, other pairs
colnames(stats_df) <- c("dataset", "insert_size", "pairs_total", "inward_pairs", "outward_pairs", "other_pairs")

# ----------------------------------------------------------
# calculate average read pair counts in insert size windows
# ----------------------------------------------------------
# %/% is the integer division operator
window_size = 25
stats_df <- stats_df %>% mutate (window = insert_size %/% window_size)

# calculate average coverage depth in each window
df_windowed <- stats_df %>%
  group_by(dataset, window)  %>%
  summarize(pairs_sum = sum(pairs_total), .groups = "drop") %>%
  mutate(insert_size = (window*window_size) + 1) %>%
  ungroup()

# p returned here is a ggplot object
p <- ggplot(df_windowed) + 
	 geom_col(aes(x=insert_size, y=pairs_sum), fill="darkslateblue", color="black", size=0.25) + 
    theme_bw(base_size = 10) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    scale_y_log10() +
	 xlim(c(0,600)) + 
    xlab("insert size (nt)") +
    ylab ("mapped read pairs with that insert size") +
	 facet_wrap(~dataset) +
    theme(strip.text.y = element_text(angle = 0)) 
  
ggsave(paste0(output_dir, "mapping_stats_plot.pdf"), p, height=7.5, width=10.5, units="in")

