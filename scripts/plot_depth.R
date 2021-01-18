library (tidyverse)
library (openxlsx)

# This script reads in a variant table and generates some plots and summaries
#
# Mark Stenglein Oct 24, 2020

# ------------------------
# import read depth info.   
#
# Depth of coverage for all viruses at all positions.
#
# all.depth was created by running this command (in mdstengl@cctsi-104:~/datasets/bunyas/all_sequences_1_11_19)
# individual depth files were created using samtools depth
#
# cat *.depth > all.depth
# cp all.depth ../bunya_dvg_analyses/data
#
# ------------------------
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  # TODO: check CLAs
  r_bindir=args[1]
  depth_file_name=args[2]
} else {
  # if running via RStudio
  r_bindir = "." 
  depth_file_name = "../results/all.depth"
}

depth_df <- read.delim(depth_file_name, sep="\t", header=FALSE)
colnames(depth_df) <- c("dataset", "virus", "position", "depth")

# rename datasets 

source(paste0(r_bindir, "/process_dataset_names.R"))
depth_df <- process_dataset_names(depth_df)

# calculated median depth of (total) coverage for each virus in each dataset and store it in a new df
median_depths <- depth_df %>% group_by(dataset, virus) %>% summarize(median_depth = median(depth))

wb <- createWorkbook("Median_depths.xlsx")
addWorksheet(wb, "median_depth")
writeData(wb, "median_depth", median_depths)
saveWorkbook(wb, "Median_depths.xlsx", overwrite = TRUE)

# higlight coverage below a certain limit
# TODO: Parameterize this
min_depth_highlight <- 100
depth_df <- depth_df %>% mutate(above_highlight = if_else(depth > min_depth_highlight, TRUE, FALSE))

# max position for drawing a red rectangle showing low coverage
max_position = max(depth_df$position)

# plot depth for all viruses, all datasets
ggplot(depth_df) +
  geom_line(aes(x=position, y=depth), size = 0.5) +
  # see: https://stackoverflow.com/questions/17521438/geom-rect-and-alpha-does-this-work-with-hard-coded-values
  annotate("rect", xmin=0, xmax=max_position, ymin=1, ymax=min_depth_highlight, alpha=0.05, fill="red") +
  scale_color_manual(values = c("red", "black")) +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        strip.text.y = element_text(angle=0)) +
  # facet_grid(nice_dataset_name ~ virus, scales="free_y") +
  facet_grid(dataset ~ virus, scales="free_y") +
  scale_y_log10() +
  xlab ("position in genome (nt)") +
  ylab (paste0("coverage depth\n", "coverage below ", min_depth_highlight, " in red"))

depth_df

ggsave("coverage_plot.pdf", units="in", width=10, height = 6.5)

ggplot(depth_df) +
  geom_boxplot(aes(x=dataset, y=depth)) +
  # annotate("rect", xmin=0, xmax=max_position, ymin=1, ymax=min_depth_highlight, alpha=0.05, fill="red") +
  # scale_color_manual(values = c("red", "black")) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ggtitle("Even high median coverage across datasets\n2020-11-27") +
  scale_y_log10() +
  ylab (paste0("boxplots of coverage depth\nat all genome positions in indicated datasets"))

ggsave("coverage_plot_2.pdf", units="in", width=10, height = 6.5)
