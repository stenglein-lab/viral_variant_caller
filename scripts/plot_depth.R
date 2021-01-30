library (tidyverse)
library (openxlsx)
library (pdftools)

# This script reads in a variant table and generates some plots and summaries
#
# Mark Stenglein Oct 24, 2020

# ---------------------------------------------------------------
# import read depth info.   
# Depth of coverage for all reference_sequence at all positions.
# ---------------------------------------------------------------
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
colnames(depth_df) <- c("dataset", "reference_sequence", "position", "depth")

# calculated median depth of (total) coverage for each ref seq in each dataset and store it in a new df
median_depths <- depth_df %>% group_by(dataset, reference_sequence) %>% summarize(median_depth = median(depth))

wb <- createWorkbook("Median_depths.xlsx")
addWorksheet(wb, "median_depth")
writeData(wb, "median_depth", median_depths)
saveWorkbook(wb, "Median_depths.xlsx", overwrite = TRUE)

# higlight coverage below a certain limit
# TODO: Parameterize this
min_depth_highlight <- 30
depth_df <- depth_df %>% mutate(above_highlight = if_else(depth > min_depth_highlight, TRUE, FALSE))

# max position for drawing a red rectangle showing low coverage
max_position = max(depth_df$position)

# ----------------------------------
# calculate average depth in windows
# ----------------------------------
# %/% is the integer division operator
window_size = 10
depth_df <- depth_df %>% mutate (window = position %/% window_size)

# calculate average coverage depth in each window
df_windowed <- depth_df %>% 
  group_by(dataset, reference_sequence, window)  %>% 
  summarize(depth = mean(depth), .groups = "drop") %>% 
  mutate(position = (window*window_size) + 1) %>% 
  ungroup()

##now plot coverage data on multiple pdf pages

# these are the datasets
datasets <- depth_df %>% group_by(dataset) %>% summarize(.groups="drop") %>% pull()

plot_some_datasets <- function(datasets){

  datasets_per_page <- 12

  page_number <- 1

  # iterate through the datasets, doing up to 12 per page
  for (i in seq(1, length(datasets), datasets_per_page)) {

    # plots_per_page at a time
    subset_datasets <- datasets[i:(i+(datasets_per_page-1))]

    # output to console which ones we're doing
    # print(paste0(subset_datasets))

    pdf_name <- paste0("coverage_plot_page_", page_number, ".pdf")

    # generate & print plot
    plot_datasets(subset_datasets, pdf_name)

    page_number = page_number + 1
  }

}

plot_datasets <- function(dataset_names, pdf_name){
  
  # subset the main dataframes to get the data just for these reference_sequence
  subset_df <- df_windowed %>% filter(dataset %in% dataset_names)

  # convert any depth of 0 into 1, since plotting on a log10 y scale...
  subset_df <- subset_df %>% mutate(depth = if_else(depth == 0, 1, depth))
  
  # output the reference_sequence names to the console
  # print(paste0(dataset_names))
  
  # p returned here is a ggplot object
  p <- ggplot(subset_df) + 
    geom_line(aes(x=position, y=depth), size=0.5) +
    geom_area(aes(x=position, y=depth), fill="lightgrey") +
    annotate("rect", xmin=0, xmax=max_position, ymin=1, ymax=min_depth_highlight, alpha=0.05, fill="red") +
    scale_color_manual(values = c("red", "black")) +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    scale_y_log10() +
    xlab("genome position (nt)") +
    ylab (paste0("coverage depth\n", "coverage below ", min_depth_highlight, " in red"))+
    facet_grid(dataset~reference_sequence, scales="free", space="free_x")
  
  ggsave(pdf_name, p, height=10.5, width=7.5, units="in")

    # print p will make the plots appear on viewer
    # print(p)
}

plot_some_datasets(datasets)

pdf_combine(c(list.files(pattern="pdf$")), output="coverage_plot.pdf")
