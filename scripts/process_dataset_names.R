
# parse out metadata from dataset names
# example dataset ID: N19_2_20_INF_pass2_dpi3_flask1_tubeid54_R1.variant_alleles.txt
process_dataset_names <- function(df) {

  df$passage <- str_match(df$dataset, "_pass(\\d+)_")[,2]
  df$flask   <- str_match(df$dataset, "_flask(\\d+)_")[,2]
  df$tube    <- str_match(df$dataset, "_tubeid(\\d+)_")[,2]
  df$lot     <- str_match(df$dataset, "(N19_[12])")[,2]

  # make a nicer dataset name
  # df <- df %>% mutate(nice_dataset_name = paste0(lot, "_P", passage, "_flask_", flask, "_(", tube,")")) %>% arrange(passage,flask,tube) %>% select(-passage, -tube, -flask)
  df <- df %>% mutate(nice_dataset_name = paste0(lot, "_P", passage, "_flask_", flask, "_(", tube,")")) %>% arrange(passage,flask,tube) 
  
  return(df)

}
