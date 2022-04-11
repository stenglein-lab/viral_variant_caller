# This function pulls out metadata (sample date) from sample IDs
#
# Mark Stenglein 3/29/2022

# this function pulls out date and returns it in a GISAID-compatible format
extract_date_from_dataset_id <- function(dataset_id){
  
  # pull out date in GISAID format
  # parse out date metadata from sample IDs
  # expected to be in the sample ID (dataset name)
  # in the format: YYYYMMDD or YYYY
  date_regex <- "_(20[0-9]{2})([0-9]{2})([0-9]{2})"
  date_fields <- str_match(dataset_id, date_regex)
  
  whole_match = date_fields[,1]
  year        = date_fields[,2]
  month       = date_fields[,3]
  day         = date_fields[,4]
  if (is.na(year)) {
    
    message(paste0("WARNING: could not parse year out of sample ID: ", dataset_id))
    message(       "         this will result in invalid GISAID submission data.")
    gisaid_date <- NA_character_
    
  } else if (!is.na(year) & !is.na(month) & !is.na(day)) {
    
    # the field covv_collection_date should be in the format:
    # 2021-03-31
    gisaid_date <- paste0(year, "-", month, "-", day)
    
  }
  else {
    
    # or if no month/day info, just the year
    # 2021
    gisaid_date <- year
    
  }
  
  gisaid_date
}

