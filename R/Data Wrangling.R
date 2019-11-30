#Purpose of this function grouping is to assist with data wrangling

#' Tibble CSV reader that prints what is being loaded
#'
#' @param path is the location of the file to be loaded.
#' @return tibble object
#' @examples
#' csv_reader("C:/Users/Nathan/Dropbox/@Dana Farber/CyCif/108 - ProjectID_28_21_BTIL_TNBC_RODIG/data_UNET/Parameter_Exploring/original/JENNTNBC01.csv")
#' @export
csv_reader <- function(path) {
  print(path)
  tibble <- read_csv(path)
  return(tibble)
}

#' Load in all CSV files assuming from mcmicro output
#'
#' @param location The path location where all .csv files are stored
#' @return list of tibble objects
#' @examples
#' load_all_csv_files("C:/Users/Nathan/Dropbox/@Dana Farber/CyCif/108 - ProjectID_28_21_BTIL_TNBC_RODIG/data_UNET/Parameter_Exploring/original/")
#' @export
load_all_csv_files <- function(location) {
  #grabs only csv files
  files <- list.files(paste(c(location,'*.csv'),sep='/'))
  #tibble read in csv reader
  output <- lapply(paste(location,files,sep="/"), csv_reader)
  return(output)
}

#' Modify Column Name for CSV file assuming mcmicro HistoCat output
#' Takes the type of data, sample name and removes it from the column name for ease of reading and compatiability with Ajit's IMAAP
#' Adds column for type of data and replacing the sample id column with names
#' default for histocat's columns are: "Cell_[imagename]marker name"
#' rename columns to reflect:
#' ImageId = [imagename]
#' Column = marker name
#' Add column for what kind = Cell
#' assumption is marker ids start at location 3
#'
#' @param data The tible object
#' @return the modified tible object
#' @examples
#' column_cleaner(tibbleobject)
#' @export
column_cleaner <- function(tib) {
  # assumption is marker ids start at location 3
  # using the first id, grab all marker names
  # assumption is location 3 until last column with Type == marker names
  # [WARNING] if features every change need to update

  # grab what type
  Type = names(tib)[3] %>% strsplit('_') %>% unlist()
  Type = Type[1]
  # id markers
  marker_channels <- names(tib %>% select(contains(paste(Type,'_',sep=''))))
  # remove Type from marker name
  tmp <- marker_channels %>% strsplit(paste(Type,'_',sep='')) %>% unlist()
  #keep the end of the split every other to keep
  tmp <- tmp[c(FALSE,TRUE)]
  # id sample name by calculating the longest substring
  samplename<-PTXQC::LCSn(tmp, min_LCS_length = 0)
  # remove sample name from marker name
  tmp <- tmp %>% strsplit(samplename) %>% unlist()
  tmp <- tmp[c(FALSE,TRUE)]
  # replace column names with new marker name
  data.table::setnames(tib,old=marker_channels,new=tmp)
  #replace ImageID with string ID
  tib$ImageId <- samplename
  #add column for type
  tib$Type <- Type
  return(tib)
}

#' Run column cleaner on all tibbles loaded in assuming from mcmicro output
#' See: column_cleaner for explaination
#'
#' @param tibble_list list of tibble objects
#' @return cleaned list of tibble objects
#' @examples
#' clean_all_column_cleaner(tibble_list)
#' @export
clean_all_column_cleaner <- function(tibble_list) {
  #tibble read in csv reader
  output <- lapply(tibble_list, column_cleaner)
  return(output)
}
