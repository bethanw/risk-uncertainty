save_dataset <- function(dataset, filepath = "default"){
  if (filename == "default"){
    filename_location <- paste(deparse(substitute(dataset)), ".csv", sep = "")
  } else {
    filename_location <- filepath
  }
  write.csv(dataset, filename_location, quote = FALSE, row.names = FALSE)
}
