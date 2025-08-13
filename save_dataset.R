save_dataset <- function(dataset, path = "~/OneDrive/DEV-Models-Batching-4 copy/data/", filename = "default"){
  if (filename == "default"){
    filename_location <- paste(path, deparse(substitute(dataset)), ".csv", sep = "")
  } else {
    filename_location <- paste(path, filename, ".csv", sep = "")
  }
  write.csv(dataset, filename_location, quote = FALSE, row.names = FALSE)
}
