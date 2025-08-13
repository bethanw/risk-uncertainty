categorical_to_binary_major_gene <- function(dataset){
  # takes dataset with categorical_gene field, add major gene fields
  if ("categorical_gene" %in% colnames(dataset)){
    dataset$BRCA1r <- rep(NA, nrow(dataset))
    dataset$BRCA2r <- rep(NA, nrow(dataset))
    dataset$PALB2r <- rep(NA, nrow(dataset))
    dataset$CHEK2r <- rep(NA, nrow(dataset))
    dataset$ATMr <- rep(NA, nrow(dataset))
    dataset$BARD1r <- rep(NA, nrow(dataset))
    dataset$RAD51Cr <- rep(NA, nrow(dataset))
    dataset$RAD51Dr <- rep(NA, nrow(dataset))
    
    dataset[!is.na(dataset$categorical_gene), c("BRCA1r", "BRCA2r", "PALB2r", "CHEK2r", "ATMr", "BARD1r", "RAD51Cr", "RAD51Dr")] <- rep("N", length(which(!is.na(dataset$categorical_gene)))*8)
    
    dataset$BRCA1r[dataset$categorical_gene == 8] <- "P"
    dataset$BRCA2r[dataset$categorical_gene == 7] <- "P"
    dataset$PALB2r[dataset$categorical_gene == 6] <- "P"
    dataset$CHEK2r[dataset$categorical_gene == 5] <- "P"
    dataset$ATMr[dataset$categorical_gene == 4] <- "P"
    dataset$BARD1r[dataset$categorical_gene == 3] <- "P"
    dataset$RAD51Cr[dataset$categorical_gene == 2] <- "P"
    dataset$RAD51Dr[dataset$categorical_gene == 1] <- "P"
    
    return(dataset)
  }
  else{
    stop("categorical_gene column not present")
  }
}
