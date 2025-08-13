add_column_categorical_major_gene <- function(dataset){
  if (!(all(c("BRCA1r", "BRCA2r", "PALB2r", "ATMr", "CHEK2r", "BARD1r", "RAD51Cr", "RAD51Dr") %in% colnames(dataset)))){
    stop(paste0(deparse(substitute(dataset)), " must include gene fields: BRCA1r, BRCA2r, PALB2r, ATMr, CHEK2r, BARD1r, RAD51Cr, RAD51Dr"))
  }
  dataset$categorical_gene <- rep(NA, nrow(dataset))
  dataset$categorical_gene[which(dataset$RAD51Dr == "P")] <- rep(1, length(which(dataset$RAD51Dr == "P")))
  dataset$categorical_gene[which(dataset$RAD51Cr == "P")] <- rep(2, length(which(dataset$RAD51Cr == "P")))
  dataset$categorical_gene[which(dataset$BARD1r == "P")] <- rep(3, length(which(dataset$BARD1r == "P")))
  dataset$categorical_gene[which(dataset$ATMr == "P")] <- rep(4, length(which(dataset$ATMr == "P")))
  dataset$categorical_gene[which(dataset$CHEK2r == "P")] <- rep(5, length(which(dataset$CHEK2r == "P")))
  dataset$categorical_gene[which(dataset$PALB2r == "P")] <- rep(6, length(which(dataset$PALB2r == "P")))
  dataset$categorical_gene[which(dataset$BRCA2r == "P")] <- rep(7, length(which(dataset$BRCA2r == "P")))
  dataset$categorical_gene[which(dataset$BRCA1r == "P")] <- rep(8, length(which(dataset$BRCA1r == "P")))
  dataset$categorical_gene[which(dataset$BRCA1r == "N" & dataset$BRCA2r == "N" & 
                                             dataset$PALB2r == "N" & dataset$ATMr == "N" &
                                             dataset$CHEK2r == "N" & dataset$BARD1r == "N" &
                                             dataset$RAD51Cr == "N" & dataset$RAD51Dr == "N")] <- rep(0, length(which(dataset$BRCA1r == "N" & dataset$BRCA2r == "N" & dataset$PALB2r == "N" & dataset$ATMr == "N" & dataset$CHEK2r == "N" & dataset$BARD1r == "N" & dataset$RAD51Cr == "N" & dataset$RAD51Dr == "N")))
  dataset$categorical_gene <- as.factor(dataset$categorical_gene)
  
  return(dataset)
}
