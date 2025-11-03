reformat_output_dataset <- function(dataset_path){
  dataset_read <- read.csv(dataset_path)
  dataset_read <- dataset_read[which(dataset_read$Proband == 1),]
  
  dataset_read$Parity[dataset_read$Parity == 3] <- rep(">2", length(which(dataset_read$Parity == 3)))
  
  dataset_read$First_live_birth[!is.na(dataset_read$First_live_birth) & dataset_read$First_live_birth == 18] <- rep("<20", length(which(!is.na(dataset_read$First_live_birth) & dataset_read$First_live_birth == 18)))
  dataset_read$First_live_birth[!is.na(dataset_read$First_live_birth) & dataset_read$First_live_birth == 23] <- rep("20-24", length(which(!is.na(dataset_read$First_live_birth) & dataset_read$First_live_birth == 23)))
  dataset_read$First_live_birth[!is.na(dataset_read$First_live_birth) & dataset_read$First_live_birth == 28] <- rep("25-29", length(which(!is.na(dataset_read$First_live_birth) & dataset_read$First_live_birth == 28)))
  dataset_read$First_live_birth[!is.na(dataset_read$First_live_birth) & dataset_read$First_live_birth == 33] <- rep(">29", length(which(!is.na(dataset_read$First_live_birth) & dataset_read$First_live_birth == 33)))
  
  dataset_read$OC_Use[is.na(dataset_read$OC_Use)] <- rep("N", length(which(is.na(dataset_read$OC_Use))))
  
  dataset_read$MHT_Use[is.na(dataset_read$MHT_Use)] <- rep("N", length(which(is.na(dataset_read$MHT_Use))))
  
  dataset_read$Menopause[is.na(dataset_read$Menopause)] <- rep("N", length(which(is.na(dataset_read$Menopause))))
  dataset_read$Menopause[dataset_read$Menopause == 37] <- rep("<40", length(which(dataset_read$Menopause == 37)))
  dataset_read$Menopause[dataset_read$Menopause == 42] <- rep("40-44", length(which(dataset_read$Menopause == 42)))
  dataset_read$Menopause[dataset_read$Menopause == 47] <- rep("45-49", length(which(dataset_read$Menopause == 47)))
  dataset_read$Menopause[dataset_read$Menopause == 52] <- rep("50-54", length(which(dataset_read$Menopause == 52)))
  dataset_read$Menopause[dataset_read$Menopause == 57] <- rep(">54", length(which(dataset_read$Menopause == 57)))
  
  return(dataset_read)
}
