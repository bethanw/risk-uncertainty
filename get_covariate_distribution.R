# Impute covariates for BOADICEA risk prediction

# Parameters:

# data_stack_for_mice_path: file containing dataset, csv file in BOADICEA format, to use as a reference to impute data
# gene_list: List of major genes in model
# n_chain: Number of mice chains
# data_stack_for_mice_colclasses: Column classes of reference dataset. Needs to be in the same order as reference dataset.
# variable_rfs: "all" if all risk factors are to be imputed. if only a subset of risk factors are to be imputed, set to a vector of the risk factors to impute
# blank_row_path: File containing blank row
# flb_path: File containing "Age at first live birth" data. Separate, because this is imputed after "Parity" is imputed.
# fams_path: File containing pedigree information for the reference data
# given_rfs: List of given risk factors
# rf_values: Values of given risk factors
# quickpred_no: Number of individuals to use if the "quickpred" feature of mice is used. Left as NA otherwise.
# quickpred_cor: Minimum correlation to use if the "quickpred" feature of mice is used. Left as NA otherwise.
# n_repeats: Number of imputations
# n_iter: Number of iterations in mice
# indiv_yob: Year of birth for the individual
# indiv_age: Age of the individual
# censoring_age: Age to calculate BOADICEA risk up to
# bc_prs_alpha: Alpha value for the PGS
# save_filepath: File name to save results to

get_covariate_distribution <- function(data_stack_for_mice_path, 
                                       gene_list = c("BRCA1r", "BRCA2r", "PALB2r", "CHEK2r", "ATMr", "RAD51Cr", "RAD51Dr", "BARD1r"), 
                                       n_chain = 1,
                                       data_stack_for_mice_colclasses, 
                                       variable_rfs = "all", 
                                       blank_row_path, 
                                       flb_path, 
                                       fams_path, 
                                       given_rfs = c(), 
                                       rf_values = c(), 
                                       quickpred_no = NA, 
                                       quickpred_cor = NA, 
                                       n_repeats = 1000, 
                                       n_iter = 10, 
                                       indiv_yob, 
                                       indiv_age = 40, 
                                       censoring_age = 50, 
                                       bc_prs_alpha = 0.44, 
                                       save_filepath = "default") {
  
  # Read in "blank_row" with same structure as the reference dataset
  blank_row <- read.csv(blank_row_path, colClasses = rep("character", 51))
  
  # Add in columns for imputing family history
  blank_row$mother.breast_v0 <- NA
  blank_row$sibling.breast_v0 <- NA
  blank_row$sibling.prostate_v0 <- NA
  blank_row$father.prostate_v0 <- NA
  
  # Put given values into "blank_row"
  for (i in 1:length(given_rfs)){
    blank_row[,given_rfs[i]] <- rf_values[i]
  }
  
  # Make stack of blank rows for imputation
  blank_row_stack <- blank_row
  for (i in 2:n_repeats){
    blank_row_stack[i,] <- blank_row
  }
  
  # Make individual ID and family ID different for each repeat
  blank_row_stack$IndivID <- 1:n_repeats
  blank_row_stack$FamID <- 1:n_repeats
  
  # Read in reference dataset
  data_stack_for_mice <- read.csv(data_stack_for_mice_path)
  data_stack_for_mice <- data_stack_for_mice[,colnames(data_stack_for_mice) != "X"]
  
  # Add column for "First_live_birth" if it is not present (here we are assuming that all other relevant columns are present)
  if (!("First_live_birth" %in% colnames(blank_row_stack))){
    blank_row_stack$First_live_birth <- rep(NA, nrow(blank_row_stack))
  }
  blank_stack_flb <- blank_row_stack$First_live_birth
  
  # Adding column for "categorical_gene" if needed
  if (!("categorical_gene" %in% colnames(blank_row_stack))){
    blank_row_stack$categorical_gene <- rep(NA, nrow(blank_row_stack))
  }
  blank_row_stack_mice_cols <- blank_row_stack[,colnames(data_stack_for_mice)]
  stack_blank_data <- rbind(blank_row_stack_mice_cols, data_stack_for_mice)
  
  # Sets correct classes for imputation - "data_stack_for_mice_colclasses" needs to be in the same order as reference dataset
  colclasses <- data_stack_for_mice_colclasses
  for (i in 1:length(colclasses)){
    if (colclasses[i] %in% c("factor", "F", "FACTOR", "Factor", "f")){
      stack_blank_data[,i] <- as.factor(stack_blank_data[,i])
    }
    if (colclasses[i] %in% c("numeric", "number", "n", "N", "NUMBER", "NUMERIC", "Number", "Numeric")){
      stack_blank_data[,i] <- as.numeric(stack_blank_data[,i])
    }
  }
  
  # Adds correct values in column "categorical_gene" according to pathogenic variant status
  stack_blank_data <- add_column_categorical_major_gene(stack_blank_data)
  stack_blank_data <- stack_blank_data[,which(!(colnames(stack_blank_data) %in% gene_list))]
  
  # Sets "variable_rfs" to the set of risk factors that are to be imputed in the individual's data
  vrf_all <- 0
  if (length(variable_rfs) == 1){
    if (variable_rfs == "all"){
      vrf_all <- 1
    }
  }
  if (vrf_all == 1){
    variable_rfs <- colnames(stack_blank_data)[!(colnames(stack_blank_data) %in% given_rfs)]
  }
  
  # Sets "mice_colnames" to the set of risk factors to use in imputation
  mice_colnames <- c(variable_rfs,given_rfs)
  
  # Log transform for Volpara percent mammographic density 
  if ("Volpara" %in% mice_colnames){
    stack_blank_data$Volpara <- log(stack_blank_data$Volpara)
  }
  
  # Log transform for alcohol, and makes 0/1 column for additional probability mass at 0
  if ("Alcohol" %in% variable_rfs | vrf_all == 1){
    mice_colnames <- mice_colnames[which(mice_colnames != "Alcohol")]
    mice_colnames <- c(mice_colnames, "Alcohol_01")
    alcohol_01 <- stack_blank_data$Alcohol
    alcohol_01[which(!is.na(alcohol_01) & alcohol_01 == 0)] <- 0
    alcohol_01[which(!is.na(alcohol_01) & alcohol_01 > 0)] <- 1
    stack_blank_data$Alcohol_01 <- as.factor(alcohol_01)
  }
  
  # Log transform for BMI
  if ("BMI" %in% mice_colnames){
    stack_blank_data$BMI <- log(stack_blank_data$BMI)
  }
  
  # Sets "quickpred_no" and "quickpred_cor" values when they are not to be used
  if (is.na(quickpred_no)){
    quickpred_no <- nrow(data_stack_for_mice)
  }
  if (is.na(quickpred_cor)){
    quickpred_cor <- 0
  }
  
  # Implements MICE
  part_mice_output <- mice(stack_blank_data[1:(quickpred_no + n_repeats),mice_colnames], pred=quickpred(stack_blank_data[1:(quickpred_no + 1),mice_colnames], mincor=quickpred_cor), m = n_chain, defaultMethod = c("norm", "logreg", "polyreg", "polr"), ignore = c(rep(TRUE, n_repeats), rep(FALSE, quickpred_no)), maxit = n_iter)
  part_mice_long <- complete(part_mice_output, action = "long")
  part_mice <- part_mice_long[which(part_mice_long$.imp == max(part_mice_long$.imp)),]
  
  # Rounding age at menarche to match up with reference
  if ("Menarche" %in% variable_rfs){
    for (i in 1:nrow(part_mice)){
      if (!is.na(part_mice$Menarche[i])){
        part_mice$Menarche[i] <- round(part_mice$Menarche[i])
      }
    }
  }
  
  # MICE repeat for nonzero alcohol values
  if ("Alcohol" %in% variable_rfs | vrf_all == 1){
    part_mice$Alcohol <- stack_blank_data$Alcohol
    part_mice_nonzero_alcohol <- subset(part_mice, !is.na(Alcohol_01) & Alcohol_01 == 1)
    n_indiv <- nrow(subset(part_mice_nonzero_alcohol, .id<=n_repeats))
    n_data <- nrow(subset(part_mice_nonzero_alcohol, .id>n_repeats))
    
    part_mice_nonzero_alcohol$Alcohol <- log(part_mice_nonzero_alcohol$Alcohol)
    
    part_mice_alcohol_output <- mice(part_mice_nonzero_alcohol[,!(colnames(part_mice_nonzero_alcohol) %in% c(".imp", ".id"))], m = 1, defaultMethod = c("norm", "logreg", "polyreg", "polr"), ignore = c(rep(TRUE, n_indiv), rep(FALSE, n_data)), maxit = 1)
    
    part_mice_alcohol_long <- complete(part_mice_alcohol_output, action = "long")
    part_mice_alcohol <- part_mice_alcohol_long[which(part_mice_alcohol_long$.imp == max(part_mice_alcohol_long$.imp)),]
    
    part_mice_alcohol$Alcohol <- exp(part_mice_alcohol$Alcohol)
    
    part_mice[which(part_mice$Alcohol_01 == 1),]$Alcohol <- part_mice_alcohol$Alcohol
    part_mice[which(part_mice$Alcohol_01 == 0),]$Alcohol <- rep(0, length(which(part_mice$Alcohol_01 == 0)))
  }
  
  # MICE repeat for age at first live birth, for those with (potentially imputed) parity > 0
  if ("First_live_birth" %in% variable_rfs | vrf_all == 1){
    flb_data <- read.csv(flb_path, colClasses = c("factor"))
    part_mice$First_live_birth <- rep(NA, nrow(part_mice))
    part_mice$First_live_birth[(n_repeats+1):nrow(part_mice)] <- as.character(flb_data$First_live_birth[1:quickpred_no])
    part_mice$First_live_birth[1:n_repeats] <- blank_stack_flb
    if ("First_live_birth" %in% given_rfs){
      part_mice$First_live_birth[1:n_repeats] <- rep(rf_values[which(given_rfs == "First_live_birth")], n_repeats)
    }
    part_mice$First_live_birth <- as.character(part_mice$First_live_birth)
    part_mice$First_live_birth <- as.factor(part_mice$First_live_birth)
    
    part_mice_nonzero_parity <- subset(part_mice, Parity != 0)
    
    n_indiv <- nrow(subset(part_mice_nonzero_parity, .id<=n_repeats))
    n_data <- nrow(subset(part_mice_nonzero_parity, .id>n_repeats))
    
    part_mice_flb_output <- mice(part_mice_nonzero_parity[,!(colnames(part_mice_nonzero_parity) %in% c(".imp", ".id", "Alcohol"))], m = 1, defaultMethod = c("norm", "logreg", "polyreg", "polr"), ignore = c(rep(TRUE, n_indiv), rep(FALSE, n_data)), maxit = 1)
    
    part_mice_flb_long <- complete(part_mice_flb_output, action = "long")
    part_mice_flb <- part_mice_flb_long[which(part_mice_flb_long$.id %in% 1:n_repeats & part_mice_flb_long$.imp == max(part_mice_flb_long$.imp)),]
    part_mice_flb$.id <- part_mice_nonzero_parity$.id[1:n_indiv]
    
    part_mice_all_parity <- subset(part_mice, .id %in% 1:n_repeats)
    for (i in 1:nrow(part_mice_all_parity)){
      if (part_mice_all_parity$.id[i] %in% part_mice_flb$.id){
        part_mice_all_parity$First_live_birth[i] <- part_mice_flb$First_live_birth[which(part_mice_flb$.id == part_mice_all_parity$.id[i])]
      }
    }
  } else {
    part_mice_all_parity <- part_mice
  }
  
  full_mice_out <- part_mice_all_parity
  colnames(full_mice_out)[colnames(full_mice_out) == ".id"] <- "IndivID"
  
  # Sets values of non-variable columns
  full_mice_out$FamID <- 1:n_repeats
  full_mice_out$MothID <- rep(NA,n_repeats)
  full_mice_out$FathID <- rep(NA,n_repeats)
  full_mice_out$Proband <- rep(1,n_repeats)
  full_mice_out$Sex <- rep("F",n_repeats)
  full_mice_out$Age <- indiv_age
  full_mice_out$Yob <- indiv_yob
  full_mice_out$Censoring_Age <- censoring_age
  full_mice_out$BRCA1t <- rep("S", n_repeats)
  full_mice_out$BRCA2t <- rep("S", n_repeats)
  full_mice_out$PALB2t <- rep("S", n_repeats)
  full_mice_out$ATMt <- rep("S", n_repeats)
  full_mice_out$CHEK2t <- rep("S", n_repeats)
  full_mice_out$RAD51Ct <- rep("S", n_repeats)
  full_mice_out$RAD51Dt <- rep("S", n_repeats)
  full_mice_out$BARD1t <- rep("S", n_repeats)
  full_mice_out$BC_PRS_alpha <- rep(bc_prs_alpha, n_repeats)
  full_mice_out$Ashkn <- rep(NA, n_repeats)
  full_mice_out$MZtwin <- rep(NA, n_repeats)
  full_mice_out$BrCa_1st <- rep(NA, n_repeats)
  full_mice_out$BrCa_2nd <- rep(NA, n_repeats)
  full_mice_out$OvCa <- rep(NA, n_repeats)
  full_mice_out$PanCa <- rep(NA, n_repeats)
  full_mice_out$ProCa <- rep(NA, n_repeats)
  full_mice_out$ER <- rep(NA, n_repeats)
  full_mice_out$PR <- rep(NA, n_repeats)
  full_mice_out$HER2 <- rep(NA, n_repeats)
  full_mice_out$CK14 <- rep(NA, n_repeats)
  full_mice_out$CK56 <- rep(NA, n_repeats)
  
  # Transforms Volpara percent mammographic density back from the log value
  if ("Volpara" %in% mice_colnames){
    full_mice_out$Volpara <- exp(full_mice_out$Volpara)
  }
  
  # Transforms BMI back from the log value
  if ("BMI" %in% mice_colnames){
    full_mice_out$BMI <- exp(full_mice_out$BMI)
  }
  
  # Replaces columns that were not used in imputation
  unused_cols <- colnames(stack_blank_data)[!(colnames(stack_blank_data) %in% given_rfs) & !(colnames(stack_blank_data) %in% variable_rfs)]
  for (i in unused_cols){
    full_mice_out[,i] <- rep(NA, nrow(full_mice_out))
  }
  
  full_mice_out_boadicea <- full_mice_out
  
  # Reformatting parity column for BOADICEA input
  full_mice_out_boadicea$Parity <- as.character(full_mice_out_boadicea$Parity)
  for (i in 1:n_repeats){
    if (!is.na(full_mice_out_boadicea$Parity[i]) & full_mice_out_boadicea$Parity[i] == ">2"){
      full_mice_out_boadicea$Parity[i] <- 3
    }
  }
  
  # Reformatting OC use column for BOADICEA input
  full_mice_out_boadicea$OC_Use <- as.character(full_mice_out_boadicea$OC_Use)
  for (i in 1:n_repeats){
    if (full_mice_out_boadicea$OC_Use[i] == "N"){
      full_mice_out_boadicea$OC_Use[i] <- NA
    }
  }
  
  # Reformatting menopausal HRT use column for BOADICEA input
  full_mice_out_boadicea$MHT_Use <- as.character(full_mice_out_boadicea$MHT_Use)
  for (i in 1:n_repeats){
    if (full_mice_out_boadicea$MHT_Use[i] == "N"){
      full_mice_out_boadicea$MHT_Use[i] <- NA
    }
  }
  
  # Reformatting age / status of menopause column for BOADICEA input
  full_mice_out_boadicea$Menopause <- as.character(full_mice_out_boadicea$Menopause)
  for (i in 1:n_repeats){
    if (full_mice_out_boadicea$Menopause[i] == "N"){
      full_mice_out_boadicea$Menopause[i] <- NA
    } else if (full_mice_out_boadicea$Menopause[i] == "<40"){
      full_mice_out_boadicea$Menopause[i] <- 37
    } else if (full_mice_out_boadicea$Menopause[i] == "40-44"){
      full_mice_out_boadicea$Menopause[i] <- 42
    } else if (full_mice_out_boadicea$Menopause[i] == "45-49"){
      full_mice_out_boadicea$Menopause[i] <- 47
    } else if (full_mice_out_boadicea$Menopause[i] == "50-54"){
      full_mice_out_boadicea$Menopause[i] <- 52
    } else if (full_mice_out_boadicea$Menopause[i] == ">54"){
      full_mice_out_boadicea$Menopause[i] <- 57
    }
  }
  
  # Reformatting age at first live birth column for BOADICEA input
  full_mice_out_boadicea$First_live_birth <- as.character(full_mice_out_boadicea$First_live_birth)
  for (i in 1:n_repeats){
    if (!is.na(full_mice_out_boadicea$First_live_birth[i])){
      if (full_mice_out_boadicea$First_live_birth[i] == "<20"){
        full_mice_out_boadicea$First_live_birth[i] <- 18
      } else if (full_mice_out_boadicea$First_live_birth[i] == "20-24"){
        full_mice_out_boadicea$First_live_birth[i] <- 23
      } else if (full_mice_out_boadicea$First_live_birth[i] == "25-29"){
        full_mice_out_boadicea$First_live_birth[i] <- 28
      } else if (full_mice_out_boadicea$First_live_birth[i] == ">29"){
        full_mice_out_boadicea$First_live_birth[i] <- 33
      }
    }
  }
  
  # Reformatting major gene columns for BOADICEA input, if was imputed as one column
  if (gene_imp_type == "categorical"){
    full_mice_out_boadicea <- categorical_to_binary_major_gene(full_mice_out_boadicea)
  }
  
  # Reformatting data columns for BOADICEA input
  rownames(full_mice_out_boadicea) <- 1:n_repeats
  full_mice_out_boadicea <- full_mice_out_boadicea[,c("FamID", "IndivID", "FathID", 
                                                      "MothID", "Proband", "Sex", 
                                                      "MZtwin", "Age", "Yob", 
                                                      "Ashkn", "BrCa_1st", "BrCa_2nd", 
                                                      "OvCa", "ProCa", "PanCa", 
                                                      "Censoring_Age", "BRCA1t", "BRCA1r", 
                                                      "BRCA2t", "BRCA2r", "PALB2t", 
                                                      "PALB2r", "ATMt", "ATMr", 
                                                      "CHEK2t", "CHEK2r", "RAD51Ct", 
                                                      "RAD51Cr", "RAD51Dt", "RAD51Dr", 
                                                      "BARD1t", "BARD1r", "ER", 
                                                      "PR", "HER2", "CK14", 
                                                      "CK56", "Menarche", "Parity", 
                                                      "First_live_birth", "OC_Use", "MHT_Use", 
                                                      "BMI", "Alcohol", "Menopause", 
                                                      "Height", "BC_PRS_alpha", 
                                                      "BC_PRS_z", "mother.breast_v0", "sibling.breast_v0", 
                                                      "sibling.prostate_v0", "father.prostate_v0", "Volpara")]
  
  full_mice_out_boadicea$FathID <- rep("1.f", n_repeats)
  full_mice_out_boadicea$MothID <- rep("1.m", n_repeats)
  
  # Adding family history from pedigree structures in reference data 
  family_reference_dataset <- read.csv(fams_path, colClasses = rep("character", ncol(read.csv(fams_path))))
  family_reference_dataset <- family_reference_dataset[,colnames(family_reference_dataset) != "X"]
  
  family_numbers <- rep(0,n_repeats)
  issues <- c()
  for (i in 1:n_repeats){
    choice_subset <- subset(family_reference_dataset, Proband == 1 & mother.breast_v0 == full_mice_out_boadicea$mother.breast_v0[i] & sibling.breast_v0 == full_mice_out_boadicea$sibling.breast_v0[i] & sibling.prostate_v0 == full_mice_out_boadicea$sibling.prostate_v0[i] & father.prostate_v0 == full_mice_out_boadicea$father.prostate_v0[i])
    if (nrow(choice_subset) == 0){
      issues <- c(issues, i)
      choice_subset <- subset(family_reference_dataset, Proband == 1 & mother.breast_v0 == full_mice_out_boadicea$mother.breast_v0[i] & sibling.breast_v0 == full_mice_out_boadicea$sibling.breast_v0[i] & father.prostate_v0 == full_mice_out_boadicea$father.prostate_v0[i])
    }
    family_numbers[i] <- sample(x = choice_subset$FamID, size = 1)
  }
  
  if(length(issues) != 0){
    warning(paste0("Issue family: ", issues))
  }
  
  for(i in 1:n_repeats){
    family_to_add <- subset(family_reference_dataset, FamID == family_numbers[i] & (is.na(Proband) | Proband != 1))
    family_to_add$FamID <- rep(full_mice_out_boadicea$FamID[i], nrow(family_to_add))
    full_mice_out_boadicea[(nrow(full_mice_out_boadicea)+1):(nrow(full_mice_out_boadicea)+nrow(family_to_add)),] <- family_to_add[,c("FamID", "IndivID", "FathID", 
                                                                                                                                     "MothID", "Proband", "Sex", 
                                                                                                                                     "MZtwin", "Age", "Yob", 
                                                                                                                                     "Ashkn", "BrCa_1st", "BrCa_2nd", 
                                                                                                                                     "OvCa", "ProCa", "PanCa", 
                                                                                                                                     "Censoring_Age", "BRCA1t", "BRCA1r", 
                                                                                                                                     "BRCA2t", "BRCA2r", "PALB2t", 
                                                                                                                                     "PALB2r", "ATMt", "ATMr", 
                                                                                                                                     "CHEK2t", "CHEK2r", "RAD51Ct", 
                                                                                                                                     "RAD51Cr", "RAD51Dt", "RAD51Dr", 
                                                                                                                                     "BARD1t", "BARD1r", "ER", 
                                                                                                                                     "PR", "HER2", "CK14", 
                                                                                                                                     "CK56", "Menarche", "Parity", 
                                                                                                                                     "First_live_birth", "OC_Use", "MHT_Use", 
                                                                                                                                     "BMI", "Alcohol", "Menopause", 
                                                                                                                                     "Height", "BC_PRS_alpha", 
                                                                                                                                     "BC_PRS_z", "mother.breast_v0", "sibling.breast_v0", 
                                                                                                                                     "sibling.prostate_v0", "father.prostate_v0", "Volpara")]
  }
  
  full_mice_out_boadicea <- full_mice_out_boadicea[order(full_mice_out_boadicea$FamID),]
  
  # Reformatting Volpara percent mammographic density for BOADICEA input
  mamm_density <- as.numeric(full_mice_out_boadicea$Volpara)
  mamm_density <- (mamm_density / 100) + rep(20, length(mamm_density))
  menopause_status <- rep(NA, length(mamm_density))
  menopause_status[which(is.na(full_mice_out_boadicea$Menopause) & full_mice_out_boadicea$Proband == 1)] <- rep(1, length(which(is.na(full_mice_out_boadicea$Menopause) & full_mice_out_boadicea$Proband == 1)))
  menopause_status[which(!is.na(full_mice_out_boadicea$Menopause))] <- rep(2, length(which(!is.na(full_mice_out_boadicea$Menopause))))
  mamm_density <- mamm_density + menopause_status
  
  full_mice_out_boadicea$Mamm_density <- mamm_density
  
  # Save imputed data
  save_dataset(full_mice_out_boadicea, filepath = save_filepath)
  
}
