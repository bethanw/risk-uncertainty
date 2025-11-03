get_covariate_distribution <- function(data_stack_for_mice_path, 
                                                fh_in_br = FALSE, 
                                                gene_list = c("BRCA1r", "BRCA2r", "PALB2r", "CHEK2r", "ATMr", "RAD51Cr", "RAD51Dr", "BARD1r"), 
                                                n_chain = 1, 
                                                make_blank_stack = TRUE, 
                                                data_stack_for_mice_colclasses, 
                                                variable_rfs = "all", 
                                                blank_row_path, 
                                                blank_row_colclasses, 
                                                flb_path, 
                                                fams_path, 
                                                given_rfs = c(), 
                                                rf_values = c(), 
                                                quickpred_no = 50000, 
                                                quickpred_cor = 1, 
                                                inc_family = TRUE, 
                                                n_repeats = 1000, 
                                                n_iter = 10, 
                                                indiv_yob, 
                                                indiv_age, 
                                                censoring_age, 
                                                bc_prs_alpha = 0.44, 
                                                save_filepath = "default",
                                                gene_imp_type = "categorical", 
                                                output_diagnostics = FALSE, 
                                                plot_density = FALSE) {
  
  # certain combinations of missing risk factors currently not working
  
  require(mice)
  
  if ("pedigree" %in% given_rfs){
    impute_pedigree <- 0
    fh_given <- read.csv(rf_values[which(given_rfs == "pedigree")])
  } else {
    impute_pedigree <- 1
  }
  
  
  # get blank row
  blank_row <- read.csv(blank_row_path, colClasses = rep("character", 51))
  
  if (impute_pedigree == 1 & fh_in_br == FALSE){
    blank_row$mother.breast_v0 <- NA
    blank_row$sibling.breast_v0 <- NA
    blank_row$sibling.prostate_v0 <- NA
    blank_row$father.prostate_v0 <- NA
  } #else {
  #   blank_row$mother.breast_v0 <- fh_given$mother.breast_v0
  #   blank_row$sibling.breast_v0 <- fh_given$sibling.breast_v0
  #   blank_row$sibling.prostate_v0 <- fh_given$sibling.prostate_v0
  #   blank_row$father.prostate_v0 <- fh_given$father.prostate_v0
  # }
  
  
  # add known information
  
  for (i in 1:length(given_rfs)){
    blank_row[,given_rfs[i]] <- rf_values[i]
  }
  
  if (make_blank_stack == TRUE){
    blank_row_stack <- blank_row
    for (i in 2:n_repeats){
      blank_row_stack[i,] <- blank_row
    }
  } else {
    blank_row_stack <- blank_row
  }
  
  blank_row_stack$IndivID <- 1:n_repeats
  blank_row_stack$FamID <- 1:n_repeats
  
  

  # introduce the data
  # data here is ukbiobank mixed with karma - relevant cols
  
  data_stack_for_mice <- read.csv(data_stack_for_mice_path)
  data_stack_for_mice <- data_stack_for_mice[,colnames(data_stack_for_mice) != "X"]
  
  if (!("First_live_birth" %in% colnames(blank_row_stack))){
    blank_row_stack$First_live_birth <- rep(NA, nrow(blank_row_stack))
  }
  blank_stack_flb <- blank_row_stack$First_live_birth
  
  if (!("categorical_gene" %in% colnames(blank_row_stack))){
    blank_row_stack$categorical_gene <- rep(NA, nrow(blank_row_stack))
  }
  blank_row_stack_mice_cols <- blank_row_stack[,colnames(data_stack_for_mice)]
  stack_blank_data <- rbind(blank_row_stack_mice_cols, data_stack_for_mice)
  
  colclasses <- data_stack_for_mice_colclasses
  for (i in 1:length(colclasses)){
    if (colclasses[i] %in% c("factor", "F", "FACTOR", "Factor", "f")){
      stack_blank_data[,i] <- as.factor(stack_blank_data[,i])
    }
    if (colclasses[i] %in% c("numeric", "number", "n", "N", "NUMBER", "NUMERIC", "Number", "Numeric")){
      stack_blank_data[,i] <- as.numeric(stack_blank_data[,i])
    }
  }
  
  
  
  if (gene_imp_type == "categorical"){
    stack_blank_data <- add_column_categorical_major_gene(stack_blank_data)
    stack_blank_data <- stack_blank_data[,which(!(colnames(stack_blank_data) %in% c("BRCA1r", "BRCA2r", "PALB2r", "ATMr", "CHEK2r", "RAD51Cr", "RAD51Dr", "BARD1r")))]
  } else {
    stack_blank_data <- stack_blank_data[,which(colnames(stack_blank_data) != "categorical_gene")]
  }
  
  
  
  vrf_all <- 0
  if (length(variable_rfs) == 1){
    if (variable_rfs == "all"){
      vrf_all <- 1
    }
  }
  
  if (vrf_all == 1){
    variable_rfs <- colnames(stack_blank_data)[!(colnames(stack_blank_data) %in% given_rfs)]
  } else if (all(!(variable_rfs %in% colnames(stack_blank_data)))){
    stop("variable_rfs not all valid rfs for imputation")
  }
  
  # can add in 'post' (?) - positive values only for some cols
  mice_colnames <- c(variable_rfs,given_rfs)
  if (gene_imp_type == "categorical"){
    mice_colnames <- mice_colnames[which(!(mice_colnames %in% gene_list))]
  } else {
    mice_colnames <- mice_colnames[which(mice_colnames != "categorical_gene")]
  }
  
  if ("Volpara" %in% mice_colnames){
    stack_blank_data$Volpara <- log(stack_blank_data$Volpara)
  }
  
  if ("Alcohol" %in% variable_rfs | vrf_all == 1){
    mice_colnames <- mice_colnames[which(mice_colnames != "Alcohol")]
    mice_colnames <- c(mice_colnames, "Alcohol_01")
    
    alcohol_01 <- stack_blank_data$Alcohol
    alcohol_01[which(!is.na(alcohol_01) & alcohol_01 == 0)] <- 0
    alcohol_01[which(!is.na(alcohol_01) & alcohol_01 > 0)] <- 1
    stack_blank_data$Alcohol_01 <- as.factor(alcohol_01)
  }
  
  if ("BMI" %in% mice_colnames){
    stack_blank_data$BMI <- log(stack_blank_data$BMI)
  }
  
  # do the mice
  
  #part_mice_output <- mice(stack_blank_data[1:(quickpred_no + n_repeats),mice_colnames], pred=quickpred(stack_blank_data[1:(quickpred_no + 1),], mincor=quickpred_cor), m = 1, defaultMethod = c("norm", "logreg", "polyreg", "polr"), ignore = c(rep(TRUE, n_repeats), rep(FALSE, quickpred_no)), maxit = n_iter)
  part_mice_output <- mice(stack_blank_data[1:(quickpred_no + n_repeats),mice_colnames], pred=quickpred(stack_blank_data[1:(quickpred_no + 1),mice_colnames], mincor=quickpred_cor), m = n_chain, defaultMethod = c("norm", "logreg", "polyreg", "polr"), ignore = c(rep(TRUE, n_repeats), rep(FALSE, quickpred_no)), maxit = n_iter)
  #part_mice_output <- mice(stack_blank_data[1:(quickpred_no + n_repeats),mice_colnames], m = 1, defaultMethod = c("norm", "logreg", "polyreg", "polr"), ignore = c(rep(TRUE, n_repeats), rep(FALSE, quickpred_no)), maxit = n_iter)
  part_mice_long <- complete(part_mice_output, action = "long")
  part_mice <- part_mice_long[which(part_mice_long$.imp == max(part_mice_long$.imp)),]
  
  if (output_diagnostics == TRUE){
    print(plot(part_mice_output))
  }
  
  if (plot_density == TRUE){
    #print(densityplot(part_mice_output))
    print(densityplot(part_mice_output, ~Volpara))
    #print(densityplot(part_mice_output, ~Menopause))
    print(densityplot(part_mice_output, ~BMI))
    #print(densityplot(part_mice_output, ~Parity))
    print(densityplot(part_mice_output, ~OC_Use))
    print(densityplot(part_mice_output, ~Menarche))
    print(densityplot(part_mice_output, ~BC_PRS_z))
    #print(densityplot(part_mice_output, ~mother.breast_v0))
  }
  
  plot(convergence(part_mice_output, diagnostic = "gr")[1:10, 3])
  
  
  if ("Menarche" %in% variable_rfs){
    for (i in 1:nrow(part_mice)){
      if (!is.na(part_mice$Menarche[i])){
        part_mice$Menarche[i] <- round(part_mice$Menarche[i])
      }
    }
  }
  
  
  
  # add back alcohol values
  
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
  
  
  
  # add back flb
  
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
    
    
    
    # put dataset back together
    
    part_mice_all_parity <- subset(part_mice, .id %in% 1:n_repeats)
    for (i in 1:nrow(part_mice_all_parity)){
      if (part_mice_all_parity$.id[i] %in% part_mice_flb$.id){
        part_mice_all_parity$First_live_birth[i] <- part_mice_flb$First_live_birth[which(part_mice_flb$.id == part_mice_all_parity$.id[i])]
      }
    }
  } else {
    part_mice_all_parity <- part_mice
  }
  
  
  
  
  
  
  # add back deterministic cols
  
  full_mice_out <- part_mice_all_parity
  colnames(full_mice_out)[colnames(full_mice_out) == ".id"] <- "IndivID"
  
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
  
  if ("Volpara" %in% mice_colnames){
    full_mice_out$Volpara <- exp(full_mice_out$Volpara)
  }
  
  if ("BMI" %in% mice_colnames){
    full_mice_out$BMI <- exp(full_mice_out$BMI)
  }
  
  # add back unused cols
  
  unused_cols <- colnames(stack_blank_data)[!(colnames(stack_blank_data) %in% given_rfs) & !(colnames(stack_blank_data) %in% variable_rfs)]
  
  for (i in unused_cols){
    full_mice_out[,i] <- rep(NA, nrow(full_mice_out))
  }
  
  # formatting for boadicea input
  
  full_mice_out_boadicea <- full_mice_out
  
  full_mice_out_boadicea$Parity <- as.character(full_mice_out_boadicea$Parity)
  for (i in 1:n_repeats){
    if (!is.na(full_mice_out_boadicea$Parity[i]) & full_mice_out_boadicea$Parity[i] == ">2"){
      full_mice_out_boadicea$Parity[i] <- 3
    }
  }       # change once can incorporate more complex fh - parity corresponding to the number of children in the family tree
  
  full_mice_out_boadicea$OC_Use <- as.character(full_mice_out_boadicea$OC_Use)
  for (i in 1:n_repeats){
    if (full_mice_out_boadicea$OC_Use[i] == "N"){
      full_mice_out_boadicea$OC_Use[i] <- NA
    }
  }
  
  full_mice_out_boadicea$MHT_Use <- as.character(full_mice_out_boadicea$MHT_Use)
  for (i in 1:n_repeats){
    if (full_mice_out_boadicea$MHT_Use[i] == "N"){
      full_mice_out_boadicea$MHT_Use[i] <- NA
    }
  }
  
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
  

  if (gene_imp_type == "categorical"){
    full_mice_out_boadicea <- categorical_to_binary_major_gene(full_mice_out_boadicea)
  }
  
  
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
  
  
  
  # add fh
  
  full_mice_out_boadicea$FathID <- rep("1.f", n_repeats)
  full_mice_out_boadicea$MothID <- rep("1.m", n_repeats)
  
  if (impute_pedigree == 1){
    
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
    
  } else {              # add fh that is already specified - in file, location is in rf_values
    for (i in 1:n_repeats){
      fh_given$FamID <- rep(i, nrow(fh_given))
      full_mice_out_boadicea[(nrow(full_mice_out_boadicea)+1):(nrow(full_mice_out_boadicea)+nrow(fh_given)),] <- family_to_add[,c("FamID", "IndivID", "FathID", 
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
  }
  
  # # formatting mammographic density for volpara
  mamm_density <- as.numeric(full_mice_out_boadicea$Volpara)
  mamm_density <- (mamm_density / 100) + rep(20, length(mamm_density))
  menopause_status <- rep(NA, length(mamm_density))
  menopause_status[which(is.na(full_mice_out_boadicea$Menopause) & full_mice_out_boadicea$Proband == 1)] <- rep(1, length(which(is.na(full_mice_out_boadicea$Menopause) & full_mice_out_boadicea$Proband == 1)))
  menopause_status[which(!is.na(full_mice_out_boadicea$Menopause))] <- rep(2, length(which(!is.na(full_mice_out_boadicea$Menopause))))
  mamm_density <- mamm_density + menopause_status

  full_mice_out_boadicea$Mamm_density <- mamm_density

  save_dataset(full_mice_out_boadicea, filepath = save_filepath)

}



