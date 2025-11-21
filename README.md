Code to produce figures and impute missing covariates for work: "Modelling Individual-level Uncertainty due to Missing Data in Personalised Breast Cancer Risk Prediction"

Uses BOADICEA-formatted UK Biobank and KARMA data

R version: 2024.09.1+394 (2024.09.1+394)

Libraries needed: 
- "mice"
- "ggplot2" (to make graphs; Figures 2, 4, S1, S2)
- "ggsankey" (to make sankey plots; Figure 3)

Functions: 

- get_covariate_distribution.R: main function that outputs imputed BOADICEA-format data for an individual with present and missing information

- add_column_categorical_major_gene.R: changes format of major gene columns in dataframe from individual columns of binary variables for each gene to single categorical variable
- categorical_to_binary_major_gene.R: changes format of major gene columns in dataframe from single categorical variable to individual columns of binary variables for each gene
- confidence_interval.R: calculates a confidence interval for the risk uncertainty distribution
- decimalplaces.R: returns the number of decimal places a given number has
- make_graph_out.R: makes a graph of the uncertainty distribution from the BOADICEA output from the imputed individuals
- make_sankey_from_results_values.R: produces a sankey diagram from given BOADICEA outputs
- reformat_output_dataset.R: puts data in the correct format for BOADICEA, from the imputation-formatted data
- save_dataset.R: saves data for use in BOADICEA
