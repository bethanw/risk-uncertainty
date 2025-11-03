make_sankey_from_results_values <- function(results1, results2, results3 = "blank", original_category = NA, inc_low_risk = FALSE, low_risk_val = 0.005, add_flow_labels = FALSE, node_labels = c(), cramp_cutoff = 35, cramp_shift = 25, cramp_x_shift = 0.1, cramp_tb_y_shift = 10){
  require(ggsankey)
  
  results3_present <- TRUE
  if (length(results3) == 1) {
    if (results3 == "blank") {
      results3_present <- FALSE
    }
  }
  
  pop_risk_1 <- which(results1 < 0.03)
  mod_risk_1 <- which(results1 >= 0.03 & results1 < 0.08)
  hig_risk_1 <- which(results1 >= 0.08)
  categories_1 <- rep(NA, length(results1))
  categories_1[pop_risk_1] <- "Near-population Risk"
  categories_1[mod_risk_1] <- "Moderate Risk"
  categories_1[hig_risk_1] <- "High Risk"
  
  pop_risk_2 <- which(results2 < 0.03)
  mod_risk_2 <- which(results2 >= 0.03 & results2 < 0.08)
  hig_risk_2 <- which(results2 >= 0.08)
  categories_2 <- rep(NA, length(results2))
  categories_2[pop_risk_2] <- "Near-population Risk"
  categories_2[mod_risk_2] <- "Moderate Risk"
  categories_2[hig_risk_2] <- "High Risk"
  
  if (inc_low_risk == TRUE){
    low_risk_1 <- which(results1 < low_risk_val)
    low_risk_2 <- which(results2 < low_risk_val)
    categories_1[low_risk_1] <- "Low Risk"
    categories_2[low_risk_2] <- "Low Risk"
  }
  
  if (results3_present){
    pop_risk_3 <- which(results3 < 0.03)
    mod_risk_3 <- which(results3 >= 0.03 & results3 < 0.08)
    hig_risk_3 <- which(results3 >= 0.08)
    categories_3 <- rep(NA, length(results3))
    categories_3[pop_risk_3] <- "Near-population Risk"
    categories_3[mod_risk_3] <- "Moderate Risk"
    categories_3[hig_risk_3] <- "High Risk"
    
    if (inc_low_risk == TRUE){
      low_risk_3 <- which(results3 < low_risk_val)
      categories_3[low_risk_3] <- "Low Risk"
    }
  }
  
  if (is.na(original_category)){
    if (mean(results1) < 0.03){
      original_category <- "Near-population Risk"
    } else if (mean(results1) < 0.08){
      original_category <- "Moderate Risk"
    } else {
      original_category <- "High Risk"
    }
  }
  
  Original_Categories <- rep(original_category, length(results1))
  
  if(!results3_present){
    df <- data.frame(Original_Categories, categories_1, categories_2) %>%
      make_long("Original_Categories", "categories_1", "categories_2")
  } else {
    df <- data.frame(Original_Categories, categories_1, categories_2, categories_3) %>%
      make_long("Original_Categories", "categories_1", "categories_2", "categories_3")
  }
  
  
  factor_risk_levels <- c("Near-population Risk", "Moderate Risk", "High Risk")
  if (inc_low_risk == TRUE){
    factor_risk_levels <- c("Low Risk", factor_risk_levels)
  }
  
  df$node <- factor(df$node, levels = factor_risk_levels)
  df$next_node <- factor(df$next_node, levels = factor_risk_levels)
  
  colours_fill <- c("High Risk" = "red", "Moderate Risk" = "orange", "Near-population Risk" = "royalblue")
  if (inc_low_risk == TRUE){
    colours_fill <- c(colours_fill, "Low Risk" = "yellowgreen")
  }
  
  sankey_out <- ggplot(df, aes(x = x, 
                               next_x = next_x, 
                               node = node, 
                               next_node = next_node,
                               fill = factor(node))) +
    geom_sankey(flow.alpha = 0.5, node.color = 1, flow.color = "#101010") +
    scale_fill_manual(values = colours_fill) +
    theme_sankey(base_size = 16) +
    guides(fill = guide_legend(title = "Risk Category")) +
    labs(x = "") + 
    scale_x_discrete(labels=node_labels)
    
  
  
  
  # get flow size for each flow
  flow_labels <- df %>% group_by(x, node, next_x, next_node) %>% tally() %>% drop_na()
  
  # get corresponding positions of flows from the sankey plot
  flow_info <- layer_data(sankey_out) %>% select(xmax, flow_start_ymax, flow_start_ymin) %>% distinct() # get flow related key positions related from the plot
  flow_info <- flow_info[with(flow_info, order(xmax, flow_start_ymax)), ] # order the positions to match the order of flow_labels
  rownames(flow_info) <- NULL # drop the row indexes
  flow_info <- cbind(as.data.frame(flow_info), as.data.frame(flow_labels)) # bind the flow positions and the corresponding labels
  
  print(flow_info[,c("x","node", "next_x", "next_node", "n")])
  
  flow_info$y_lab_coord <- (flow_info$flow_start_ymin + flow_info$flow_start_ymax)/2
  flow_info$cramp <- rep(0, nrow(flow_info))
  
  flow_info$n_flows_from_node <- rep(NA, nrow(flow_info))
  for (i in 1:nrow(flow_info)){
    flow_info$n_flows_from_node[i] <- length(which(flow_info$x == flow_info$x[i] & flow_info$node == flow_info$node[i]))
  }
  
  flow_info$cramp_direction <- rep(NA, nrow(flow_info))
  for (i in 1:nrow(flow_info)){
    if (flow_info$n_flows_from_node[i] == 1){
      flow_info$cramp_direction[i] <- 0
    }
    
    if (flow_info$n_flows_from_node[i] == 2){
      if (flow_info$flow_start_ymax[i] > flow_info$flow_start_ymax[which(flow_info$x == flow_info$x[i] & flow_info$node == flow_info$node[i] & flow_info$flow_start_ymax != flow_info$flow_start_ymax[i])]){
        flow_info$cramp_direction[i] <- 1
      } else {
        flow_info$cramp_direction[i] <- -1
      }
    }
    
    if (flow_info$n_flows_from_node[i] == 3){
      flow_subset <- subset(flow_info, x == flow_info$x[i] & node == flow_info$node[i])
      if (flow_subset$flow_start_ymax[1] == flow_info$flow_start_ymax[i]){
        flow_info$cramp_direction[i] <- -1
      } else if (flow_subset$flow_start_ymax[2] == flow_info$flow_start_ymax[i]) {
        flow_info$cramp_direction[i] <- 0
      } else {
        flow_info$cramp_direction[i] <- 1
      }
    }
    
  }
  
  for (i in 1:nrow(flow_info)){
    if (flow_info$flow_start_ymax[i] - flow_info$flow_start_ymin[i] < cramp_cutoff){
      flow_info$cramp[i] <- 1
      if (flow_info$cramp_direction[i] == 1){
        flow_info$y_lab_coord[i] <- flow_info$flow_start_ymax[i] + (cramp_shift)
      } else if (flow_info$cramp_direction[i] == -1){
        flow_info$y_lab_coord[i] <- flow_info$flow_start_ymin[i] - (cramp_shift)
      }
      
      if (flow_info$n_flows_from_node[i] == 3 & flow_info$cramp_direction[i] == 0){
        flow_info$xmax[i] <- flow_info$xmax[i] + cramp_x_shift
        if (flow_info$flow_start_ymax[i+1] == max(flow_info$flow_start_ymax)){
          flow_info$y_lab_coord[i] <- flow_info$flow_start_ymin[i] - cramp_tb_y_shift
        }
        if (flow_info$flow_start_ymin[i-1] == min(flow_info$flow_start_ymin)){
          flow_info$y_lab_coord[i] <- flow_info$flow_start_ymin[i] + cramp_tb_y_shift
        }
      }
      
      if (flow_info$n_flows_from_node[i] == 1 & flow_info$cramp_direction[i] == 0){
        flow_info$y_lab_coord[i] <- flow_info$flow_start_ymax[i] + cramp_tb_y_shift*5
      }
      
    }
  }
  
  if (add_flow_labels){
    # add labels to the flows
    for (i in 1:nrow(flow_info)){
      sankey_out <- sankey_out + annotate("text", x = flow_info$xmax[i],
                                          y = flow_info$y_lab_coord[i],
                                          label = sprintf("%d", flow_info$n[i]), hjust = -1,
                                          size = 2) + ylab("")
    }
  }
  return(sankey_out)
}
