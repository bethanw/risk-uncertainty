confidence_interval <- function(results, round = "dp", dp = "auto", format = "decimal", alpha = 0.05){
  sort_results <- sort(results)
  n_results <- length(results)
  if (format == "percentage"){
    if (dp == "auto"){
      dp <- 1
    }
    
    lower <- round(sort_results[round((alpha/2)*n_results)]*100, digits = dp)
    upper <- round(sort_results[round((1-(alpha/2))*n_results)]*100, digits = dp)
    
    if(decimalplaces(lower) != dp){
      lower_0s <- paste0(c(rep(0, dp - decimalplaces(lower))))
    } else {
      lower_0s <- ""
    }
    if(decimalplaces(upper) != dp){
      upper_0s <- paste0(c(rep(0, dp - decimalplaces(upper))))
    } else {
      upper_0s <- ""
    }
    
    if (decimalplaces(lower) == 0 & dp != 0){
      lower_point <- "."
    } else {
      lower_point <- ""
    }
    
    if (decimalplaces(upper) == 0 & dp != 0){
      upper_point <- "."
    } else {
      upper_point <- ""
    }
    
    return(noquote(paste0("(", lower, lower_point, lower_0s, "%, ", upper, upper_point, upper_0s, "%)")))
  } 
  
  if (format == "decimal") {
    if(dp == "auto"){
      dp <- 3
    }
    
    lower <- round(sort_results[round(0.025*n_results)], digits = dp)
    upper <- round(sort_results[round(0.975*n_results)], digits = dp)
    
    if(decimalplaces(lower) != dp){
      lower_0s <- paste0(c(rep(0, dp - decimalplaces(lower))))
    } else {
      lower_0s <- ""
    }
    if(decimalplaces(upper) != dp){
      upper_0s <- paste0(c(rep(0, dp - decimalplaces(upper))))
    } else {
      upper_0s <- ""
    }
    
    if (decimalplaces(lower) == 0 & dp != 0){
      lower_point <- "."
    } else {
      lower_point <- ""
    }
    
    if (decimalplaces(upper) == 0 & dp != 0){
      upper_point <- "."
    } else {
      upper_point <- ""
    }
    
    return(noquote(paste0("(", lower, lower_point, lower_0s, ", ", upper, upper_point, upper_0s, ")")))
  }
}
