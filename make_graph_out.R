make_graph_out <- function(out, bw, x_low = 0, x_high, include_legend = TRUE, add_smooth = FALSE, point_estimate = NA){
  out_df <- data.frame(out)
  out_df$colours <- ifelse(out_df$out < 0.03, "Near-population risk", ifelse(out_df$out < 0.08, "Moderate risk", "High risk"))
  
  plot <- ggplot(out_df, aes(x = out)) + geom_histogram(data = out_df, binwidth = bw, boundary = bw, aes(fill = colours)) + 
    #labs(x = "10-year risk from 40", y = "Count", fill = "") + scale_fill_manual(values = c("Near-population risk" = "lightskyblue", "Moderate risk" = "orange", "High risk" = "red")) + 
    labs(x = "10-year risk from 40", y = "Count", fill = "") + scale_fill_manual(values = c("Near-population risk" = "royalblue", "Moderate risk" = "orange", "High risk" = "red")) + 
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', colour = "white"), 
          panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid', colour = "white"),
          axis.line = element_line(color = "black", size = 0.2)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    # scale_x_continuous(expand = c(0, 0), limits = c(x_low,x_high))
    scale_x_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(x_low,x_high))
  
  if (include_legend == FALSE){
    plot <- plot + theme(legend.position="none")
  }
  
  if (add_smooth == TRUE){
    nbin <- (x_high - x_low)/bw
    bin_heights <- rep(0, nbin)
    for (i in 1:nbin){
      bin_heights[i] <- length(which(out >= x_low+((i-1)*bw) & out < x_low+(i*bw)))
    }
    plot <- plot + 
    stat_spline(data = data.frame(x_coord = seq(from = x_low+(bw/2), to = x_high-(bw/2), length.out = nbin), y_coord = bin_heights), mapping = aes(x = x_coord, y = y_coord), linetype="dashed", color = "black")
  }
  
  if (!is.na(point_estimate)){
    plot <- plot + geom_vline(xintercept = point_estimate, linetype = "solid", color = "black", size = 0.5)
  }
  
  return(plot)
}
