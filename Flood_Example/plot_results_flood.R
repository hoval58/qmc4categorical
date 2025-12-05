
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

#plots the variance for of the 5 estimators considered in the paper: MC, plain RQMC, RQMC with powers of 2, per stratum RQMC, importance adjusted RQMC estimators.

plot_variance <- function(var_lists, sample_sizes, file = NULL) {
  
  # Combine all variance lists into a single data frame
  df_combined <- do.call(rbind, lapply(names(var_lists), function(name) {
    data.frame(
      SampleSize = sample_sizes,
      Variance = unlist(var_lists[[name]]),
      Method = name
    )
  }))
  
  # Order methods in the legend
  desired_order <- c(
    "MC",
    "RQMC",
    "RQMC powers of 2",
    "RQMC L independent",
    "RQMC adjusted"
  )
  df_combined$Method <- factor(df_combined$Method, levels = desired_order)
  
  # Create the plot
  p <- ggplot(df_combined, aes(
    x = SampleSize,
    y = Variance,
    color = Method,
    shape = Method
  )) +
    geom_point(size = 3) +
    geom_line() +
    scale_x_continuous(trans = 'log2', breaks = 2^seq(3, 12, by = 1)) +
    scale_y_continuous(
      trans = 'log2',
      breaks = 10^seq(-9, -1, by = 2),
      labels = scales::label_scientific(digits = 3)
    ) +
    scale_shape_manual(values = c(
      "MC" = 18,
      "RQMC" = 15,
      "RQMC powers of 2" = 9,
      "RQMC L independent" = 19,
      "RQMC adjusted" = 17
    )) +
    labs(x = "Sample Size", y = "Variance", color = NULL, shape = NULL) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12),
      legend.position = c(0.02, 0.02),
      legend.justification = c("left", "bottom"),
      legend.text = element_text(size = 12),
      legend.key.size = unit(0.8, "lines"),
      legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3),
      legend.margin = margin(6, 6, 6, 6),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    guides(
      color = guide_legend(nrow = 5, byrow = TRUE),
      shape = guide_legend(nrow = 5, byrow = TRUE)
    )
  
  
  print(p)
  if (!is.null(file)) {
    ggsave(file, plot = p, width = 8 , height = 6, units = "in") #save file
  }
  
  return(p)
}

#plots the variance for the importance adjusted estimators using rho=1,2,3 and infinity(equal allocation)

plot_rho_variance <- function(
    var_lists,
    sample_sizes,
    slope = -1.8,
    file = NULL
) {
  
  
  # ---- Combine simulation variances into one dataframe ----
  df_combined <- do.call(rbind, lapply(names(var_lists), function(name) {
    data.frame(
      SampleSize = sample_sizes,
      Variance = unlist(var_lists[[name]]),
      Method = name
    )
  }))
  
  # ---- adds a reference line ----
  
  first_point <- df_combined[df_combined$Method == "rho=2", ][1, ]
  x0 <- first_point$SampleSize
  y0 <- first_point$Variance
  ref_line <- data.frame(
    SampleSize = sample_sizes,
    Variance = y0 * (sample_sizes / x0)^(-1.8),
    Method = "Slope -1.8"
  )
  
  #combine them
  df_all <- rbind(df_combined, ref_line)
  
 p<-ggplot(
   df_all,
   aes(
     x = SampleSize,
     y = Variance,
     color = Method,
     linetype = Method,
     shape = Method          
   )
 ) +
   geom_point(data = df_combined, size = 3) +
   geom_line() +
   scale_x_continuous(trans = 'log2', breaks = 2^seq(3, 12, by = 1)) +
   scale_y_continuous(
     trans = 'log2',
     breaks = 10^seq(-9, -1, by = 2),
     labels = scales::label_scientific(digits = 3)
   ) +
   
   
 scale_shape_manual(values = c(
   "rho=1" = 19,        # solid circle
   "rho=2" = 17,        # triangle
   "rho=3" = 15,        # square
   "rho=infinity" = 18, # diamond
   "Slope -1.8" = 3     # plus sign
 )) +
   
 
 scale_linetype_manual(values = c(
   "rho=1" = "solid",
   "rho=2" = "solid",
   "rho=3" = "solid",
   "rho=infinity" = "solid",
   "Slope -1.8" = "dashed"
 )) +
   scale_color_manual(values = c(
     "rho=1" = "#F8766D",
     "rho=2" = "#7CAE00",
     "rho=3" = "#00BFC4",
     "rho=infinity" = "#C77CFF",
     "Slope -1.8" = "black"
   )) +
   labs(
     x = "Sample Size",
     y = "Variance",
     color = NULL,
     linetype = NULL,
     shape = NULL            
   ) +
   theme_minimal() +
   theme(
     axis.text = element_text(size = 11),
     axis.title = element_text(size = 12),
     legend.position = c(0.02, 0.02),
     legend.justification = c("left", "bottom"),
     legend.text = element_text(size = 12),
     legend.key.size = unit(0.8, "lines"),
     legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.5),
     legend.margin = margin(6, 6, 6, 6),
     plot.margin = margin(10, 10, 10, 10)
   ) +
   guides(
     color = guide_legend(nrow = 5, byrow = TRUE),
     linetype = guide_legend(nrow = 5, byrow = TRUE),
     shape = guide_legend(nrow = 5, byrow = TRUE)
   )
  print(p)
  
  if (!is.null(file)) {
    ggsave(file, plot = p, width = 8, height = 6, units = "in")#save file
  }
  
  return(p)
}
