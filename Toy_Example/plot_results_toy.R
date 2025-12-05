
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

plot_variance <- function(var_mc_list,var_plain_rqmc_list,var_rqmc_pow2_list,var_rqmc_adjusted_list,var_rqmc_indep_list, sample_sizes, file = NULL) {
  
  # Construct all dataframes
  df_mc <- data.frame(SampleSize = sample_sizes, Variance = unlist(var_mc_list), Method = "MC")
  df_qmc <- data.frame(SampleSize = sample_sizes, Variance = unlist(var_plain_rqmc_list), Method = "RQMC")
  df_qmc_alloc <- data.frame(SampleSize = sample_sizes, Variance = unlist(var_rqmc_pow2_list), Method = "RQMC powers of 2")
  df_rqmc_adjusted <- data.frame(SampleSize = sample_sizes, Variance = unlist(var_rqmc_adjusted_list), Method = "RQMC adjusted")
  df_rqmc_indep <- data.frame(SampleSize = sample_sizes, Variance = unlist(var_rqmc_indep_list), Method = "RQMC L independent")
  
  # Reference lines
  C1 <- var_mc_list[[1]] * sample_sizes[1]
  df_ref1 <- data.frame(SampleSize = sample_sizes, Variance = C1 * sample_sizes^(-1), Method = "Slope = -1")
  
  C2 <- var_plain_rqmc_list[[10]] * sample_sizes[10]^2
  df_ref2 <- data.frame(SampleSize = sample_sizes, Variance = C2 * sample_sizes^(-2), Method = "Slope = -2")
  
  C3 <- var_rqmc_pow2_list[[10]] * sample_sizes[10]^3
  df_ref3 <- data.frame(SampleSize = sample_sizes, Variance = C3 * sample_sizes^(-3), Method = "Slope = -3")
  
  # Combine all
  df_all <- rbind(df_mc, df_qmc, df_qmc_alloc, df_rqmc_adjusted, df_rqmc_indep, df_ref1, df_ref2, df_ref3)
  method_levels <- c(
    "MC",
    "RQMC",
    "RQMC adjusted",
    "RQMC L independent",
    "RQMC powers of 2",
    "Slope = -1",
    "Slope = -2",
    "Slope = -3"
  )
  df_all$Method <- factor(df_all$Method, levels = method_levels)
  
  # Plot
  p <- ggplot(df_all, aes(x = SampleSize, y = Variance, group = Method)) +
    # Main curves with points
    geom_point(
      data = subset(df_all, !grepl("^Slope", Method)),
      aes(color = Method, shape = Method),
      size = 3
    ) +
    geom_line(
      data = subset(df_all, !grepl("^Slope", Method)),
      aes(color = Method)
    ) +
    # Reference slope lines
    geom_line(
      data = subset(df_all, grepl("^Slope", Method)),
      aes(color = Method, linetype = Method),
      linewidth = 1.2
    ) +
    
    # Scales
    scale_x_continuous(trans = 'log2', breaks = 2^seq(3, 12, by = 1)) +
    scale_y_continuous(
      trans = 'log2',
      breaks = 10^seq(-9, -1, by = 2),
      labels = scales::label_scientific(digits = 3)
    ) +
    scale_color_manual(values = c(
      "MC" = "blue",
      "RQMC" = "chartreuse4",
      "RQMC adjusted" = "brown4",
      "RQMC L independent" = "darkviolet",
      "RQMC powers of 2" = "deeppink",
      "Slope = -1" = "steelblue",
      "Slope = -2" = "green",
      "Slope = -3" = "pink"
    )) +
    scale_linetype_manual(values = c(
      "Slope = -1" = "dashed",
      "Slope = -2" = "dashed",
      "Slope = -3" = "dashed"
    )) +
    scale_shape_manual(values = c(
      "MC" = 19,                 # solid circle
      "RQMC" = 17,               # triangle
      "RQMC adjusted" = 8,      # square
      "RQMC L independent" = 18, # diamond
      "RQMC powers of 2" = 15     # plus
    )) +
    
    # Labels and guides
    labs(
      x = "Sample Size",
      y = "Variance",
      color = NULL,
      shape = NULL,
      linetype = NULL
    ) +
    guides(
      color = guide_legend(override.aes = list(
        linetype = c("solid", "solid", "solid", "solid", "solid", "dashed", "dashed", "dashed"),
        shape = c(19, 17, 8, 18, 15, NA, NA, NA)
      )),
      shape = "none",   # remove duplicate legend
      linetype = "none"
    ) +
    
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
    )
  print(p)
  if (!is.null(file)) {
    ggsave(file, plot = p, width = 8 , height = 6, units = "in") #save file
  }
  
  return(p)
  
  
}