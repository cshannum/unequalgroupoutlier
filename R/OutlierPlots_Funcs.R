# add knots lines to plots as the enter back into Mixed_Model_Spline_Fit function

outliers_plots <- function(Y, coeff){

  # remove ID from matrix
  coeff_only <- coeff[, -1]
  # Find unique ID for identifying outliers
  unique_ID <- unique(coeff$ID)
  N <- length(unique_ID)
  # count number of each ID fo creating depth and index vectors
  count_ID <- as.numeric(table(Y$ID))
  depths <- DepthProc::depth(coeff_only, coeff_only, method = "LP")
  depths <- as.numeric(depths@.Data)

  # used for plotting
  d_order <- order(order(depths))
  # so depths are plotted in increasing order
  index <- rep(d_order, count_ID)
  depth_col <- rep(depths, count_ID)

  # add to Y dataframe
  Y$index <- index
  Y$Depth <- depth_col
  Y$grouping <- paste(abs(log(Y$Depth)) , Y$ID, sep = "")

  # rainbow plot of fits - colors depend on depths
  rb_plot <- rainbowplot(Y, N)

  # Determing outliers
  mo <- mrfDepth::outlyingness(robustbase::fullRank(coeff_only))
  outliers <- which(!mo$flagX)

  outlier_index <- which(Y$ID %in% unique_ID[outliers])
  outlier_col <- rep("Normal", dim(Y)[1])
  outlier_col[outlier_index] <- "Anomaly"

  Y$outlier <- outlier_col

  # rainbow plot of depth values in ascending order with outliers as different point shape
  rb_depth_plot <- rainbowplot.depth(Y, N)

  # plots 10 "deepest" and then the anomalies
  M <- ifelse(N < 10, N, 10)
  Mdeepest <- tail(sort(depths), tail = M)
  small_df <- Y %>%
   dplyr::filter(Y$Depth %in% Mdeepest | Y$outlier == "Anomaly")

  outlier_plot <- plot_Mdeep_outliers(small_df, M)

  plot_list <- list("rainbow_plot" = rb_plot,
                    "depth_plot" = rb_depth_plot,
                    "outliers_plot" = outlier_plot,
                    "outliers_ID" = unique_ID[outliers])

  gridExtra::grid.arrange(rb_plot, rb_depth_plot)

  return(plot_list)
}

rainbowplot.depth <- function(Y, N){

  rb_depth <-ggplot2::ggplot(Y, ggplot2::aes(x = index, y = Depth, group = ID, col = Depth, shape = as.factor(outlier))) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_gradientn(colours = grDevices::rainbow(N, start = 0, end = 5/6), guide = F) +
    ggplot2::scale_shape_manual(values = c(3, 16)) +
    ggplot2::labs(title = "L2 Depth", hjust = .5, shape = "Type") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(hjust = .5, size = ggplot2::rel(2)))

  return(rb_depth)

}

rainbowplot <- function(Y, N){

  rb_plot <- ggplot2::ggplot() +
    ggplot2::geom_line(data = Y,
                       ggplot2::aes(x = x,  y = fitted, group = grouping, col = Depth), size = 1) +
    ggplot2::scale_color_gradientn(colours = grDevices::rainbow(N, start = 0, end = 5/6)) +
    ggplot2::labs(title = "Rainbow Plot of L2 Depth", hjust = .5) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(hjust = .5, size= rel(2)))

  return(rb_plot)

}

plot_Mdeep_outliers <- function(Y, M){
  title <- paste(M, "Deepest Curves and All Outlier Curves")

  o_plot <- ggplot2::ggplot(Y, ggplot2::aes(x = x, y = fitted, group = grouping, col = as.factor(outlier))) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = c("blue3", "firebrick3")) +
    ggplot2::labs(title = title, hjust = .5, color = "Type") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(hjust = .5, size= ggplot2::rel(2)))
}
