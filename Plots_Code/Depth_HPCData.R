rm(list = ls())
library(tidyverse)
library(DepthProc, quietly = T)
library(unequalgroupoutlier)

#folder_location <- "/Users/channum/Documents/Project/Code/Data/" # for NREL comp
folder_location <- "/Users/cshannum/Documents/CSM/NRELCollab/Data/" # for home comp

options(warn = 2)
#-------------------------------------------------------------------------------
# Access metadata to get jobs of interest --------------------------------------

load(paste(folder_location, "10k_jobs_metadata.RData", sep = ""))
load(file = paste(folder_location, "User_Power_Data.RData", sep = ""))
# load(file = paste(folder_location, "Y.Rdata", sep=""))
# load(file = paste(folder_location, "x.Rdata", sep=""))

user_jobs <- df_jobs %>%
  filter(user_name == "c49b64ce0319f7718fa1af317745dfd4e1771331") #some annomous user

#-------------------------------------------------------------------------------
# # Get filenames for data -----------------------------------------------------
#
# #notice that I append the filepath to point to where the files are stored
# tic <- proc.time()
# user_power <- data_frame(File_Short = unique(user_jobs$File_Short)) %>%
#   #mutate(File = file.path( "/Users/channum/Documents/Project/Code/Data/10k_processed_data", File_Short)) %>%
#   mutate(File = file.path(paste(folder_location, "10k_processed_data", sep = ""), File_Short)) %>%
#   split(.$File) %>%
#   map_df(.f = function(z) {
#     cur_jobs <- user_jobs %>%
#       filter(File_Short == z$File_Short)
#
#     load(z$File)
#     cur_power <- df.power %>%
#       filter(Job %in% cur_jobs$Job)
#     rm(list = c("df.ganglia", "df.power", "df.job"))
#     return(cur_power)
#   })
#
# toc <- proc.time()
# elapsed_time <- toc - tic
# elapsed_time
#
# #save(user_power, file = paste(folder_location, "User_Power_Data.RData", sep = ""))

#-------------------------------------------------------------------------------
# # Take a look at lengths of jobs ---------------------------------------------
# #Rel time is measured in milliseconds
# #RelTime.Min is measured in minutes
#
# user_runtimes <- user_power %>%
#   group_by(Job) %>%
#   summarise(Runtime_Min = diff(range(RelTime.Min)))
#
# user_runtimes %>%
#   ggplot(., aes(x = Runtime_Min)) +
#   geom_histogram(fill = "dodgerblue", colour = "black", binwidth = 300) +
#   labs(x = "Runtime (Minutes)") +
#   theme_bw()
#
# user_power %>%
#   ggplot(., aes(x = RelTime.Min, y = value, group = JobHost)) +
#   geom_line(alpha = .1)

#-------------------------------------------------------------------------------
# Put into usable format -------------------------------------------------------

# identify all unique Job's
unique_host <- unique(user_power$JobHost)

# Put JobHost, value and RelTime in order by Jobhost then RelTime
Order_user_power <- user_power %>%
  dplyr::select(JobHost, value, RelTime) %>%
  arrange(JobHost, RelTime)

# How many of each unique job there is
Job_table <- table(Order_user_power$JobHost)
# create vector containing length of each unique job dataset
MyIndex <- c()
for(i in 1:length(Job_table)){
  MyIndex <- append(MyIndex, rep(as.numeric(Job_table[i]), Job_table[i]))
}

# append vector to Order_user_power
Order_user_power$length <- MyIndex

# Filter to grab only ones with length greater than 1000
length_index <- which(Job_table > 1000 & Job_table < 6000)
Order_user_power <- Order_user_power %>%
  filter(length > 1000)

# how many datasets do we have
N <- length(unique(Order_user_power$JobHost)) # 100

RelTime_Sec <- round(Order_user_power$RelTime/1000, digits = 0)

# append to Order_user_power our time vector
Order_user_power$RelTime_Sec <- RelTime_Sec

# Final version
Y_df <- Order_user_power %>%
  dplyr::select(JobHost, value, RelTime_Sec)

#-------------------------------------------------------------------------------
# Choose knots -----------------------------------------------------------------

# Change columns names for input
Y <- Y_df %>%
  rename(ID = JobHost, y = value, x = RelTime_Sec)

# scale response variable so variance is closer to 1
Y_Sd <- Y %>% group_by(ID) %>% summarise(sd=sd(y))
norm_sd <- mean(Y_Sd$sd)
Y$y <- Y$y/norm_sd

# For later use
unique_ID <- unique(Y$ID)

# Choose knots
K <- 60
knots <- choose_knots(Y, K)
# Compute range parameter
theta <- 3*sd(knots)*K^(-1/3)

# plot data and basis
basis <- function(x, k, c){exp(-abs(x - k)/c)} # radial basis
Y_basis <- Y %>%
  mutate(basis1 = basis(x, knots[2], theta)) %>% # choose middle knots
  mutate(basis2 = basis(x, knots[20], theta)) %>%
  mutate(basis3 = basis(x, knots[30], theta)) %>%
  mutate(basis4 = basis(x, knots[40], theta))

png("/Users/cshannum/Documents/CSM/NRELCollab/Package/unequalgroupoutlier/Plots_Code/CPU_Power_Plot.png")
Plot_y <- ggplot(Y_basis, aes(x = x, y = y, group = ID)) +
  #geom_vline(xintercept = knots, alpha = .3, linetype = "dashed") +
  geom_line(alpha = .5) +
  labs(title = "CPU Power", hjust = .5, shape = "Type") +
  theme_bw() +
  theme(panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = .5, size = ggplot2::rel(2)))
print(Plot_y)
dev.off()

# png("/Users/cshannum/Documents/CSM/NRELCollab/Package/unequalgroupoutlier/Plots_Code/CPU_Power_Basis_Plot.png")
# Plot_y_basis <- ggplot() +
#   geom_vline(xintercept = knots, alpha = .3, linetype = "dashed") +
#   geom_line(data = Y_basis, aes(x = x, y = basis1*100, group = ID), col = "deeppink3", size = 2) +
#   geom_line(data = Y_basis, aes(x = x, y = basis2*100, group = ID), col = "blue4", size = 2) +
#   geom_line(data = Y_basis, aes(x = x, y = basis3*100, group = ID), col = "chocolate2", size = 2) +
#   geom_line(data = Y_basis, aes(x = x, y = basis4*100, group = ID), col = "chartreuse4", size = 2) +
#     labs(title = "Radial Basis", hjust = .5, shape = "Type") +
#   theme_bw() +
#   theme(panel.grid.major = ggplot2::element_blank(),
#         panel.grid.minor = ggplot2::element_blank(),
#         plot.title = ggplot2::element_text(hjust = .5, size = ggplot2::rel(2)))
# print(Plot_y_basis)
# dev.off()

#-------------------------------------------------------------------------------
# Fit the data ####

fit_all_values <- Mixed_Model_Spline_Fit(Y, knots, theta = theta)

save(fit_all_values, file = "/Users/cshannum/Documents/CSM/NRELCollab/Package/unequalgroupoutlier/Plots_Code/FitData.Rdata")
fit_all_values <- load(file = "/Users/cshannum/Documents/CSM/NRELCollab/Package/unequalgroupoutlier/Plots_Code/FitData.Rdata")

coeff <- inner_join(fit_all_values$coeff_fixed, fit_all_values$coeff_random, by = "ID")

coeff

# png("/Users/cshannum/Documents/CSM/NRELCollab/Package/unequalgroupoutlier/Plots_Code/CPU_Power_Fitted_Plot.png")
# Plot_fitted <- ggplot() +
#   geom_vline(xintercept = knots, alpha = .3, linetype = "dashed") +
#   geom_line(data = fit_all_values$Y, aes(x = x, y = fitted, group = ID), alpha = .5) +
#   labs(title = "Fitted CPU Power", hjust = .5, shape = "Type") +
#   theme_bw() +
#   theme(panel.grid.major = ggplot2::element_blank(),
#         panel.grid.minor = ggplot2::element_blank(),
#         plot.title = ggplot2::element_text(hjust = .5, size = ggplot2::rel(2)))
# print(Plot_fitted)
# dev.off()

#-------------------------------------------------------------------------------
# Put coefficients into matrix ####
coeff <- inner_join(fit_all_values$coeff_fixed, fit_all_values$coeff_random, by = "ID")

out <- outliers_plots(Y = fit_all_values$Y, coeff = coeff)
print(out$outliers_ID)

png("/Users/cshannum/Documents/CSM/NRELCollab/Package/unequalgroupoutlier/Plots_Code/CPU_Power_Depth_Plot.png")
print(out$depth_plot)
dev.off()

png("/Users/cshannum/Documents/CSM/NRELCollab/Package/unequalgroupoutlier/Plots_Code/CPU_Power_Rainbow_Plot.png")
print(out$rainbow_plot)
dev.off()

png("/Users/cshannum/Documents/CSM/NRELCollab/Package/unequalgroupoutlier/Plots_Code/CPU_Power_Outlier_Plot.png")
print(out$outliers_plot)
dev.off()

# #-------------------------------------------------------------------------------
# # Rainbow curve ####
# N <- length(unique_ID)
# depths <- DepthProc::depth(coeff, coeff, method = "LP", pdim = 2)
# plot(sort(depths))
#
# depth_order <- order(depths, decreasing = T)
# ID_table <- table(fit_all_values$Y$ID)
# depth_col <- c()
# index <- c()
# loc <- c()
# for(i in 1:N){
#   J <- as.numeric(ID_table[i])
#   depth_col <- append(depth_col, rep(depths[i], J))
#   index <- append(index, rep(which(depth_order == i) , J))
#   t <- ifelse(i %% 2 == 0, -5, 5)
#   loc <- append(loc, rep(t, J))
# }
#
# fit_all_values$Y$Depth <- depth_col
# depth_rainbow_df <- fit_all_values$Y %>%
#   dplyr::select(ID, x, fitted, Depth) %>%
#   arrange(Depth)
# depth_rainbow_df$order <- apply(depth_rainbow_df[,c("ID", "Depth")], 2, function(x) paste(x, collapse = " "))
#
# ggplot() +
#   geom_vline(xintercept = knots, alpha = .1, linetype = "dashed") +
#   geom_line(data = depth_rainbow_df,
#             aes(x = x,  y = fitted, col = Depth, group = ID)) +
#   scale_color_gradientn(colours = rainbow(N, start = 0, end = 5/6)) +
#   labs(title = "Rainbow Plot of Depth", hjust = .5) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         plot.title = element_text(hjust = .5, size= rel(2)))
#
# # ------------------------------------------------------------------------------
# # HDR --------------------------------------------------------------------------
#
# hdr_density <- hdrcde::hdr.den(depths,
#                                plot.lines = T,
#                                prob = c(50, 90, 95),
#                                col = c("dodgerblue3", "chartreuse3", "firebrick3"),
#                                main = "HDR 1D: all outliers",
#                                h = sd(depths)*length(depths)^(-1/3))
# points(sort(depths), rep(0, length(depths)),
#        col = outliers_col[order(depths)],
#        pch = pch_type[order(depths)])
#
# num_boxes <- dim(hdr_density$hdr)[2] -1
# index_seq <- seq(1, num_boxes, 2)
# boxes <- c(-Inf, hdr_density$hdr[1, ])
# outliers <- c()
# for(i in index_seq){
#   outliers <- append(outliers,
#                      which(depths >= boxes[i] &
#                              depths <= boxes[i + 1]))
# }
#
# loc <- lapply(outliers, function(x) which(x == order(depths)))
# loc <- unlist(loc)
# loc <- sort(loc)
# diff_loc <- diff(loc)
# if(any(diff_loc > 5)){
#   t <- which(diff_loc > 5)
#   low <- loc[t] + 1
#   up <- loc[t + 1] - 1
#   group <- c(low:up)
#   loc_group <- order(depths)[group]
# }
#
# outliers_df <- depth_rainbow_df %>%
#   dplyr::select(ID, x, fitted) %>%
#   mutate(outliers = 1)
# index <- which(outliers_df$ID %in% unique_ID[outliers])
# index2 <- which(outliers_df$ID %in% unique_ID[loc_group])
# outliers_df$outliers[index] <- 3
# outliers_df$outliers[index2] <- 2
#
# ordered_names <- paste(outliers_df$outliers, outliers_df$ID, sep = "")
# outliers_df$order <- ordered_names
#
# ggplot() +
#   geom_vline(xintercept = knots, col = "pink", linetype = "dashed", alpha = .1) +
#   geom_line(data = outliers_df, aes(x = x, y = fitted, group = order, color = as.factor(outliers))) +
#   scale_colour_manual(values = c("black", "dodgerblue3", "firebrick3"))
#
# outliers_col <- rep("black", N)
# outliers_col[outliers] <- "firebrick3"
# outliers_col[loc_group] <- "dodgerblue3"
# pch_type <- rep(16, N)
# pch_type[outliers] <- 8
#
# pca <- prcomp(t(coeff), center = F, scale. = F)
# plot(pca$rotation[, 1:2], col = outliers_col, pch = pch_type)
# plot(depths[order(depths)], col = outliers_col[order(depths)], pch = pch_type[order(depths)])
#
# # ------------------------------------------------------------------------------
# # kMeans clustering ------------------------------------------------------------
#
# cluster <- kmeans(depths, centers = 3)
# cluster$cluster
#
# plot(depths[order(depths)],
#      col = cluster$cluster[order(depths)],
#      pch = pch_type[order(depths)])
# plot(pca$rotation[, 1:2], col = cluster$cluster, pch = pch_type)
# # ------------------------------------------------------------------------------
# # Density Based Cluster --------------------------------------------------------
# library(dbscan)
#
# depth_diff <- diff(sort(depths))
#
# dbcluster <- dbscan(as.matrix(depths),
#                     eps = quantile(depth_diff, probs = .95),
#                     minPts = 3)
# max(dbcluster$cluster)
# plot(depths[order(depths)],
#      col = dbcluster$cluster[order(depths)],
#      pch = pch_type[order(depths)])
#
# db_color <- c()
# for(i in 1:N){
#   J <- as.numeric(ID_table[i])
#   db_color <- append(db_color, rep(dbcluster$cluster[i], J))
# }
#
# db_plot_df <- depth_rainbow_df %>%
#   dplyr::select(ID, x, fitted)
# db_plot_df$db_color <- db_color
#
# ggplot() +
#   geom_vline(xintercept = knots, alpha = .1, linetype = "dashed") +
#   geom_line(data = db_plot_df,
#             aes(x = x,  y = fitted, col = as.factor(db_color), group = ID)) +
#               scale_color_manual(values = c("#030303", "#8B6914",
#                                             "#228B22", "#CD2626",
#                                             "#1874CD", "#FF1493",
#                                             "#EE7600", "#008B8B",
#                                             "#EEEE00")) +
#               labs(title = "Different db clusters", hjust = .5) +
#               theme_bw() +
#               theme(panel.grid.major = element_blank(),
#                     panel.grid.minor = element_blank(),
#                     plot.title = element_text(hjust = .5, size= rel(2)))
#

##########









