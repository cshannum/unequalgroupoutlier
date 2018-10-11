rm(list = ls())
library(tidyverse)


folder_location <- "/Users/channum/Documents/Project/Code/Data/" # for NREL comp
#folder_location <- "/Users/cshannum/Documents/CSM/NRELCollab/Data/" # for home comp


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

# version for plotting
Yplot <- Order_user_power %>%
  dplyr::select(JobHost, value, RelTime_Sec)

# create plot but don't plot it
plotY <- Yplot %>%
  ggplot(., aes(x = RelTime_Sec, y = value, group = JobHost)) +
  geom_line(alpha = .5)

#-------------------------------------------------------------------------------
# Choose knots -----------------------------------------------------------------

library(unequalgroupoutlier)
K <- 60
Y <- Yplot %>%
  rename(ID = JobHost, y = value, x = RelTime_Sec)
Y_Sd <- Y %>% group_by(ID) %>% summarise(sd=sd(y))
norm_c <- mean(Y_Sd$sd)
Y$y <- Y$y/norm_c


unique_ID <- unique(Y$ID)
knots <- choose_knots(Y, K)

# plot

basis <- function(x, k, c){exp(-abs(x - k)/c)}
Y_new <- Y %>%
  mutate(basis1 = basis(x, knots[30], theta))
1 - length(which(Y_new$basis1 < .05))/dim(Y_new)[1]

Plot_y <- ggplot(Y_new, aes(x = x, y = y, group = ID)) +
  geom_vline(xintercept = knots, alpha = .5, col = "pink", linetype = "dashed") +
  geom_line(alpha = .5)

Plot_y_basis <- Plot_y +
  geom_line(aes(x = x, y = basis1*75 + 25, group = ID), col = "hotpink2", size = 2) #+
  #geom_vline(xintercept = Y_new$x[which(Y_new$basis1 < .05)]) # if you want to see where the basis applies weight larger than .05
#print(Plot_y_basis)

# Fit the data ####
theta <- sd(knots)*K^(-1/3)
fit_all_values <- Mixed_Model_Spline_Fit(Y, knots, theta = theta)

Plot_fitted <- ggplot(fit_all_values$Y, aes(x = x, y = fitted, group = ID)) +
  geom_vline(xintercept = knots, alpha = .5, col = "pink", linetype = "dashed") +
  geom_line(alpha = .5)
print(Plot_fitted)

# Put coefficients into matrix ####
coeff <- inner_join(fit_all_values$coeff_fixed, fit_all_values$coeff_random, by = "ID")
coeff <- as.matrix(coeff[,-1])

# Bag Plot ####
library(aplpack)
bagplot_outliers <- function(coeff){
  # first do pca on coeff
  pca <- prcomp(coeff, center = F, scale. = F, retx = F)

  # Make sure dkmethod = 1 because we already are using the principal components
  outliers_bag <- list(pxy.outer = NULL)
  t <- 2
  while(length(outliers_bag$pxy.outer) == 0){
    t <- t + 1
    outliers_bag <- compute.bagplot(pca$rotation[, 1:2],
                                    factor = t,
                                    dkmethod = 2)
  }

  outliers_bag <- bagplot(pca$rotation[, 1:2],
                          factor = t,
                          pch = 16,
                          cex = 1,
                          main = "Bivariate Bagplot",
                          dkmethod = 1,
                          show.whiskers = F)
  outliers_bag$factor <- as.list(t)

  outliers <- which(pca$rotation[, 1] %in% outliers_bag$pxy.outlier[, 1])
  outliers_bag$outliers_index <- outliers
  return(outliers_bag)
}

pca <- prcomp(t(coeff), scale. = F, center = F)

# bagplot(pca$rotation[, 1:2],
#         factor = 9,
#         pch = 16,
#         cex = 1,
#         main = "Bivariate Bagplot",
#         dkmethod = 1,
#         show.whiskers = F,
#         show.looppoints = F,
#         show.bagpoints = F)

#bag_out <- bagplot_outliers(coeff = t(coeff))
bag_out <- compute.bagplot(pca$rotation[, 1:2],
                   factor = 9,
                   dkmethod = 1)
plot(bag_out,
        pch = 16,
        cex = 1,
        main = "Bivariate Bagplot")
outliers_bag <- which(pca$rotation[, 1] %in% bag_out$pxy.outlier[, 1])
bag_out$outliers_index <- outliers_bag
unique_host <- unique(Y$ID)
outlier_index <- bag_out$outliers_index
ID_outlier_index <- which(Y$ID %in% unique_host[outlier_index])

col_coding <- rep(1, dim(Y)[1])
col_coding[ID_outlier_index] <- 2

fit_all_values$Y$outliers <- col_coding

ggplot(fit_all_values$Y, aes(x = x, y = fitted, group = ID, col = as.factor(outliers))) +
  geom_line(alpha = .5)

# Rainbow curve ####
library(DepthProc, quietly = T)
N <- length(unique_ID)
depths <- DepthProc::depth(coeff, coeff, method = "LP")

depth_order <- order(depths, decreasing = T)
ID_table <- table(fit_all_values$Y$ID)
depth_col <- c()
index <- c()
loc <- c()
for(i in 1:N){
  J <- as.numeric(ID_table[i])
  depth_col <- append(depth_col, rep(depths[i], J))
  index <- append(index, rep(which(depth_order == i) , J))
  t <- ifelse(i %% 2 == 0, -5, 5)
  loc <- append(loc, rep(t, J))
}

fit_all_values$Y$Depth <- depth_col
depth_rainbow_df <- fit_all_values$Y %>%
  dplyr::select(ID, x, fitted, Depth) %>%
  arrange(Depth)
depth_rainbow_df$order <- apply(depth_rainbow_df[,c("ID", "Depth")], 2, function(x) paste(x, collapse = " "))

ggplot() +
  geom_vline(xintercept = knots, alpha = .1, linetype = "dashed") +
  geom_line(data = depth_rainbow_df,
            aes(x = x,  y = fitted, col = Depth, group = as.factor(desc(Depth))), size = 1) +
  scale_color_gradientn(colours = rainbow(N, start = 0, end = 5/6)) +
  labs(title = "Rainbow Plot of Depth", hjust = .5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5, size= rel(2)))

fit_all_values$Y$index <- index
fit_all_values$Y$Names <- fit_all_values$Y$ID
id_outliers <- fit_all_values$Y$Names %in% as.character(outliers_bag)
fit_all_values$Y$Names[!id_outliers] <- ""

ggplot(fit_all_values$Y, aes(x = index, y = Depth, group = ID, col = Depth)) +
  geom_point() +
  scale_color_gradientn(colours = rainbow(N, start = 0, end = 5/6)) +
  #geom_text(aes(x = index + loc, y = Depth, group = ID, label = Names), size = 4, hjust = .5) +
  labs(title = "LP Depth", hjust = .5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5, size= rel(2)))

# Box Plot of Depths ####

box_depths <- boxplot(depths)
outliers_index <- which(depths %in% box_depths$out)
outliers_ID <- which(fit_all_values$Y$ID %in% unique_host[outliers_index])
outliers_col <- rep(1, dim(Y)[1])
outliers_col[outliers_ID] <- 2
fit_all_values$Y$depth_outliers <- outliers_col

ggplot(fit_all_values$Y, aes(x = x, y = fitted, group = ID, col = as.factor(depth_outliers)))+
  geom_line(alpha = .5)

# Adjusted Box plot for Skewed Data ####
library(robustbase)

adj_box_depth <- adjbox(depths, range = .4)
outliers_adj <- which(depths %in% adj_box_depth$out)
outliers_ID <- which(fit_all_values$Y$ID %in% unique_host[outliers_adj])
outliers_col <- rep(1, dim(Y)[1])
outliers_col[outliers_ID] <- 2
fit_all_values$Y$depth_outliers_adj <- outliers_col

ggplot(fit_all_values$Y, aes(x = x, y = fitted, group = ID, col = as.factor(depth_outliers_adj)))+
  geom_line(alpha = .5)


# HDR 1D ####
library(hdrcde)
library(plotrix)

cluster <- kmeans(coeff, centers = 4)
cluster_col <- cluster$cluster

plot(pca$rotation[, 1:2], col = cluster_col)

hdr_1d <- hdr(x = depths, lambda = 1, all.modes = T)
hdr_density <- hdrcde::hdr.den(depths, plot.lines = T, prob = c(50, 90, 95), col = c("blue", "green", "red"), main = "HDR 1D")
points(sort(depths), rep(0, length(depths)), pch = 16)
hdr_density$hdr
#points(hdr_density$hdr[1, ], rep(10, 22), col = "magenta", pch = 8)

new_depths <- (depths - mean(depths))/sd(depths)
bw <- .5
den <- density(new_depths, bw = bw)
hdr_density <- hdrcde::hdr.den(den = den,
                               plot.lines = T,
                               prob = c(50, 90, 95),
                               col = c("blue", "green", "red"))
outliers <- c()
num_box <- sum(!is.na(hdr_density$hdr[1, ]))/2
for(i in 1:num_box){
  outliers <- append(outliers,
                     which(new_depths < hdr_density$hdr[1, i] |
                             new_depths > hdr_density$hdr[1, (i+1)]))
}

which(depths <= quantile(depths, prob = .01))

col_simple <- rep("black", 272)
col_simple[outliers] <- "red"
points(new_depths[order(depths)], seq(0, .3, len = 272), col = col_simple[order(depths)], pch = 16)

depth_outliers <- rep(1, dim(fit_all_values$Y)[1])
outlier_index <- which(fit_all_values$Y$ID %in% unique_ID[outliers])
depth_outliers[outlier_index] <- 2
fit_all_values$Y$depth_outliers <- depth_outliers
Plot_depth_outliers <- fit_all_values$Y %>%
  select(x, fitted, ID, depth_outliers)
sub_plot <- Plot_depth_outliers %>%
  filter(depth_outliers == 2)
ggplot()+
  geom_line(data = Plot_depth_outliers,
            aes(x = x, y = fitted, group = ID,
                col = as.factor(depth_outliers)),
            alpha = .5) +
  scale_color_manual(values = c("black", "red")) +
  geom_line(data = sub_plot,
            aes(x = x, y = fitted, group = ID), col = "red") +
  labs(title = "Outliers from HDR of Depths")



# HDR Boxplot 2D ####


hdr_coeff <- hdr.2d(pca$rotation[, 1], pca$rotation[, 2], prob = c(50, 95, 99))
#hdr_coeff$fxy
outliers <- which(hdr_coeff$fxy < hdr_coeff$falpha[2])


hdr.boxplot.2d(pca$rotation[, 1], pca$rotation[, 2],
               prob = c(50, 95, 99),
               show.points = t,
               main = "HDR Bivariate Boxplot",
               xlab = "PC1",
               ylab = "PC2",
               shadecols = c(gray((3:1)/(3 + 1))))
points(pca$rotation[outliers, ],
       pch = 8,
       col = "green3",
       lwd = 5)
#es <- emptyspace(pca$rotation[,1], pca$rotation[, 2])
legend(es,
       c("50%", "95%", "99%", "Identified Outliers"),
       col = c(gray((3:1)/(3 + 1)), "green3"),
       pch = c(15, 15, 15, 8), pt.cex = c(1.5, 1.5, 1.5, 1), pt.lwd = 5)


outliers_hdr <- outliers
outliers_ID <- which(fit_all_values$Y$ID %in% unique_host[outliers_hdr])
outliers_col <- rep(1, dim(Y)[1])
outliers_col[outliers_ID] <- 2
fit_all_values$Y$hdr_outliers <- outliers_col


# All plots ####

ggplot(fit_all_values$Y, aes(x = x, y = fitted, group = ID, col = as.factor(outliers)))+
  geom_line(alpha = .5) +
  scale_color_manual(values = c("black", "red")) +
  labs(title = "Outliers from Bagplot")

ggplot(fit_all_values$Y, aes(x = x, y = fitted, group = ID, col = as.factor(depth_outliers)))+
  geom_line(alpha = .5) +
  scale_color_manual(values = c("black", "red")) +
  labs(title = "Outliers from Boxplot of Depths")

ggplot(fit_all_values$Y, aes(x = x, y = fitted, group = ID, col = as.factor(depth_outliers_adj)))+
  geom_line(alpha = .5) +
  scale_color_manual(values = c("black", "red")) +
  labs(title = "Outliers from Adjusted Boxplot of Depths")

ggplot(fit_all_values$Y, aes(x = x, y = fitted, group = ID, col = as.factor(hdr_outliers)))+
  geom_line(alpha = .5) +
  scale_color_manual(values = c("black", "red")) +
  labs(title = "Outliers from HDR")

ggplot(fit_all_values$Y, aes(x = index, y = Depth, group = ID, col = Depth)) +
  geom_point() +
  scale_color_gradientn(colours = rainbow(N, start = 0, end = 5/6)) +
  #geom_text(aes(x = index + loc, y = Depth, group = ID, label = Names), size = 4, hjust = .5) +
  labs(title = "LP Depth", hjust = .5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5, size= rel(2)))

ggplot() +
  geom_vline(xintercept = knots, alpha = .1, linetype = "dashed") +
  geom_line(data = fit_all_values$Y,
            aes(x = x,  y = fitted, group = ID, col = Depth), size = 1) +
  scale_color_gradientn(colours = rainbow(N, start = 0, end = 5/6)) +
  labs(title = "Rainbow Plot of Depth", hjust = .5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5, size= rel(2)))





