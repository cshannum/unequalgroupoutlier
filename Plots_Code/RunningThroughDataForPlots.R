rm(list = ls())
library(tidyverse)
library(ggplot2)
library(scatterplot3d)
library(unequalgroupoutlier)

#folder_location <- "/Users/channum/Documents/Project/Code/Data/" # for NREL comp
folder_location <- "/Users/cshannum/Documents/CSM/NRELCollab/Data/" # for home comp


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
summary(MyIndex)

# append vector to Order_user_power
Order_user_power$length <- MyIndex

# Filter to grab only ones with length greater than 1000
length_index <- which(Job_table > 1000 & Job_table < 6000)
Order_user_power <- Order_user_power

#%>%
 # filter(length > 1000 & length < 6000)

# how many datasets do we have
N <- length(unique(Order_user_power$JobHost)) # 100
RelTime_Sec <- round(Order_user_power$RelTime/1000, digits = 0)

# append to Order_user_power our time vector
Order_user_power$RelTime_Sec <- RelTime_Sec

# Put all into useable matrix form
# Y <- Order_user_power %>%
#   dplyr::select(JobHost, value, RelTime_Sec) %>%
#   spread(key = RelTime_Sec, value = value)

# version for plotting
Y <- Order_user_power %>%
  dplyr::select(JobHost, value, RelTime_Sec)

# create plot but don't plot it
plotY <- Y %>%
  ggplot(., aes(x = RelTime_Sec, y = value, group = JobHost)) +
  geom_line(alpha = .5)


#-------------------------------------------------------------------------------
# Choose knots -----------------------------------------------------------------

# choose knots
K <- 100
P <- 1
Y <- Y %>%
  rename(ID = JobHost, y = value, x = RelTime_Sec)
knots <- choose_knots(Y, K)
# add knots so really short functions have more knots than before
knots_add <- seq(2, 500, len = 50)
knots <- c(knots_add, knots)
# Reduce to only jobs that go across at least 2 knots
# keep_Jobs <- unique_host[-(which(Job_table <= knots[2]))] # will be left with 297 datasets
# N <- length(keep_Jobs)
# Y <- Y %>%
#   dplyr::filter(ID %in% keep_Jobs)

# plot
# par(mfrow = c(1, 1))
# plotY +
#   geom_vline(xintercept = knots, alpha = .1, col = "red")

out <- Mixed_Model_Spline_Fit(Y, knots, theta = 10000)
#
# too_many_knots_rm <- c(11, 12, 15)
ggplot(out$Y, aes(x = x, y = fitted, group = ID)) +
  geom_line(alpha = .5)

save(out, file = paste(folder_location, "out.Rdata"))
#
# P <- 1
# coeff_mat <- matrix(NA, nrow = N, ncol = (K + P + 1))
# for(i in 1:N){
#   coeff_mat[i, ] <- as.numeric(c(out$coeff_fixed[i, -1], out$coeff_random[i, -1]))
# }
#
# pca <- prcomp(t(coeff_mat))
# plot(pca$rotation)
#
# index <- which(pca$rotation[,1] < -.08 & pca$rotation[,2] > .15)
#
# foo <- out$Y %>%
#   dplyr::filter(ID %in% unique_host[index])
#
# ggplot()+
#   geom_line(data = foo, aes(x = x, y = fitted, group = ID))
#
#
# plot(foo$x[which(foo$ID == unique_host[index[1]])],
#      foo$fitted[which(foo$ID == unique_host[index[1]])],
#      type = "l", lwd = 5)
# lines(foo$x[which(foo$ID == unique_host[index[1]])],
#       foo$y[which(foo$ID == unique_host[index[1]])],
#       col ="pink", lwd = 5)
# abline(v = knots)
#
# summary(as.numeric(Job_table))
#
# knots2 <- choose_knots(Y, 1000)
# knots2[1:5]
# length(which(as.numeric(Job_table) <= knots2[2]))
#
# scatterplot3d(pca$rotation[,1:3], type = "h")
# scores <- pca$rotation[,1:3]    # scores for first three PC's
#
# # k-means clustering [assume 3 clusters]
# km     <- kmeans(scores, centers=3, nstart=5)
# ggdata <- data.frame(scores, Cluster=km$cluster)
#
# # stat_ellipse is not part of the base ggplot package
# #source("https://raw.github.com/low-decarie/FAAV/master/r/stat-ellipse.R")
#
# ggplot(ggdata) +
#   scale_colour_manual(values = c("Black", "Blue", "Red", "green")) +
#   geom_point(aes(x=PC1, y=PC2, color=factor(Cluster)), size=5, shape=20) +
#   stat_ellipse(aes(x=PC1,y=PC2,fill=factor(Cluster)),
#                geom="polygon", level=0.95, alpha=0.2) +
#   guides(color=guide_legend("Cluster"),
#          fill=guide_legend("Cluster"))
# outliers1 <- which(km$cluster == 2)
# outliers2 <- which(km$cluster == 3)
#
# ID_outliers1 <- unique(Y$ID)[outliers1]
# ID_outliers2 <- unique(Y$ID)[outliers2]
#
# Y_outliers1 <- out$Y %>%
#   dplyr::select(ID, x, fitted) %>%
#   dplyr::filter(ID %in% ID_outliers1)
# Y_outliers2 <- out$Y %>%
#   dplyr::select(ID, x, fitted) %>%
#   dplyr::filter(ID %in% ID_outliers2)
#
#
# ggplot() +
#   geom_line(data = out$Y, aes(x = x , y = fitted, group = ID)) +
#   geom_line(data = Y_outliers1,
#             aes(x = x, y = fitted, group = ID),
#             color = "red",
#             alpha = 1) +
#   geom_line(data = Y_outliers2,
#             aes(x = x, y = fitted, group = ID),
#             color = "blue",
#             alpha = 1)
#
#
#
#
#
