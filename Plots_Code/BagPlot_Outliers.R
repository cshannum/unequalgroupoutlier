# This is an example script

rm(list = ls())
library(unequalgroupoutlier)
library(tidyverse)
library(plotrix)
library(knitr)

# get toy data ####
ToyData <- Unequal_ToyData(N = 1000)
dat <- ToyData$data #+ c(runif(100, 0, .5), runif(5, .5, .7), runif(5, -.5, 0))
x <- ToyData$x
endpoints <- ToyData$endpoints
mycol <- c(rep("black", 100), rep("red", 5), rep("blue", 5))

# plot toy data
a <- min(dat, na.rm = T)
b <- max(dat, na.rm = T)
plot(x, dat[1, ], type = "l", ylim = c(a, b))
for(i in 1:110){
  lines(x, dat[i, ], type ="l", col = mycol[i])
}

y <- c()
newx <- c()
ID <- c()
mycol <-  c()
for(i in 1:100){
  y <- append(y, dat[i, 1:endpoints[i] ])
  newx <- append(newx, x[1:endpoints[i]])
  ID <- append(ID, rep(i, endpoints[i]))
  mycol <- append(mycol, rep("black", endpoints[i]))
}
for(i in 101:105){
  y <- append(y, dat[i, 1:endpoints[i] ])
  newx <- append(newx, x[1:endpoints[i]])
  ID <- append(ID, rep(i, endpoints[i]))
  mycol <- append(mycol, rep("red", endpoints[i]))
}
for(i in 106:110){
  y <- append(y, dat[i, 1:endpoints[i] ])
  newx <- append(newx, x[1:endpoints[i]])
  ID <- append(ID, rep(i, endpoints[i]))
  mycol <- append(mycol, rep("blue", endpoints[i]))
}

Y <- tibble(ID = ID, x = newx, y = y)


# Determine knots ####
K <- 60
knots <- choose_knots(Y, K = K)
knots <- seq(0, 10, len = 80)
abline(v = knots, lty = 2)

# Compute the fit for radial basis ####
fit_all_R <- Mixed_Model_Spline_Fit.Keep_all_Knots(Y, knots = knots, theta = 2)

# plot the fitted data ####
fit_all_R$Y$col <- mycol
Data <- tibble(x = ToyData$x, fx = ToyData$f_x, fx.A1 = ToyData$f_x.A1, fx.A2 = ToyData$f_x.A2)
plot_fit_data <- ggplot() +
  geom_line(data = fit_all_R$Y,
            aes(x = x,  y = fitted, group = ID, col = as.factor(col)), size = 2) +
  scale_color_manual(values = c("black", "blue", "red")) +
  geom_line(data = Data,
            aes(x = x, y = fx), size = 1, col = "grey", linetype = "dashed") +
  geom_line(data = Data,
            aes(x = x, y = fx.A1), size = 1, col = "pink", linetype = "dashed") +
  geom_line(data = Data,
            aes(x = x, y = fx.A2), size = 1, col = "lightblue", linetype = "dashed") +
  geom_vline(xintercept = knots,
             linetype = "dashed",
             alpha = .1) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot_fit_data

# basis <- function(x, k, c){exp(-abs(x-k)^2/c)}
# fit_all_R$Y <- fit_all_R$Y %>%
#   mutate(basis1 = basis(x, knots[80], 6))
# plot_fit_data + geom_line(data = fit_all_R$Y, aes(x = x, y = basis1), col = "magenta") +
#   geom_vline(xintercept =fit_all_R$Y$x[which(as.numeric(fit_all_R$Y$basis1) < .01)])
#
# 1 - length(which(as.numeric(fit_all_R$Y$basis1) < .01))/dim(fit_all_R$Y)[1]

# create coefficients matrix ####

coeff <- inner_join(fit_all_R$coeff_fixed, fit_all_R$coeff_random, by = "ID")
coeff <- as.matrix(coeff)[,-1]

# Bagplot ####
library(aplpack)

# first do pca on coeff
pca <- prcomp(t(coeff), center = F, scale. = F, retx = F)

# Make sure dkmethod = 1 because we already are using the principal components
outliers_bag <- bagplot(pca$rotation[, 1:2],
                        factor = 6,
                        pch = 16,
                        cex = 1,
                        main = "Bivariate Bagplot", dkmethod = 2)
outliers <- which(pca$rotation[, 1] %in% outliers_bag$pxy.outlier[, 1])
plot(outliers_bag$hdepths[order(outliers_bag$hdepths)], col = rainbow(110))


mycol <-  c()
for(i in 1:110){
  if(i %in% outliers){
    mycol <- append(mycol, rep("green3", endpoints[i]))}else{
    mycol <- append(mycol, rep("black", endpoints[i]))
  }
}

fit_all_R$Y$col <- mycol
ggplot() +
  geom_line(data = fit_all_R$Y,
            aes(x = x,  y = fitted, group = ID, col = as.factor(col)), size = 2) +
  scale_color_manual(values = c("black", "green3"))


# Rainbow curve ####
library(DepthProc)

depths <- DepthProc::depth(coeff, coeff, method = "LP")
graphics_col <- c(rep("black", 100), rep("red", 5), rep("blue", 5))
#plot(depths[order(depths)], col = brewer.pal(n=11, "Paired") )

box_depths <- boxplot(depths)
which(depths %in% box_depths$out)
box_depths$names
col_simple <- rep("black", 110)
col_simple[101:110] <- "red"
plot(depths[order(depths)],
     col = col_simple[order(depths)],
     pch = 16, ylab = "Depth")

depth_order <- order(depths, decreasing = T)
mycol <-  c()

depth_col <- c()
index <- c()
loc <- c()
for(i in 1:110){
  depth_col <- append(depth_col, rep(depths[i], endpoints[i]))
  index <- append(index, rep(which(depth_order == i) , endpoints[i]))
  t <- ifelse(i %% 2 == 0, -5, 5)
  loc <- append(loc, rep(t, endpoints[i]))
}

fit_all_R$Y$Depth <- depth_col
ggplot() +
  geom_vline(xintercept = knots, alpha = .1, linetype = "dashed") +
  geom_line(data = fit_all_R$Y,
            aes(x = x,  y = fitted, group = ID, col = Depth), size = 1) +
  scale_color_gradientn(colours = rainbow(110, start = 0, end = 5/6)) +
  labs(title = "Rainbow Plot of Depth", hjust = .5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5, size= rel(2)))

fit_all_R$Y$index <- index
fit_all_R$Y$Names <- fit_all_R$Y$ID
id_outliers <- fit_all_R$Y$Names %in% as.character(outliers)
fit_all_R$Y$Names[!id_outliers] <- ""

ggplot(fit_all_R$Y, aes(x = index, y = Depth, group = ID, col = Depth)) +
  geom_point() +
  scale_color_gradientn(colours = rainbow(110, start = 0, end = 5/6)) +
  #geom_text(aes(x = index + loc, y = Depth, group = ID, label = Names), size = 4, hjust = .5) +
  labs(title = "LP Depth", hjust = .5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5, size= rel(2)))



# HDR Boxplot ####
library(hdrcde)
hdr_coeff <- hdr.2d(pca$rotation[, 1], pca$rotation[, 2], prob = c(50, 95, 99))
#hdr_coeff$fxy
outliers <- which(hdr_coeff$fxy < hdr_coeff$falpha[2])


hdr.boxplot.2d(pca$rotation[,1], pca$rotation[, 2],
               prob = c(50, 95, 99),
               show.points = t,
               main = "HDR Bivariate Boxplot",
               xlab = "PC1",
               ylab = "PC2")
points(pca$rotation[outliers, ],
       pch = 8,
       col = "green3",
       lwd = 5)
points(pca$rotation,
       pch = 16,
       col = c(rep("black", 100), rep("red", 10)),
       lwd = 5)
es <- emptyspace(pca$rotation[,1], pca$rotation[, 2])
legend(es,
       c("50%", "95%", "99%","True Outliers", "Identified Outliers"),
       col = c(gray((3:1)/(3 + 1)), "red", "green3"),
       pch = c(15, 15, 15, 16, 8), pt.cex = c(1.5, 1.5, 1.5, 1, 1), pt.lwd = 5)

# Rainbow plot order by highest density ####

hdr_order <- order(hdr_coeff$fxy)

hdr_depth <-  c()
index_hdr <- c()
rainbow_col <- rev(rainbow(110, start = 0, end = 5/6))
for(i in 1:110){
  hdr_depth <- append(hdr_depth, rep(hdr_coeff$fxy[i], endpoints[i]))
  index_hdr <- append(index_hdr, rep(which(hdr_order == i), endpoints[i]))
}

fit_all_R$Y$hdr_Depth <- hdr_depth
fit_all_R$Y$index_hdr <- index_hdr
ggplot() +
  geom_vline(xintercept = knots, alpha = .1, linetype = "dashed") +
  geom_line(data = fit_all_R$Y,
            aes(x = x,  y = fitted, group = ID, col = hdr_Depth), size = 1) +
  scale_color_gradientn(colours = rainbow(110, start = 0, end = 5/6)) +
  labs(title = "Rainbow Plot of Kernel Density", hjust = .5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5, size= rel(2)))
plot(hdr_coeff$fxy[hdr_order], pch = 16,
     col = rainbow(110, start = 0, end = 5/6),
     ylab = "Kernel Denisty Estimate")
abline(h = hdr_coeff$falpha, lty = "dashed")
text(108, hdr_coeff$falpha[1]+1, labels = "1%" )
text(108, hdr_coeff$falpha[2]+1, labels = "5%" )
text(108, hdr_coeff$falpha[3]+1, labels = "50%" )


ggplot(fit_all_R$Y, aes(x = index_hdr, y = hdr_Depth, group = ID, col = hdr_Depth)) +
  geom_point() +
  scale_color_gradientn(colours = rainbow(110, start = 0, end = 5/6)) +
  #geom_text(aes(x = index + loc, y = Depth, group = ID, label = Names), size = 4, hjust = .5) +
  geom_hline(yintercept = hdr_coeff$falpha, linetype = "dashed") +
  annotate("text", x = 100, y = hdr_coeff$falpha[1] - 1.5, label = "1%") +
  annotate("text", x = 100, y = hdr_coeff$falpha[2] + 1.5, label = "5%") +
  annotate("text", x = 100, y = hdr_coeff$falpha[3] + 1, label = "50%") +
  labs(title = "HDR Depth", hjust = .5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5, size= rel(2)))




