
rm(list = ls())
n <- 1000
knots <- seq(1, n, len = 40)
x <- c(1:n)

Y <-rmvnorm(100, x, x^3*diag(n)) + rmvnorm(100, x, 300*diag(n))
plot(x, Y[1,])

endpoints <- sample(300:n, n, replace = T)
for(i in 1:100){
  Y[i, c(endpoints[i]:n)] <- NA
}

fit <- Mixed_Spline_Fit_Full(Y, x, spline = "Radial", P = 3)

par(mfrow = c(2, 2))
plot(x, Y[1,])
lines(x, fit$fit[1, ], lwd = 4, col = "red")
abline(v = fit$knots, lty = 2)
plot(x, Y[2,])
lines(x, fit$fit[2, ], lwd = 4, col = "red")
abline(v = fit$knots, lty = 2)
plot(x, Y[3,])
lines(x, fit$fit[3, ], lwd = 4, col = "red")
abline(v = fit$knots, lty = 2)
plot(x, Y[4,])
lines(x, fit$fit[4, ], lwd = 4, col = "red")
abline(v = fit$knots, lty = 2)

