# This script creates unequal toy data

#' Create Unequal functional data
#'
#' @param M Number of non-anomalous functional datasets
#' @param A Total number of anomalous datasets
#' @param N Number of total possible observations in each dataset.
#'
#' @return A list containing the toy data \code{data} each row is a dataset.
#' \code{label} for if the dataset is normal ("N") or anomalous ("A1" or "A2")
#' Three functions used to create the functional data: \code{f_x}, the function used to create the normal data,
#' \code{f_x.A1} the function used to create one set of anomalous data, and \code{f_x.A2} the
#' function used to create the second set of anomalous data.
Unequal_ToyData <- function(M = 100, A = 10, N = 1000){


  if(N <= 500){stop("N must be greater than 500")}
  #Start by creating data of same length
  sigma_eps <- .1
  x <- seq(0, 10, len = N)
  y <- matrix(NA, nrow = M, ncol = N)
  A1 <- ceiling(A/2)
  A2 <- floor(A/2)
  label <- c(rep("N", M), rep("A1", A1), rep("A2", A2))

  for(i in 1:M){
    f_x1 <- dnorm(x, mean = 2, sd = .25) +  2*dnorm(x, mean = 3, sd = 2) -
      dnorm(x, mean = 7, sd = .75)
    y[i,] <- f_x1 + rnorm(N, sd = sigma_eps)
  }

  # create outlier data
  y.out <- matrix(NA, nrow = A1, ncol = N)
  for(i in 1:A1){
    f_x2 <- 2*dnorm(x, mean = 3, sd = 2) -
      dnorm(x, mean = 7, sd = .75)
    y.out[i,] <- f_x2 + rnorm(N, sd = sigma_eps)
  }

  y.out2 <- matrix(NA, nrow = A2, ncol = N)
  for(i in 1:A2){
    f_x3 <- 2*dnorm(x, mean = 5, sd = 1) + 2*dnorm(x, mean = 3, sd = 2)
    y.out2[i,] <- f_x3 + rnorm(N, sd = sigma_eps)
  }

  # combine all data into one matrix
  data <- rbind(y, y.out, y.out2)

  # Now cut off some data
  endpoints <- sample((N-500):N, dim(data)[1], replace = T)
  for(i in 1:dim(data)[1]){
    data[i, ] <- c(data[i ,c(1:endpoints[i])], rep(NA, N - endpoints[i]))
  }

  mylist <- list("data" = data,
                 "x" = x,
                 "endpoints" = endpoints,
                 "labal" = label,
                 "f_x" = f_x1,
                 "f_x.A1" = f_x2,
                 "f_x.A2" = f_x3)

  return(mylist)
}













