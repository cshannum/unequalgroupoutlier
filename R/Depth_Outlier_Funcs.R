
# This script contains 2 functions to boostrap a cutoff value C for the depth plots
# methods from paper: MultivariateFuncDepth_NOx pgs 8-9

#' Bootstrap a cutoff value to identify anomalies
#'
#' @param coeff A dataframe of coefficients of interest. The first column is \code{ID} identifier.
#'  The rest of the columns are for the parameter to be estimate. Each row is the estimated parameters
#'  fore each curve.
#' @param d.method A character string determining the depth function to use: "LP", "Projection",
#' "Mahalanobis", or "Euclidean". It is suggested to not use "Tukey" due to singularity in
#' coefficient matrix. For details see \code{\link[DepthProc]{depth}}
#' @param c.method A character string determining the method to estimate the cutoff value.
#' This can be "depth" or "alpha".
#' @param alpha A value determining the percentage of rows to remove from \code{coeff}.
#' \code{alpha} should be between (0, 1) with a suggested value of 0.05. Do not need to
#' identify if \code{c.method} = "depth".
#' @param B A value determining how many bootstrap datasets should be made to estimate
#' the cutoff value with a suggested rate of 1000.
#'
#' @details The function starts by computing the depths for each parameter set by \code{d.method}.
#'
#' The "alpha" \code{c.method} removes the alpha percent least deep coefficients. The rest of the
#' coefficients are bootstrapped and new depths are computed for each new bootstrapped set. The
#' 1% empirical percentile of the depths is saved. The cuttoff value is the median of these
#' 1% empirical percentile of the depths.
#'
#' The "depth" \code{c.method} bootstraps the coefficients with probability related to the
#' original depth values.  New depths are computed for each new bootstrapped set. The
#' 1% empirical percentile of the depths is saved. The cuttoff value is the median of these
#' 1% empirical percentile of the depths.
#'
#' @seealso \code{\link[DepthProc]{depth}}, \code{\link{bootstrap_C.alpha}}, and
#'  \code{\link{bootstrap_C.depth}}
#'
#' @return \code{$d} the depths computed by \code{d.method} over all coefficients.
#' \code{$Cb} the cutoff value; depths below cutoff may be anomalous.
bootstrap_C <- function(coeff, d.method, c.method, alpha, B){
  # B is the number of times to bootstrap (paper uses B = 200)
  if(missing(coeff)){stop("Coefficients not supplied")}

  d <- DepthProc::depth(coeff, coeff, method = d.method)
  d <- d@.Data

  if(c.method == "alpha"){
    if(missing(alpha)){stop("Must supply alpha value")}
    Cb <- bootstrap_C.alpha(coeff = coeff, d = d, B = B, alpha = alpha, d.method = d.method)
  }else{
    Cb <- bootstrap_C.depth(coeff = coeff, d = d, B = B, d.method = d.method)
  }

  d.list <- list("depths" = d, "Cb" = Cb)
  return(d.list)

}

#' Alpha method to bootstrap a cutoff value to identify anomalies
#'
#' @param coeff A dataframe of coefficients of interest. Each column is for the parameter to be estimate. Each row is the estimated parameters
#'  fore each curve. (No ID column here)
#'@param B A value determining how many bootstrap datasets should be made to estimate
#' the cutoff value with a suggested rate of 1000.
#' @param d A vector of depths for all of the coefficients.
#' @param alpha A value determining the percentage of rows to remove from \code{coeff}.
#' \code{alpha} should be between (0, 1) with a suggested value of 0.05. Do not need to
#' identify if \code{c.method} = "depth".
#' @param d.method A character string determining the depth function to use: "LP", "Projection",
#' "Mahalanobis", or "Euclidean". It is suggested to not use "Tukey" due to singularity in
#' coefficient matrix. For details see \code{\link[DepthProc]{depth}}
#'
#' @details The "alpha" \code{c.method} removes the alpha percent least deep coefficients. The rest of the
#' coefficients are bootstrapped and new depths are computed for each new bootstrapped set. The
#' 1% empirical percentile of the depths is saved. The cuttoff value is the median of these
#' 1% empirical percentile of the depths.
#'
#' @seealso \code{\link[DepthProc]{depth}}, \code{\link{bootstrap_C.alpha}}, and
#'  \code{\link{bootstrap_C.depth}}
#'
#' @return \code{$Cb} the cutoff value
bootstrap_C.alpha <- function(coeff, d, B, alpha, d.method){
  n <- length(d) # number of curves

  # remove the alpha % least deep curves
  t <- max(floor(n*.01), 1)
  rm.index <- order(d, decreasing = F)[1:t]

  # bootstrap remaining depths
  B.index <- sample(c(1:n)[-rm.index], n*B, replace = T)
  # put into easy format
  B.index <- matrix(B.index, nr = n, nc = B)

  cutoff <- rep(0, B)


  for(i in 1:B){
    if(class(coeff) == "matrix"){
      B.coeff <- coeff[B.index[, i], ]
    }else{
      B.coeff <- coeff[B.index[, i]]
    }
    # compute new depths for new dataset
    d.new <-  DepthProc::depth(B.coeff, method = d.method)@.Data
    # compute empirical 1% of each bootstrapped depth set
    cutoff[i] <- sort(d.new)[t]
  }
  # B.coeff <- apply(B.index, 2, function(x) list(coeff[x, ]))
  # # compute new depths for new dataset
  # d.new <- lapply(B.coeff, function(x) DepthProc::depth(x[[1]], x[[1]], method = d.method)@.Data)
  # compute empirical 1% of each bootstrapped depth set
  #cutoff <- unlist(lapply(d.new, function(x) x[order(x, decreasing = T)[floor(n*.01)]]))

  # final cutoff is the median of all bootstrapped cutoffs
  Cb <- median(cutoff)

  return(Cb)

}

#' Bootstrap a cutoff value to identify anomalies
#'
#' @param coeff A dataframe of coefficients of interest. Each column is for the parameter to be estimate. Each row is the estimated parameters
#'  fore each curve. (No ID column here)
#' @param B A value determining how many bootstrap datasets should be made to estimate
#' the cutoff value with a suggested rate of 1000.
#' @param d A vector of depths for all of the coefficients.
#' @param d.method A character string determining the depth function to use: "LP", "Projection",
#' "Mahalanobis", or "Euclidean". It is suggested to not use "Tukey" due to singularity in
#' coefficient matrix. For details see \code{\link[DepthProc]{depth}}
#'
#' @details The "depth" \code{c.method} bootstraps the coefficients with probability related to the
#' original depth values.  New depths are computed for each new bootstrapped set. The
#' 1% empirical percentile of the depths is saved. The cuttoff value is the median of these
#' 1% empirical percentile of the depths.
#'
#' @seealso \code{\link[DepthProc]{depth}}, \code{\link{bootstrap_C.alpha}}, and
#'  \code{\link{bootstrap_C.depth}}
#'
#' @return \code{$Cb} the cutoff value.
bootstrap_C.depth <- function(coeff, d, B, d.method){
  n <- length(d)

  # compute probability for each depth
  d.prob <- d/sum(d)
  # bootstrap depths with probability related to depth
  B.index <- sample(c(1:n), n*B, prob = d.prob, replace = T)
  # put into easy format
  B.index <- matrix(B.index, nr = n, nc = B)
  cutoff <- rep(0, B)
  t <- max(floor(n*.01), 1)
  for(i in 1:B){
    if(class(coeff) == "matrix"){
      B.coeff <- coeff[B.index[, i], ]
    }else{
      B.coeff <- coeff[B.index[, i]]
    }
    # compute new depths for new dataset
    d.new <-  DepthProc::depth(B.coeff, B.coeff, method = d.method)@.Data
    # compute empirical 1% of each bootstrapped depth set
    cutoff[i] <- sort(d.new)[t]
  }
  # final cutoff is the median of all bootstrapped cutoffs
  Cb <- median(cutoff)
  return(Cb)
}

#' Identify anomalous coefficients using depth
#'
#' @param coeff A dataframe of coefficients of interest. The first column is \code{ID} identifier.
#'  The rest of the columns are for the parameter to be estimate. Each row is the estimated parameters
#'  fore each curve.
#' @param d.method A character string determining the depth function to use: "LP", "Projection",
#' "Mahalanobis", or "Euclidean". It is suggested to not use "Tukey" due to singularity in
#' coefficient matrix. For details see \code{\link[DepthProc]{depth}}
#' @param c.method A character string determining the method to estimate the cutoff value.
#' This can be "depth" or "alpha".
#' @param alpha A value determining the percentage of rows to remove from \code{coeff}.
#' \code{alpha} should be between (0, 1) with a suggested value of 0.05. Do not need to
#' identify if \code{c.method} = "depth".
#' @param B A value determining how many bootstrap datasets should be made to estimate
#' the cutoff value with a suggested rate of 1000.
#'
#' @details The function uses a bootstrap method to estimate a cutoff depth value.
#' Depths below this cutoff depth value are flagged as anomalous.
#'
#' The "alpha" \code{c.method} removes the alpha percent least deep coefficients. The rest of the
#' coefficients are bootstrapped and new depths are computed for each new bootstrapped set. The
#' 1% empirical percentile of the depths is saved. The cuttoff value is the median of these
#' 1% empirical percentile of the depths.
#'
#' The "depth" \code{c.method} bootstraps the coefficients with probability related to the
#' original depth values.  New depths are computed for each new bootstrapped set. The
#' 1 percent empirical percentile of the depths is saved. The cuttoff value is the median of these
#' 1 percent empirical percentile of the depths.
#'
#' @references
#' \insertRef{febrero2008}{unequalgroupoutlier}
#'
#' @seealso \code{\link[DepthProc]{depth}}, \code{\link{bootstrap_C.alpha}}, and
#'  \code{\link{bootstrap_C.depth}}
#'
#' @return A list with \code{$outliers} is a tibble containing the ID, outlier identification, depth value.
#' True identifies outliers. \code{$Cb} is the cutoff value
depth_Outliers <- function(coeff, d.method = "L2", c.method = "depth", alpha = .05, B = 1000){

  if(class(coeff)[1] == "matrix"){
    d.ID <- c(1:nrow(coeff))

  }else{
    if(class(coeff)[1] == "numeric"){
      d.ID <- c(1:length(coeff))
    }else{
      index <- which(colnames(coeff) == "ID")

      d.ID <- unlist(coeff[, index])
      coeff <- coeff[, -index]
    }
  }

  d.list <- bootstrap_C(coeff = coeff, d.method = d.method, c.method = c.method, alpha = alpha, B = B)
  d.out <- d.list$depths <= d.list$Cb
  d.out <- dplyr::tibble("ID" = d.ID,
                         "outliers" = d.out,
                         "depth"  = d.list$depths)
  outlist <- list("outliers" = d.out,
                  "Cb" = d.list$Cb)

  return(outlist)
  # returns vector T = outlier; F = not outlier
}
