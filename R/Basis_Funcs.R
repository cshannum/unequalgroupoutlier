
#' Truncated Piecewise polynomial basis
#'
#' @param x A vector of a predictor variable.
#' @param int_knots A vector containing the interior knots.
#' @param p A value specifying the degree of the polynomial.
#' @return A matrix consisting of (x_i - k_j)^p if x_i > k_j or 0 otherwise for
#' i = 1,...,n and j = 1,... k. Where n equals the length of \code{x} and k equals the length of \code{int_knots}.
#' @examples
#' out <- piecwise_poly_basis_fun(x = c(1:9), k = c(2:4), p = 3)
#' plot(out[,1], type = "l") # see what your basis looks like
#'
#' @references
#' \insertRef{wand2003}{unequalgroupoutlier}
piecewise_poly_basis_fun <- function(x, int_knots, p) {
  #creates piecwise polynomial, p, basis for x with knot
  # x vector
  # int_knots vector
  # p value

  if(P < 0){stop("P cannot be negative")}

  # put into matrix form
  x <- matrix(x, nrow = length(x), ncol = length(int_knots))
  int_knots <- matrix(int_knots, nrow = dim(x)[1], ncol = length(int_knots), byrow = T)

  Z <- (x > int_knots)*(x - int_knots)^p
  return(Z)

}

#' Radial basis
#'
#' @param x A vector of a predictor variable.
#' @param int_knots A vector containing the interior knots.
#' @param p A value specifying the degree of the polynomial.
#' @param theta OPTIONAL: A value specifying the range parameter.
#' @return \code{Z} A radial basis matrix dimensions = n x k where n equals the length of \code{x} and k equals the length of \code{int_knots}.
#' @details Let the {i, j} value of \code{Z_p} be exp( - abs(x_i - k_j)^p / theta) and the {i, j}
#'  value of \code{Omega} be exp( - abs(k_i - k_j)^p / theta), then \code{Z} = \code{Z_p}%%*%%\code{Omega}^-1/2. The correction
#'  is required to be able to use this matrix in the mixed model framework.
#'
#' @references
#' \insertRef{wand2003}{unequalgroupoutlier}
#' @examples radial_basis_fun(x = c(1:9), k = c(2:4), p = 3)
radial_basis_fun <- function(x, int_knots, p, r = 1, theta = 2){
  # Creates Z matrix for radial basis
  if(p < 0){stop("P cannot be negative")}

  # cov(||K_i - K_j||) matrix
  Ok <- matrix(int_knots, nrow = length(int_knots), ncol = length(int_knots))

  # put into matrix form
  x <- matrix(x, nrow = length(x), ncol = length(int_knots))
  int_knots <- matrix(int_knots, nrow = dim(x)[1], ncol = length(int_knots), byrow = T)

  # cov(||x - K_i||) matrix
  Z_k <- exp(-abs(x - int_knots)^r / theta)

  # matrix to adjust Z_k to be used in mixed model
  Omega_k <- exp(-abs(Ok - t(Ok))^r / theta)

  # compute svd to take squareroot of matrix
  svd_O <- svd(Omega_k)
  # need Omega ^ -1/2
  neg_sqrt_Omega_k <- solve(svd_O$u %*% diag(sqrt(svd_O$d)) %*% t(svd_O$v))

  Z <- Z_k %*% neg_sqrt_Omega_k
}

# B_spline_basis_fun <- function(x, knots, p, drv = 1){
#   # Creates Z matrix for B-splines
#
#   range.x <- c(1.05 * min(x) - 0.05 * max(x), 1.05 * max(x) - 0.05 * min(x))
#
#   allKnots <- c(rep(range.x[1], 4), knots, rep(range.x[2], 4))
#
#   B <- spline.des(allKnots, x, ord= P + 1, derivs = drv,
#                    outer.ok = TRUE)$design
#   return(B)
# }
