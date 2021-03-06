% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Mixed_Model_Funcs.R
\name{Mixed_Model_Spline_Fit}
\alias{Mixed_Model_Spline_Fit}
\title{Fitting Mixed Models with Splines}
\usage{
Mixed_Model_Spline_Fit(Y, knots, K, P = 1, spline = "Radial", r = 1,
  theta = 2)
}
\arguments{
\item{Y}{A tibble containing \code{ID}, \code{Y}, and \code{x}. Where \code{ID} is
the distinction used between datasets, \code{Y} is the response variable and
\code{x} is the predictor variable.}

\item{knots}{OPTIONAL: A vector of knots to be used to compute random effects in
mixed model. SUPPLY ONLY \code{knots} OR \code{K}, NOT BOTH.}

\item{K}{OPTIONAL: A value for how many number knots the user would like: DO NOT
SUPPLY BOTH \code{knots} AND \code{K}}

\item{P}{OPTIONAL: A value for the degree of the polynomial. Default P = 1 for
linear trend.}

\item{spline}{OPTIONAL: A character indicating the type of spline basis function
to use: "Radial" or "Truncated Poly".}

\item{r}{OPTIONAL: A value for the degree of the radial basis function. Not
used in the truncated polynomial basis function.}

\item{theta}{OPTIONAL: A value for the range of the radial basis function. Not
used in the truncated polynomial basis function.}
}
\value{
A large list is returned containing:
\code{Y}, a tibble containing the inputted \code{ID}, \code{y}, \code{x}, and estimated fit of the data \code{fitted}.
\code{coeff_fixed} a tibble containing \code{ID}, and the estimated fixed effects coefficients labeled \code{beta0},
\code{beta1}, ..., \code{betaP}.
\code{coeff_rand} a tibble containing \code{ID}, and the estimated random effects coefficients labeled \code{u1},
\code{u2}, ..., \code{uK}.
}
\description{
This function fits a mixed model with minimal user input.
}
\details{
The trend function determined by P, is a
polynomial: beta0 + beta1x + beta2x^2 + ... + betaPx^P. This forms the fixed effects part of the model. The bases at
knots create the random effects part of the mixed model. The final model looks like:

y_i = beta0 + beta1x_i + beta2x_i^2 + ... + betaPx_i^P + sum_{j=1}^K (u_j B_j(x_i)) + epsilon_i.
}
\references{
\insertRef{wand2003}{unequalgroupoutlier}

 @examples
library(tidyverse)
x <- rep(c(1:1000), 100)
Y <- rnorm(100*1000, x, 2500)
ID <- rep(1:100, 1000)
ID <- ID[order(ID)]

Y <- tibble(ID = ID, y = Y, x = x)
knots <- choose_knots(Y, 40)

fit <- Mixed_Model_Spline_Fit(Y, knots, theta = 10)

ggplot(fit$Y, aes(x = x, y = fitted, group = ID))+
  geom_line(alpha = .5)
}
\seealso{
\code{\link{Mixed_Spline_Fit_Single}} to fit a mixed model with splines on individual datasets.
}
