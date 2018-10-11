
# > library(devtools)
# > document()


#' Create matrices for the mixed model with splines
#'
#' @description This function creates the two matrices for the fixed and
#' random components of the mixed model. The matrix corresponding to the
#' random components is a made from a spline basis.
#'
#' @param x A vector of a predictor variable.
#' @param knots A vector containing the knots for the basis.
#' @param P OPTIONAL: A value specifying the degree of the polynomial.
#' @param spline OPTIONAL: The basis to be used for the spline model.
#' This can either be "Truncated Poly" for the truncated polynomial
#' basis function or "Radial" for the radial basis function.
#' @param r OPTIONAL: A value for the degree of the radial basis function. Not
#' used in the truncated polynomial basis function.
#' @param theta OPTIONAL: A value for the range of the radial basis function. Not
#' used in the truncated polynomial basis function.
#'
#' @return \code{Z_p} is a basis matrix dimensions = n x k
#' where n equals the length of \code{x} and k equals the length of \code{int_knots}.
#' This matrix is used for the random part of the mixed model. \code{X_p} is the matrix
#' with dimensions n x P. The columns are: the intercept, x, x^2, ..., x^P. This models the mean function of
#' the data. \code{X_p} is used for the fixed part of the mixed model.
#'
#' @details The truncated polynomial basis function is of the form: (x - k)+^P.
#' The radial basis function is of the form exp(-abs(x - k)^r / theta).#'
#' Set P = 1 for linear trend.
#'
#' @seealso \code{\link{radial_basis_fun}} and \code{\link{piecewise_poly_basis_fun}}
#'
#' @references
#' \insertRef{wand2003}{unequalgroupoutlier}
#'
#' @examples
#' knots <- seq(1, 1000, len = 40)
#' x <- c(1:1000)
#' y <- rnorm(1000, x, 25)
#' model <- Mixed_Spline_Model(x, knots, P = 1)
Mixed_Spline_Model <- function(x, knots, P = 1, spline = "Radial", r = 1, theta = 2){
  # creates Xp and Zp for x
  # with K knots, for P polynomial basis
  # basis can = "Truncated poly", "Radial"

  if(spline != "Radial" && spline != "Truncated Poly"){
    stop("spline must either be Radial or Truncated Poly")
  }
  if(P < 0){stop("P cannot be negative")}

  N <- length(x)
  K <- length(knots)

  # create matrix for random coefficients
  if(spline == "Truncated Poly"){
    basis = piecewise_poly_basis_fun
    # Compute basis matrix
    Z_p <- basis(x = x, int_knots = knots, p = P)

  }else{
    basis = radial_basis_fun
    # Compute basis matrix
    Z_p <-  basis(x = x, int_knots = knots, p = P, r = r, theta = theta)
  }

  # Create the fixed effect part
  P <- c(0:P)
  X_p<- sapply(P, function(t) x^t)

  return(list("Z_p" = Z_p, "X_p" = X_p))
}

#' Fitting single dataset with mixed models and splines
#'
#' @description This function fits a mixed model.
#'
#' @param x A vector of a predictor variable.
#' @param y A vector of dependent variable
#' @param model A list containing \code{Z_p}, the matrix formed from a spline basis,
#' and \code{X_p}, the matrix formed to fit the trend. \code{Z_p} is associated
#' with the random effects and \code{X_p} is associated with the fixed effects.
#'
#' @return \code{fit} is the vector of fitted values.
#' \code{coeff_fix} is the vector of fixed effects coefficients.
#' \code{coeff_rand} is the vector of random effects coefficients.
#'
#' @details Must provide a model.
#'
#' @seealso \code{\link{Mixed_Spline_Model}} for how to create a model.

#' @references
#' \insertRef{wand2003}{unequalgroupoutlier}
#'
#' @examples
#' knots <- seq(1, 1000, len = 40)
#' x <- c(1:1000)
#' y <- rnorm(1000, x, 25)
#' model <- Mixed_Spline_Model(x, knots, P = 1)
#' fit <- Mixed_Splinefit(x, y , model)
#' plot(fit$fit, type = "l")
Mixed_Spline_Fit_Single<-function(x, y, model){
  # fits mixed model for Xp and Zp on x and y
  Zp <<- model$Z_p
  Xp <<- model$X_p

  if(is.null(Zp)){stop("Missing Z_p in model")}
  if(is.null(Xp)){stop("Missing X_p in model")}

  group.x <<- rep(1, length(x))
  group_df <- nlme::groupedData(y ~ x|group.x, data = data.frame(x, y))

  # lme(y ~ -1 + Xp,
  #     #random = nlme::pdIdent( ~ -1),
  #     data = group_df)

  mod <- nlme::lme(y ~ -1 + Xp,
                   random = nlme::pdIdent( ~ -1 + Zp),
                   data = group_df, method = "ML",
                   control = nlme::lmeControl(opt = "optim",
                                              msVerbose = F ,
                                              tolerance = 10^-3))#,
                                              #sigma = Ysd))
  #print(mod$sigma)

  # mod <- nlme::lme(y ~ -1 + Xp,
  #                  random = nlme::pdIdent( ~ -1 + Zp),
  #                  data = group_df, method = "REML",
  #                  control = nlme::lmeControl(opt = "nlminb", tolerance = 10^-3))

  #Getting Estimates
  beta_hat <- mod$coefficients$fixed
  b_hat <- unlist(mod$coefficients$random)
  if(length(beta_hat) == 1){
    fit <- Xp * beta_hat + Zp %*% b_hat
  }else{
    if(length(b_hat) == 1){
      fit <- Xp %*% beta_hat + Zp * b_hat
    }else{
      if(length(b_hat) == 1 && length(beta_hat == 1)){
        fit <- Xp * beta_hat + Zp * b_hat
      }else{
        #Find smooth function
        fit <- Xp %*% beta_hat + Zp %*% b_hat
      }
    }
  }

  #rm(Zp, Xp, group.x)
  return(list("fit" = fit, "coeff_fix" = beta_hat, "coeff_rand" = b_hat))
}

#' Fitting Mixed Models with Splines
#'
#' @description This function fits a mixed model with minimal user input.
#'
#' @param Y A tibble containing \code{ID}, \code{Y}, and \code{x}. Where \code{ID} is
#'  the distinction used between datasets, \code{Y} is the response variable and
#'  \code{x} is the predictor variable.
#' @param knots OPTIONAL: A vector of knots to be used to compute random effects in
#' mixed model. SUPPLY ONLY \code{knots} OR \code{K}, NOT BOTH.
#' @param K OPTIONAL: A value for how many number knots the user would like: DO NOT
#' SUPPLY BOTH \code{knots} AND \code{K}
#' @param P OPTIONAL: A value for the degree of the polynomial. Default P = 1 for
#' linear trend.
#' @param spline OPTIONAL: A character indicating the type of spline basis function
#' to use: "Radial" or "Truncated Poly".
#' @param r OPTIONAL: A value for the degree of the radial basis function. Not
#' used in the truncated polynomial basis function.
#' @param theta OPTIONAL: A value for the range of the radial basis function. Not
#' used in the truncated polynomial basis function.
#'
#' @return A large list is returned containing:
#' \code{Y}, a tibble containing the inputted \code{ID}, \code{y}, \code{x}, and estimated fit of the data \code{fitted}.
#' \code{coeff_fixed} a tibble containing \code{ID}, and the estimated fixed effects coefficients labeled \code{beta0},
#' \code{beta1}, ..., \code{betaP}.
#' \code{coeff_rand} a tibble containing \code{ID}, and the estimated random effects coefficients labeled \code{u1},
#' \code{u2}, ..., \code{uK}.
#'
#' @details The trend function determined by P, is a
#' polynomial: beta0 + beta1x + beta2x^2 + ... + betaPx^P. This forms the fixed effects part of the model. The bases at
#' knots create the random effects part of the mixed model. The final model looks like:
#'
#' y_i = beta0 + beta1x_i + beta2x_i^2 + ... + betaPx_i^P + sum_{j=1}^K (u_j B_j(x_i)) + epsilon_i.
#'
#' @seealso \code{\link{Mixed_Spline_Fit_Single}} to fit a mixed model with splines on individual datasets.
#'
#' @references
#' \insertRef{wand2003}{unequalgroupoutlier}
#'
#'  @examples
#' library(tidyverse)
#' x <- rep(c(1:1000), 100)
#' Y <- rnorm(100*1000, x, 2500)
#' ID <- rep(1:100, 1000)
#' ID <- ID[order(ID)]
#'
#' Y <- tibble(ID = ID, y = Y, x = x)
#' knots <- choose_knots(Y, 40)
#'
#' fit <- Mixed_Model_Spline_Fit(Y, knots, theta = 10)
#'
#' ggplot(fit$Y, aes(x = x, y = fitted, group = ID))+
#'   geom_line(alpha = .5)
Mixed_Model_Spline_Fit<-function(Y, knots, K, P = 1, spline = "Radial", r = 1, theta = 2){
  user_op <- getOption("warn")
  #options(warn = 2)

  # fits mixed model for Xp and Zp on x and y
  dimY <- dim(Y)
  unique_ID <- unique(Y$ID)
  M <- length(unique_ID)
  P <<- P
  theta <<- theta

  # Errors and missing input
  if(dimY[2] != 3){
    stop("Y must be a dataframe with a column for ID, y, and x.")
  }

  if(spline != "Radial" && spline != "Truncated Poly"){
    stop("spline must = Radial or Truncated Poly")
  }
  if(P < 0){stop("P cannot be negative")}

  if(missing(knots)){
    if(missing(K)){
      K <- ceiling(min(c(N / 4, 20)))
      knots <- choose_knots(Y, K)
    }else{
      knots <- choose_knots(Y, K = K)
    }
  }else{
    K <- length(knots)
  }

  # add fitted column to dataframe
  Y <- tibble::add_column(Y, fitted = NA)

  # create tibble for coefficients
  coeff_fix_names <- paste("beta", c(0:P), sep = "")
  coeff_random_names <- paste("u", c(1:K), sep = "")

  coeff_fix <- tibble(ID = unique_ID)
  for(i in 1:(P+1)){
    coeff_fix <- coeff_fix %>% add_column(!!(coeff_fix_names[i]) := NA)
  }
  coeff_random <- tibble(ID = unique_ID)
  for(i in 1:K){
    coeff_random <- coeff_random %>% add_column(!!(coeff_random_names[i]) := NA)
  }

  # Fit the mixed model
  new_model <- list()
  for(i in 1:M){

    df <- Y %>%
     dplyr::filter(ID == unique_ID[i])

    y <- df$Y
    x <- df$x

    # for constant functions
    if(all(y == y[1])){
      fit_all_values <- list()
      fit <- lm(y~1)
      fit_all_values$fit <- fit$fitted.values
      fit_all_values$coeff_fix <- c(fit$coefficients, 0)
      fit_all_values$coeff_rand <- rep(0, K)

    }else{

      # Determine if we should remove any knots
      last_x <- tail(x, n = 1)
      Kout <- which(knots > last_x - .001)

      model <- Mixed_Spline_Model(x = x,
                                  knots = knots,
                                  P = P,
                                  spline = spline,
                                  r = r,
                                  theta = theta)

      # t <- dim(model$Z_p)[2]
      # check <- mean(model$Z_p[, t])
      # Kout <- c()
      # while(check < .1){
      #   Kout <- append(Kout, t)
      #   t <- t - 1
      #   check <- mean(model$Z_p[, t])
      # }

      if(length(Kout) != 0){
        print(paste("Removing knots", Kout[1], "through", K, "to fit dataset", i))
        model$Z_p <- model$Z_p[, -Kout]
      }


      fit_all_values <- try(Mixed_Spline_Fit_Single(x = x, y = y, model), silent = T)
      # fit_all_values$coeff_fix
      # if(fit_all_values$coeff_fix[1] < -20){
      #   print("AHHH")
      #   break
      # }
      if(is.character(fit_all_values)){
        options(warn = user_op)
        print(fit_all_values)
        stop("Error: Try increasing theta OR decreasing r.")}
      # add 0's for coefficients of knots outside of y range
      fit_all_values$coeff_rand[Kout] <- 0
       }

    # Put into nice format for output
    Y <- Y %>%
      mutate(fitted = replace(fitted, which(ID == unique_ID[i]), fit_all_values$fit))

    coeff_fix[i, -1] <- fit_all_values$coeff_fix
    coeff_random[i, -1] <- fit_all_values$coeff_rand

  }

  # Put into nice format for output

  out_tibble <- list(Y = Y, coeff_fixed = coeff_fix, coeff_random = coeff_random, knots = knots)

  options(warn = user_op)
  remove(Zp, Xp, group.x, P, pos = ".GlobalEnv")
  return(out_tibble)
}

#################################

Mixed_Model_Spline_Fit.Keep_all_Knots<-function(Y, knots, K, P = 1, spline = "Radial", r = 1, theta = 2){
  user_op <- getOption("warn")
  options(warn = 2)

  # fits mixed model for Xp and Zp on x and y
  dimY <- dim(Y)
  unique_ID <- unique(Y$ID)
  M <- length(unique_ID)
  P <<- P
  theta <<- theta

  # Errors and missing input
  if(dimY[2] != 3){
    stop("Y must be a dataframe with a column for ID, y, and x.")
  }

  if(spline != "Radial" && spline != "Truncated Poly"){
    stop("spline must = Radial or Truncated Poly")
  }
  if(P < 0){stop("P cannot be negative")}

  if(missing(knots)){
    if(missing(K)){
      K <- ceiling(min(c(N / 4, 20)))
      knots <- choose_knots(x, Y, K)
    }else{
      knots <- choose_knots(x, Y, K = K)
    }
  }else{
    K <- length(knots)
  }

  print("This may take a few minutes")
  # add fitted column to dataframe
  Y <- tibble::add_column(Y, fitted = NA)

  # create tibble for coefficients
  coeff_fix_names <- paste("beta", c(0:P), sep = "")
  coeff_random_names <- paste("u", c(1:K), sep = "")

  coeff_fix <- tibble(ID = unique_ID)
  for(i in 1:(P+1)){
    coeff_fix <- coeff_fix %>% add_column(!!(coeff_fix_names[i]) := NA)
  }
  coeff_random <- tibble(ID = unique_ID)
  for(i in 1:K){
    coeff_random <- coeff_random %>% add_column(!!(coeff_random_names[i]) := NA)
  }

  # Fit the mixed model
  new_model <- list()
  for(i in 1:M){

    df <- Y %>%
      filter(ID == unique_ID[i])

    y <- df$y
    x <- df$x

    # for constant functions
    if(all(y == y[1])){
      fit_all_values <- list()
      fit <- lm(y~1)
      fit_all_values$fit <- fit$fitted.values
      fit_all_values$coeff_fix <- c(fit$coefficients, 0)
      fit_all_values$coeff_rand <- rep(0, K)

    }else{

      # Create the fixed effects and random effects matrices
      model <- Mixed_Spline_Model(x = x,
                                  knots = knots,
                                  P = P,
                                  spline = spline,
                                  r = r,
                                  theta = theta)
      #model$Z_p <- model$Z_p[,-c(33:40)]
      #fit23_theta41000_RM <- Mixed_Spline_Fit_Single(x = x, y = y, model)$fit
      fit_all_values <- try(Mixed_Spline_Fit_Single(x = x, y = y, model), silent = T)
      if(is.character(fit_all_values)){
        options(warn = user_op)
        print(fit_all_values)
        stop("Error: Try increasing theta OR decreasing r.")
      }
    }

    # Put into nice format for output
    Y <- Y %>%
      mutate(fitted = replace(fitted, which(ID == unique_ID[i]), fit_all_values$fit))

    coeff_fix[i, -1] <- fit_all_values$coeff_fix #[, -1] to remove ID value since it is already in Y
    coeff_random[i, -1] <- fit_all_values$coeff_rand

  }

  # Put into nice format for output

  out_tibble <- list(Y = Y, coeff_fixed = coeff_fix, coeff_random = coeff_random, knots = knots)

  options(warn = user_op) # set user options back to original
  remove(Zp, Xp, group.x, P, pos = ".GlobalEnv") # remove variables from global environment
  return(out_tibble)
}
