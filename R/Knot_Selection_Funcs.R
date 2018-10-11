

#' Choose knots over data with unequal lengths
#'
#' @param Y Matrix of functional datasets. Should contain \code{$ID}
#' indicating the separate curves, \code{$x} the time points, and \code{y} the observations.
#' @param K Total number of knots.
#'
#' @details \code{x} should be the predictor variable for the longest dataset in \code{y}.
#' Where there are observations in \code{y} the predictor variable should be the same as \code{x}.
#'
#' @return The knots based on differing endpoints.
choose_knots <- function(Y, K){

  n <- max(table(Y$ID))

  if(K < min(n/4, 20)){
    stop("K too small")
  }

  x <- Y$x

  m <- K/2
  xdiff <- max(x) - min(x)
  box <- seq(min(x), max(x) -xdiff*.001 , len = m)
  out <- length(x)
  for(i in 2:length(box)){
    out[i] <- sum(x > box[i])
  }

  myprobs <- out/sum(out)
  index <- floor(myprobs * K)
  index[which(index == 0)] <- 1

  if(sum(index) != K){
    kdiff <- sum(index) - K
    if(kdiff < 0){
      index2 <- c(1:abs(kdiff))
      index[index2] <- index[index2] + 1
    }else{
      index2 <- which(index > 1)
      index_len <- length(index2)
      if(index_len == 0){
        while(sum(index) != K){
          index <- index[-length(index)]
          m <- (max(x)-min(x))/length(index)
        }
      }else{
        a <- index_len - kdiff + 1
        index[c(a:index_len)] <- index[c(a:index_len)] - 1
      }
    }
  }

  index <- index + 1
  knots <- c()
  for(i in 2:length(box)){
    ki <- seq(box[i - 1], box[i], len = index[i - 1])[-1]
    knots <- append(knots, ki)
  }

  knots <- c(min(x), knots)
  return(knots)
}
