
endpoints <- function(X, type){
  # type = "Shift", "Iso", "Shape", "Main"
  
  n <- nrow(X)
  p <- ncol(X)
  
  if(type == "Iso"){
    tmp_diff <- abs(apply(X, 1, diff))
    jump <- apply(tmp_diff, 2, function(x) max(which(x > 6)))
    jumpmat <- cbind(ifelse(jump > .1*p, jump, .1*p), p)
    endpts <- apply(jumpmat, 1, function(x) sample(x[1]:x[2], 1))
  }else{
    endpts <- sample((.1*p):p, n, replace = T)
  }

  for(i in 1:n){
    X[i, endpts[i]:p] <- NA
  }  
  
  return(X)
  
}

mat_to_tibble <- function(y, x){
  # y are typical curves, x are anomalous curves
  
  N1 <- nrow(y)
  N2 <- nrow(x)
  N <- ncol(y)
  
  if(ncol(y) != ncol(c)){stop("y and x must have some number of columns")}
  
  y_all <- c()
  yt <- c()
  IDy <- c()
  for(i in 1:N1){
    t.na <- which(is.na(y[i, ]))
    k <- min(t.na) - 1
    if(is.infinite(k)){k <- N}
    y_all <- append(y_all, y[i, -t.na])
    yt <- append(yt, t[1:k])
    IDy <- append(IDy, rep(paste("nonA", i), k))
  }
  
  # anomalous curves
  x_all <- c()
  xt <- c()
  IDx <- c()
  for(i in 1:N2){
    t.na <- which(is.na(x[i, ]))
    k <- min(t.na) - 1
    if(is.infinite(k)){k <- N}
    if(k == 0){break}
    x_all <- append(x_all, x[i, -t.na])
    xt <- append(xt, t[1:k])
    IDx <- append(IDx, rep(paste("A", i), k))
  }
  
  # ID
  ID <- c(IDy, IDx)
  tot_curve <- c(y_all, x_all)
  tot_t <- c(yt, xt)
  # Combine all into tibble
  curve_df <- tibble("ID" = ID,
                     "Y" = tot_curve,
                     "x" = tot_t)
  
  return(curve_df)
 
}

