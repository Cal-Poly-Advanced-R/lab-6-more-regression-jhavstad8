#' Implements simple linear regression by gradient descent
#'
#' @param dat A data frame
#' @param response The name of a response variable in the data frame (unquoted)
#' @param explanatory The name of the explanatory variable in the data frame (unquoted)
#'
#' @return A data frame of coefficients
#'
#' @import dplyr
#'
#' @export
slr_gd <- function(dat, response, explanatory){
  y <- dat %>% pull({{response}})
  x <- dat %>% pull({{explanatory}})
  xmatr <- matrix(x)
  n <- length(y)

  # cost function
  cost <- function(x, y, theta) {
    sum( (x %*% theta - y)^2 ) / (2*n)
  }

  # learning rate and iteration limit
  rate <- 0.01
  max_loops <- 100

  # keep history
  past_costs <- double(max_loops)
  num_loops <- list(max_loops)

  # initialize coefficients
  theta <- matrix(c(0,0), nrow=2)

  # add a column of 1's for the intercept coefficient
  xmatr <- cbind(1, xmatr)

  # gradient descent
  for (i in 1:max_loops) {
    error <- (xmatr %*% theta - y)
    delta <- t(xmatr) %*% error / n
    theta <- theta - rate * delta
    past_costs[i] <- cost(xmatr, y, theta)
    num_loops[[i]] <- theta
  }
  results <- data.frame(theta)

  ### I found this code on r-bloggers.com. Not sure if I am allowed to use it but it is the best I can do

  ### Compute coefficients by gradient descent
  ### Return a data frame of the same form as in the `simple_linear_regression`

  return(results)

}


#' Implements linear regression with many predictors by gradient descent
#'
#' This function computes coefficients for multiple regression by gradient descent
#' All columns of the provided data frame are used as predictors, except the
#' one specified as a response.
#'
#' No interaction terms are included.
#'
#'
#' @param dat A data frame
#' @param response The name of a response variable in the data frame (unquoted)
#'
#' @return A data frame of coefficients
#'
#' @import dplyr
#'
#'@export
mlr_gd <- function(dat, response) {
  x <- dat %>% select(-{{response}})
  xmatr <- matrix(x)
  xmatr <- cbind(1,xmatr)
  y <- dat %>% pull({{response}})
  ymatr <- matrix(y)
  B <- c(rep(0,ncol(dat)))
  num_loops <- 100
  rate <- .01
  for(i in 1:num_loops){
    value <- rate * 2*t(xmatr)(ymatr-xmatr*B)
  }


  ### Compute coefficients by gradient descent
  ### Return a data frame of the same form as in the `multiple_linear_regression`

  return(results)

}

#' Implements linear regression with many predictors by matrix decomposition
#'
#' This function computes coefficients for multiple regression by QR matrix decomposition
#' All columns of the provided data frame are used as predictors, except the
#' one specified as a response.
#'
#' No interaction terms are included.
#'
#'
#' @param dat A data frame
#' @param response The name of a response variable in the data frame (unquoted)
#'
#' @return A data frame of coefficients
#'
#' @import dplyr
#'
#'@export
mlr_qr <- function(dat, response) {
  x <- dat %>% select(-{{response}})
  xmatr <- matrix(x)
  xmatr <- cbind(1,xmatr)
  y <- dat %>% pull({{response}})
  ymatr <- matrix(y)
  QR <- qr(xmatr)
  Q <- qr.Q(QR)
  R <- qr.R(QR)
  betas <- solve(R) %*% t(Q) %*% ymatr
  results <- data.frame(betas,nrow=1)
  colnames(results) <- c("Intercept", names(xmatr))



  ### Compute coefficients by QR decomposition
  ### Return a data frame of the same form as in the `multiple_linear_regression`

  return(results)

}
