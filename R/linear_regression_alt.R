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

  # set values and coefs
  y <- dat %>% pull({{response}})
  x <- dat %>% pull({{explanatory}})
  xmatr <- matrix(x)
  n <- length(y)
  m = 0
  c = 0
  # learning rate and iteration limit
  rate <- 0.01
  num_loops = 0

  # gradient descent
  repeat{
  Y_pred = m*x + c  # The current predicted value of Y
  D_m = (-2/n) * sum(x * (y - Y_pred))  # Derivative wrt m
  D_c = (-2/n) * sum(y - Y_pred)  # Derivative wrt c
  m = m - rate * D_m  # Update m
  c = c - rate * D_c  # Update c
  num_loops = num_loops + 1
  if(abs(D_m) < .001 ){
    if(abs(D_c) < .001){
      break
    }
  }
  }
  results <- matrix(c(c,m),nrow = 1)
  colnames(results) <- c("Intercept", names(dat %>% select({{explanatory}})))



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
  xmatr <- as.matrix(x)
  xmatr <- cbind(1,xmatr)
  y <- dat %>% pull({{response}})
  ymatr <- as.matrix(y)
  B <- matrix(0,ncol(dat))
  num_loops <- 100000
  rate <- .01
  n <- nrow(dat)
  for(i in 1:num_loops){
    cost <- sum((xmatr %*% B - ymatr)^2) /(2*n)
    B <- B - rate * (1/n)*(t(xmatr) %*% (xmatr%*%B - ymatr))
    }

    results <- as.data.frame(t(B))
    names(results)[1] <- "Intercept"


    ### I took some of this code from an RPubs article by Mike Fang and adapted it

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
  xmatr <- as.matrix(x)
  xmatr <- cbind(1,xmatr)
  y <- dat %>% pull({{response}})
  ymatr <- as.matrix(y)
  QR <- qr(xmatr)
  Q <- qr.Q(QR)
  R <- qr.R(QR)
  betas <- solve(R) %*% t(Q) %*% ymatr
  results <- as.data.frame(t(betas),nrow=1)
  colnames(results) <- c("Intercept", names(x))



  ### Compute coefficients by QR decomposition
  ### Return a data frame of the same form as in the `multiple_linear_regression`

  return(results)

}
