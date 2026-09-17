################################################################################
#
# OPTIMIZATION ALGORITHMS OF ZHANG ET AL. (2021). THE CODE IS EXTRACTED FROM
# https://github.com/qiongzhangclemson/optimaldesign
# 
# THE CODE HAS BEEN MODIFIED TO ACCOUNT FOR UNBALANCED CLINICAL TRIALS
# MODIFICATIONS BY ALAN R. VAZQUEZ (TECNOLOGICO DE MONTERREY).
#
# Zhang, Q., Khademi, A., and Song, Y. (2021). Mini-max optimal design of 
# two-armed trials with side information. INFORMS Journal of Computing, 
# 34:165–182
#
################################################################################

opt.design.approx <- function(H, balance=FALSE, e.cnst = 6, max.time = 300, print.output = 1) {
    #########################################
    # Function to generate LB_APPROX design
    #########################################  
    R <- H%*%solve(t(H)%*%H)%*%t(H)
    Q <- R*R
    x <- matrix(mip(Q, balance, e.cnst, max.time,print.output), ncol=1)
    return(x)
}

opt.design.iter <- function(H, balance=FALSE, e.cnst = 6, MaxIT=100, tol=0.0001, max.time = 500, print.output=1) {
  #########################################
  # Function to generate EXACT designs
  ######################################### 
  # data related parameters
  n <- nrow(H)
  R <- solve(t(H)%*%H)
  R0 <- R %*% t(H)
  R1 <- H%*%R0
  Zs <- H
  Zs <- unique(Zs)
  
  
  # initial setting
  IT <- 0
  diff <- 10
  checkre <-NULL
  while(diff>tol & IT<=MaxIT) {
    results <- mip2(R, R0, R1, Zs, balance=balance, e.cnst, max.time)
    theta.hat <- results[1]
    xast <- results[-1]
    Rast <- R+R0%*%diag(xast)%*%R1%*%diag(xast)%*%t(R0)
    zast <- lp_sdp(Rast,max.time,print.output) 
    theta.bar <- t(zast)%*%Rast%*%zast
    Zs <- rbind(Zs, zast)
    diff <- theta.bar-theta.hat
    IT <- IT+1
    checkre <- rbind(checkre, c(IT, theta.hat, theta.bar))
  }
  return(matrix(xast, ncol=1))
}


lp_sdp <- function(R, max.time, print.output=1) {
  ##########################
  # Input: 
  # lp solve Semidefinite matrix R
  # Output:
  # x
  ##########################
  model <- list()
  params <- list(MIPGap=0.001, TimeLimit=max.time,OutputFlag=print.output)
  
  n <- nrow(R)
  R22 <- R[2:n, 2:n]
  R21 <- R[1, 2:n]
  model$obj   <- c(as.vector(R22), R21-rowSums(R22))
  model$modelsense <- 'max'
  
  A0 <- diag(1, (n-1)^2)
  A1 <- rbind(-A0, -A0, A0, A0)
  A22 <- matrix(0, (n-1)^2, n-1)
  A23 <- -kronecker(diag(1, n-1), matrix(1, n-1, 1))
  A24 <- -kronecker(matrix(1, n-1, 1), diag(1, n-1))
  A21 <- -A23-A24
  A2 <- rbind(A21, A22, A23, A24)
  model$A <- cbind(A1, A2)
  model$rhs   <- c(rep(1, (n-1)^2), rep(0, 3*(n-1)^2))
  model$sense <- rep('<', nrow(model$A))
  model$vtype <- c(rep('C',(n-1)^2), rep('B', n-1))
  result <-gurobi(model, params)
  # The elements in the output vector z can only take 
  # values -1 and +1
  z <- 2*result$x[-(1:(n-1)^2)]-1
  z <- c(1, z)
  return(z)
}

mip <- function(Q, balance = FALSE, e.cnst, max.time, print.output = 1) {
  ##########################
  # the MIP problem in approx optimal design
  ##########################
  n <- nrow(Q)
  model <- list()
  params <- list(MIPGap=0.02, TimeLimit=max.time,OutputFlag=print.output)
  
  if(balance){
    model$A     <- matrix(1, 1, n)
    model$rhs <- n/2  
    model$sense <- c('=')
  } else {
    model$A     <- matrix(1, 2, n)
    model$rhs <- c(n/2 + e.cnst, n/2 - e.cnst)  
    model$sense <- c('<=', '>=')
  }
  
  model$Q     <- Q
  model$obj   <- -rowSums(Q)
  model$vtype <- rep('B', n)
  result <-gurobi(model, params)
  return(2*result$x-1)
}


mip2 <- function(R, R0, R1, Zs, balance = FALSE, e.cnst, max.time, print.output = 1) {
  ##########################
  # the MIP problem in Exact optimal design
  ##########################
  n <- nrow(R1)
  model <- list()
  params <- list(TimeLimit=max.time, MIPGap=0.001,OutputFlag=print.output)
  model$modelsense <- 'min'
  
  if(balance){
    model$A     <- matrix(1, 1, n+1)
    model$A[1,1]     <- 0    
    model$rhs   <- n/2
    model$sense <- c('=')
  } else {
    model$A     <- matrix(1, 2, n+1)
    model$A[1,1]     <- 0
    model$rhs <- c(n/2 + e.cnst, n/2 - e.cnst)  
    model$sense <- c('<=', '>=')
  }
  
  model$obj   <- c(1, rep(0, n))
  model$vtype <- c('C',rep('B', n))
  model$quadcon <- list()
  for(i in 1:nrow(Zs)) {
    z <- Zs[i,]
    qc <- list()
    qc$Qc <- matrix(0, n+1, n+1)
    hatH <- R1*(t(R0)%*% outer(z, z) %*% R0)
    qc$Qc[2:(n+1), 2:(n+1)] <- 4*hatH
    qc$q <- c(-1, -4*rowSums(hatH))
    qc$rhs <- -t(z)%*%R%*%z-sum(hatH)
    model$quadcon[[i]] <- qc
  }
  result <-gurobi(model, params)
  return(c(result$x[1], 2*result$x[-1]-1))
}

