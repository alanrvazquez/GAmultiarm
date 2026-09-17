################################################################################
#
# AUXILIAR FUNCTIONS TO CONSTRUCT AND EVALUATE CLINICAL TRIALS. 
# 
#
# Author: Alan Vazquez
# Affiliation: Tecnologico de Monterrey
# Email: alanrvazquez@tec.mx
#
################################################################################

moments.matrix <- function(m, n.bin = 0, type = 'l'){
  #=============================================================================
  # moments.matrix: Compute the moments matrix for linear and quadratic models.
  #                 For linear models, the function can output the moments matrix
  #                 for a model containing continuous and binary covariates in 
  #                 {-1, +1}.
  # Inputs: 
  #       m: number of covariates.
  #       n.bin: number of binary covariates.
  #       type: 'l' for linear model, 'q' for linear+quadratic model, and 
  #              'full' for linear, quadratic and interaction model.
  #
  # Output:
  #       M: Moments matrix.
  #=============================================================================
  if (type == 'l'){
    M <- 1/3*diag(m+1)
    M[1,1] <- 1
    
    # If there are binary variables, assign them at the end.
    if (n.bin > 0){
      las.elem <- tail(1:(m+1),n.bin)
      diag(M)[las.elem] <- 1
    }
    return(M)
  } else if (type == 'full' && n.bin == 0){
    p <- (m+1)*(m+2)/2 
    M <- diag(p)
    diag(M) <- c(1, rep(1/3, m), rep(1/9, m*(m-1)/2), rep(1/5-1/9, m))
    las.elem <- tail(1:p,m)
    M[las.elem,las.elem] <- M[las.elem,las.elem] + 1/9
    M[1,las.elem] <- 1/3
    M[las.elem,1] <- 1/3
    return(M)
  } else if (type == 'q' && n.bin == 0){
    p <- 1 + 2*m 
    M <- diag(p)
    diag(M) <- c(1, rep(1/3, m), rep(1/5-1/9, m))
    las.elem <- tail(1:p,m)
    M[las.elem,las.elem] <- M[las.elem,las.elem] + 1/9
    M[1,las.elem] <- 1/3
    M[las.elem,1] <- 1/3
    return(M)
  } else {print('Undefined model')}
  
}

Ioptimality <- function(H, X, m, n.bin = 0, model.type = 'l', approx = FALSE) {
  #=============================================================================
  # Calculate I-optimality criterion for several two-arm trials in X.
  #
  # INPUTS:
  # H:      n x p covariate matrix.
  # X:      n x ndesigns matrix with two-arm designs.
  # n.bin:  number of discrete variables.
  # approx: Should we use the inverse of its approximation.
  #
  # OUTPUT:
  # ave.obj: 1 x ndesign vector with I-optimality values for the ndesigns in X.
  #
  #=============================================================================
  
  ss <- ncol(X)
  ave.obj <- rep(999, ss)
  M <- moments.matrix(m, n.bin, model.type)
  
  for (j in 1:ss) {
    x <- X[,j]
    n <- length(x)
    D <- diag(x)
    R <- t(H)%*%H
    Rinv <- solve(R)
    Rx <- t(H)%*%D%*%H
    Sigma <- R-Rx%*%Rinv%*%Rx
    
    if ( approx ){
      Sigma.approx <- Rinv+Rinv%*%Rx%*%Rinv%*%Rx%*%Rinv
      ave.obj[j] <- sum(diag(Sigma.approx %*% M))
    } else {
      # Compute the original covariance matrix in Equation (7) of Zhang et al. (2021)
      Sigmainv <- try(solve(Sigma), TRUE) 
      # If matrix is invertible, then use it to find z which maximizes the variance
      # over the whole feasible covariate space. 
      if(is.matrix(Sigmainv)) {
        ave.obj[j] <- sum(diag(Sigmainv %*% M))
      }
    }
    
  }
  return(ave.obj)
}

max_var <- function(R, n.cat = 0, max.time, print.output=1) {
  ##########################
  # Find the vector with the maximum SC variance with continuous 
  # and binary variables 
  # Input: 
  # Semidefinite matrix R
  # NOTE: The model must have an intercept, p-n.cat-1 continuous covaraites
  #       and n.cat two-level binary variables, in that order.
  # Output:
  # z
  ##########################
  
  # If all variables are categorical, then use Zhang's algorithm.
  if(n.cat == (ncol(R)-1)){
    z <- lp_sdp(R, max.time, print.output)
  } else { 
    # Otherwise, use our algorithm.
    model <- list()
    params <- list(MIPGap = 0.001, TimeLimit = max.time,
                   OutputFlag = print.output)
    
    p <- ncol(R)
    # Variance-covariate matrix for covariates.
    R22 <- R[2:p, 2:p]
    # Cavriances and covariances for the intercept.
    R21 <- R[1, 2:p]
    
    model$Q   <- R22
    model$obj   <- R21 - rep(1,p-1)%*%R22
    model$modelsense <- 'max'
    if (p - n.cat <= 1){
      model$vtype <- c(rep('C', p-1))  
    } else {
      model$vtype <- c(rep('C', p-n.cat-1), rep('B', n.cat))  
    }
    model$A <- diag(1, p-1)
    model$rhs <- rep(1, p-1)
    model$sense <- rep("<=", p-1)
    
    result <- gurobi(model, params)
    z <- 2*(result$x) - 1
    z <- c(1, z)
    
  }
  
  return(z)
}


est_lp <- function(H, X, Zs, n.cat = 0, max.time = 300, print.output = 0) {
  
  #=============================================================================
  # Calculate G-optimality criterion.
  # If the number of covariates is smaller than 50,
  # the function outputs the maximum variance over the
  # whole feasible space of covariates.
  # Otherwise, it outputs the maximum variance 
  # over the covariate vectors in matrix Zs.
  #=============================================================================
  
  p <- ncol(H)
  ss <- ncol(X)
  max.obj <- rep(999, ss)
  max.obj.approx <- rep(999, ss)
  if(p <= 50) {
    for(j in 1:ss) {
      x <- X[,j]
      n <- length(x)
      D <- diag(x)
      R <- t(H)%*%H
      Rinv <- solve(R)
      Rx <- t(H)%*%D%*%H
      Sigma <- R-Rx%*%Rinv%*%Rx
      # Compute the original covariance matrix in Equation (7) of Zhang et al. (2021)
      Sigmainv <- try(solve(Sigma), TRUE) 
      # If matrix is invertible, then use it to find z which maximizes the variance
      # over the whole feasible covariate space. 
      if(is.matrix(Sigmainv)) {
        zast <- max_var(Sigmainv, n.cat, max.time, print.output)
        max.obj[j] <- t(zast)%*% Sigmainv %*% zast
      } 
      Sigma.approx <- Rinv+Rinv%*%Rx%*%Rinv%*%Rx%*%Rinv
      zast <- max_var(Sigma.approx, n.cat, max.time, print.output)
      max.obj.approx[j] <- t(zast)%*% Sigma.approx %*% zast
    }
    re <- data.frame(obj=max.obj, obj.approx=max.obj.approx)
  }
  if(p > 50) {
    Sigmainvs <- list()
    Sigma.approx <- list()
    for(j in 1:ncol(X)) {
      x <- X[,j]
      n <- length(x)
      D <- diag(x)
      R <- t(H)%*%H
      Rinv <- solve(R)
      Rx <- t(H)%*%D%*%H
      Sigma <- R-Rx%*%Rinv%*%Rx
      Sigmainvs[[j]] <- solve(Sigma)
      Sigma.approx[[j]] <- Rinv+Rinv%*%Rx%*%Rinv%*%Rx%*%Rinv
    }
    zs <- matrix(0, 0, ncol(X))
    for(i in 1:nrow(Zs)) {
      z <- Zs[i, ]
      zs <- rbind(zs, sapply(Sigmainvs, function(U) t(z)%*%U%*%z))
    }
    max.obj <- apply(zs, 2, max)
    zs <- matrix(0, 0, ncol(X))
    for(i in 1:nrow(Zs)) {
      z <- Zs[i, ]
      zs <- rbind(zs, sapply(Sigma.approx, function(U) t(z)%*%U%*%z))
    }
    max.obj.approx <- apply(zs, 2, max)
    re <- data.frame(Gopt=max.obj, Gopt.approx=max.obj.approx)
  }
  return(re)
}

sing_Iopt <- function(x, p, n, ss, H, R, Rinv, M, approx = TRUE){
  #=============================================================================
  # Calculate I-optimality criterion for several two-arm trials in X.
  #
  # INPUTS:
  # H:      n x p covariate matrix.
  # X:      n x ndesigns matrix with two-arm designs.
  # n.bin:  number of discrete variables.
  # approx: Should we use the inverse of its approximation.
  #
  # OUTPUT:
  # ave.obj: 1 x ndesign vector with I-optimality values for the ndesigns in X.
  #
  #=============================================================================
  
  D <- diag(x)
  Rx <- t(H)%*%D%*%H
  Sigma <- R-Rx%*%Rinv%*%Rx
  
  if ( approx ){
    Sigma.approx <- Rinv+Rinv%*%Rx%*%Rinv%*%Rx%*%Rinv
    obj.val <- sum(diag(Sigma.approx %*% M))
  } else {
    # Compute the original covariance matrix in Equation (7) of Zhang et al. (2021)
    Sigmainv <- try(solve(Sigma), TRUE) 
    # If matrix is invertible, then use it to find z which maximizes the variance
    # over the whole feasble covariate space. 
    if(is.matrix(Sigmainv)) {
      obj.val <- sum(diag(Sigmainv %*% M))
    }
  }
  return(obj.val)
}

Ioptimality_multiple <- function(H, X, m, n.cat, model.type, M = NULL) {
  #=============================================================================
  # Objective function.
  # If the number of covariates is smaller than 50,
  # the function outputs the maximum variance over the
  # whole feasible space of covariates.
  # Otherwise, it outputs the maximum variance 
  # over the covariate vectors in matrix Zs.
  #=============================================================================
  
  # Load useful objects.
  p <- ncol(H)
  n <- nrow(H)
  R <- t(H)%*%H # H'H
  Rinv <- solve(R) # (H'H)^{-1}
  if(is.null(M)){
    M <- moments.matrix(m, n.cat, type = model.type) # M matrix.  
  }
  
  ss <- ncol(X)
  ave.obj.approx <- apply(X, MARGIN = 2, FUN = sing_Iopt, p=p, n=n, ss=ss, H=H, R=R, 
                   Rinv=Rinv, M=M, approx = TRUE)
  ave.obj <- apply(X, MARGIN = 2, FUN = sing_Iopt, p=p, n=n, ss=ss, H=H, R=R, 
                          Rinv=Rinv, M=M, approx = FALSE)
  
  re <- data.frame(Iopt=ave.obj, Iopt.approx=ave.obj.approx)
  return(re)
}
