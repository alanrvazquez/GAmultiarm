################################################################################
#
# OPTIMIZATION ALGORITHMS TO CONSTRUCT TWO-ARM CLINICAL TRIALS. 
# 
#
# Author: Alan Vazquez
# Affiliation: Tecnologico de Monterrey
# Email: alanrvazquez@tec.mx
#
################################################################################

library(gurobi)

CEAlgorithm <- function(H, m, n.cat = 0, max.iter = 5, tol = 0.0000001, model.type = 'l', M = NULL,
                         use.approx = TRUE, my.seed = 1818975, cnst = 1000){
  
  ##############################################################################
  # Coordinate-exchange algorithm in Section 5.  
  
  # INPUT: 
  #   H: n x p matrix with the n observations on the p predictors.
  #   n.cat: Number of categorical predictors in H. The categorical predictors
  #        should be at the end of the matrix.
  #   max.iter: Maximum number of iterations for the algorithm.
  #   tol: Tolerance.
  #   use.approx: TRUE if the approximated covariance matrix is used. 
  #             FALSE if the actual covariance matrix is used.
  #   my.seed: Set seed for reproducibility.
  #   cnst: Constant for numerical stability during the optimization. The
  #         values of the objective function are small.
  
  # OUTPUT: 
  #   Y: Best design according to the approximated average SC variance criterion 
  #      among the "max.iter" iterations.
  
  ##############################################################################
  
  # Allocate preliminary objects.
  p <- ncol(H)
  n <- nrow(H)
  R <- t(H)%*%H # H'H
  Rinv <- solve(R) # (H'H)^{-1}
  if(is.null(M)){
    M <- moments.matrix(m, n.cat, type = model.type) # M matrix.  
  }
  
  
  # Create weights for objective function.
  Cmat <- matrix(NA, ncol = n, nrow = n)
  AuxMat.Three <- Rinv%*%M
  for (u in 1:n){
    for (v in 1:n){
      Umat <- H[u,]%*%t(H[u,])
      Vmat <- H[v,]%*%t(H[v,])
      AuxMat.One <- Rinv%*%Umat # J_i in Proposition 2.
      AuxMat.Two <- Rinv%*%Vmat # J_j in Proposition 2.
      Cmat[u, v] <- sum( diag(AuxMat.One%*%AuxMat.Two%*%AuxMat.Three) )
    }  
  }
  
  # Create starting solutions.
  set.seed(my.seed)
  Y <- sample(c(-1,1), size = n*max.iter, replace = TRUE)
  Y <- matrix(Y, ncol = max.iter, nrow = n) 
  
  # Apply heuristic algorithm.
  npass <- n*20
  for (j in 1:max.iter){
    
    citer <- 1
    Copy.Y <- Y[,j] # Copy initial solution.
    obj.pass <- 10*3
    
    while (citer <= npass){
      
      for (i in 1:n){
        # Apply Equation (10).
        eval.cond <- cnst*Y[-i, j]%*%Cmat[i, -i] 
        
        if( eval.cond > tol){
          Y[i, j] <- -1
          } else {
            Y[i,j] <- 1
            }
        }
      
      # Check if we need to go over all coordinates again.
      break.cond <- sum(abs(Copy.Y - Y[,j])) # Has the solution changed?
      
      if (break.cond  > tol ){ 
        # If it changed, Copy new improved solution.
        Copy.Y <- Y[,j] 
      
        } else {
        # Otherwise, exit the while loop.
        break 
          }
      citer <- citer + 1
      
    } # end while
    
  } # end for (j in 1:max.iter)
  
  # Evaluate optimized designs.
  Iopt.values <- Ioptimality(H, Y, m, n.cat, model.type, use.approx)
  
  # Select the best design.
  lab.best <- which.min(Iopt.values)
  return(as.matrix(Y[,lab.best]))
}

IPAlgorithm <- function(H, m, n.cat = 0, e.cnst = 6, model.type = 'l', M = NULL, 
                        max.time = 60, print.output = 1, tol = 0.001){
  
  ##############################################################################
  # Integer programming algorithm in Section 4.  
  
  # INPUT: 
  #   H: n x p matrix with the n observations on the p predictors.
  #   n.cat: Number of categorical predictors in H. The categorical predictors
  #        should be at the end of the matrix.
  #   max.time: Maximum number of seconds for solving the integer programming 
  #             problem using Gurobi. 
  #   print.output: 1 if yes, 0 otherwise.
  #   tol: Tolerance for numerical stability of Gurobi.
  
  # OUTPUT: 
  #   X: Optimal design according to the approximated I-optimality criterion.
  
  ##############################################################################
  
  # Allocate preliminary objects.
  p <- ncol(H)
  n <- nrow(H)
  R <- t(H)%*%H # H'H
  Rinv <- solve(R) # (H'H)^{-1}
  n.choose.two <- choose(n,2)
  ncombinat.two <- combn(1:n, 2)
  if(is.null(M)){
    M <- moments.matrix(m, n.cat, type = model.type) # M matrix.  
  }
  
  # Create matrix of weights.
  Cmat <- matrix(NA, ncol = n, nrow = n)
  AuxMat.Three <- Rinv%*%M
  for (u in 1:n){
    for (v in 1:n){
      Umat <- H[u,]%*%t(H[u,])
      Vmat <- H[v,]%*%t(H[v,])
      AuxMat.One <- Rinv%*%Umat # J_i in Proposition 2.
      AuxMat.Two <- Rinv%*%Vmat # J_j in Proposition 2.
      Cmat[u, v] <- sum( diag(AuxMat.One%*%AuxMat.Two%*%AuxMat.Three) )
    }  
  }
  
  # Matrices for constraints 1 and 2.
  Mat.const.one <- matrix(0, ncol = n, nrow = n.choose.two)
  Mat.const.two <- matrix(0, ncol = n, nrow = n.choose.two)
  for (i in 1:n.choose.two){
    Mat.const.one[i, ncombinat.two[1,i]] <- -1
    Mat.const.two[i, ncombinat.two[2,i]] <- -1
  }
  A.mat.one <- cbind(diag(n.choose.two), Mat.const.one)
  A.mat.two <- cbind(diag(n.choose.two), Mat.const.two)
  
  # Matrix for Constraint 3.
  Mat.const.three <- matrix(0, ncol = n, nrow = n.choose.two)
  for (i in 1:n.choose.two){
    Mat.const.three[i, ncombinat.two[,i]] <- 1
  }
  A.mat.three <- cbind(-1*diag(n.choose.two), Mat.const.three)
  
  # Matrix for Constraints 4 and 5.
  A.mat.four <- c(rep(0, times = n.choose.two), rep(1, times = n))
  
  # Coefficients for constraints.
  b.vec <- c(rep(0, times = 2*n.choose.two), rep(1, times = n.choose.two))
  
  if (n %% 2 == 0){
    b.vec <- c(b.vec, e.cnst - n/2, n/2 + e.cnst)  # If multiple of 2.
  } else {
    b.vec <- c(b.vec, e.cnst - (n-1)/2, (n-1)/2 + e.cnst)  # If not multiple of 2.
  }
  
  
  const.directions <- c(rep("<=", times = 3*n.choose.two + 2))
  
  # Coefficients for objective function.
  c.vec.obj <- rep(0, times = n.choose.two + n)
  for (i in 1:n.choose.two){
    c.vec.obj[i] <- 4*Cmat[ncombinat.two[1,i], ncombinat.two[2,i]]
    c.vec.obj[n.choose.two + ncombinat.two[1,i]] <- c.vec.obj[n.choose.two + ncombinat.two[1,i]] - 2*Cmat[ncombinat.two[1,i], ncombinat.two[2,i]]
    c.vec.obj[n.choose.two + ncombinat.two[2,i]] <- c.vec.obj[n.choose.two + ncombinat.two[2,i]] - 2*Cmat[ncombinat.two[1,i], ncombinat.two[2,i]]
  }
  
  # Define Gurobi Model.
  model <- list()
  model$obj <- 100*c.vec.obj
  model$A <- rbind(A.mat.one, A.mat.two, A.mat.three, -1*A.mat.four, A.mat.four)
  model$modelsense <- 'min'
  model$rhs <- b.vec
  model$sense <- const.directions
  model$vtype <- 'B'
  
  params <- list(MIPGap = tol, TimeLimit = max.time, OutputFlag = print.output)
  
  # Solve Gurobi model.
  result <- gurobi(model, params)
  
  # Transform to binary variables to original variables with -1 and +1.
  u.var.sol <- result$x[(n.choose.two+1):(n.choose.two+n)]
  X.BL.sol <- as.matrix(2*u.var.sol - 1)
  
  return(X.BL.sol)
}
