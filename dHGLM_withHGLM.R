
#Cholesky decomp of GRM
L <- chol(Lasky$G)
Z0 <- diag(1, nrow = nrow(Lasky$G), ncol = ncol(Lasky$G))
Z <- Z0 %*% L

Y <- Lasky$Pheno$PC1 #Phenotype
std.Y <- transform2Normal(matrix(Y, nrow = length(Y), ncol = 1))

n <- nrow(Lasky$Zmat) #Number of individuals
w.mean <- rep(1,n) #Initial weights
res.var = 0 #Initial values
iter = 0
MAXiter = 20
while (abs(res.var-1) > 0.0001 & iter < MAXiter) {
  iter = iter +1
  if(iter == 1){
    hg1 <- hglm(y = std.Y, X = matrix(1, n, 1), Z = Lasky$Zmat, weights = w.mean)
  }else{
    hg1 <- hglm(y = std.Y, X = matrix(1, n, 1), Z = Lasky$Zmat, weights = w.mean, startval = startVals)
  }
  res.var <- hg1$varFix #new residual variance
  cat("Iteration:", iter,"\n", "Convergence to 1:", res.var, "\n")
  res <- Y - hg1$fv[1:n] #residuals
  hv <- hg1$hv[1:n] #hat values for the fixed part of the model?
  dev <- (res^2)/(1-hv) #response for dispersion model
  hg2 <- hglm(y = dev, X = matrix(1,n,1), Z = Lasky$Zmat,
              family = Gamma(link=log), weights = (1-hv)/2, verbose = T, maxit = 50)
  w.mean <- 1/hg2$fv #New weights
  
  startVals <- c(hg1$fixef, hg1$ranef, hg1$varRanef, hg1$varFix)
}

saveRDS(list(mean = hg1, disp = hg2), "~/Documents/Dropbox/Work/Sorghum_vQTL/dHGLM_PC1.std.rds")
