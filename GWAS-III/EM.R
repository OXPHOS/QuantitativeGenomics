A <- cov.wt(X)$cov

EM_algorithm = function(Y, X_j, A){
  # Calculate the inverse of A once since it is used repeatedly in the algorithm
  # This method is faster than solve(A)
  solve_A = chol2inv(chol(A))
  n = length(Y)
  I = diag(1, n)
  log_L = c()
  # set starting values
  sigma_sq_a = 5
  sigma_sq_e = 5
  beta = as.vector(rep(0, ncol(X_j)))
  C = A * sigma_sq_a + I * sigma_sq_e
  log_L[1] = -1/2 * determinant(C)$modulus - 1/2 * t(Y - X_j %*% beta) %*% chol2inv(chol(C)) %*% (Y - X_j %*% beta)
  iter = 2
  while(1){
    S = chol2inv(chol(I + solve_A * sigma_sq_e / sigma_sq_a))
    alpha = S %*% (Y - X_j %*% beta)
    
    V = S  * sigma_sq_e
    
    beta = chol2inv(chol(t(X_j) %*% X_j)) %*% t(X_j) %*% (Y - alpha)
    
    # use as.numeric() to make sure value is saved as a scalar
    sigma_sq_a = as.numeric(1/n * (t(alpha) %*% solve_A %*% alpha + sum(diag(solve_A %*% V))))
    
    # use as.numeric() to make sure value is saved as a scalar
    sigma_sq_e = as.numeric( 1/n * ( t(Y - X_j %*% beta - alpha) %*% (Y - X_j %*% beta - alpha) + sum(diag(V))))
    
    C = A * sigma_sq_a + I * sigma_sq_e
    log_L[iter] = -1/2 * determinant(C)$modulus - 1/2 * t(Y - X_j %*% beta) %*% chol2inv(chol(C)) %*% (Y - X_j %*% beta)
    
    if(log_L[iter] - log_L[iter-1] < 1e-5){
      break;
    }
    iter = iter + 1
  }
  
  return(list(beta = beta, 
              sigma_sq_a = sigma_sq_a, 
              sigma_sq_e = sigma_sq_e, 
              log_L = log_L[iter-1]))  
}


# Null model
log_L_null = EM_algorithm(Y, Xmu, A)$log_L
# log_L_null = Mclust(cbind(Y, One))$loglik

pvals.EM = c()
# Full model
for(j in 1:sample.size){
  X_j = cbind(1, X[j, ])
  fit = EM_algorithm(Y, X_j, A)
  pvals.EM[j] = pchisq(-2 * (log_L_null - fit$log_L), 1, lower.tail=FALSE)
  cat('j=',j,'\n')
}

SNP.EM <- which(pvals.EM<1e-4)
SNP.EM.sites <- data.frame(names(predata[SNP.EM+6]),pvals.EM[SNP.EM])
names(SNP.EM.sites) <- c('genetic markers','p-values')
save(SNP.EM.sites,file = 'candidate casual SNP sites with mixed regression model.Rdata')

png("Manhattan plot for GWAS with mixed model.png",width=800,height=500)
plot(-log10(pvals.EM[-201:-231]), main="Fig.5 Manhattan plot under mixed model", xlab="genetic markers in analysis", ylab="-log10 p-values",type='p')
SNPabline(4,0)
dev.off()

obs.mlogpval.vec.EM <- sort( -log( pvals.EM ) )
exp.mlogpval.vec.EM <- sort( -log( c(1:length(obs.mlogpval.vec.EM)) * (1/length(obs.mlogpval.vec.EM)) ) )
png("QQ plot for GWAS with mixed model.png", width=500, height=500)
plot(exp.mlogpval.vec.EM, obs.mlogpval.vec.EM, main="Fig.6 QQ plot under mixed model", xlab="expected -log pvalue assuming null hypothesis true for every test", ylab="observed -log p-values",xaxs='i',yaxs='i')
abline(0,1)
dev.off()
