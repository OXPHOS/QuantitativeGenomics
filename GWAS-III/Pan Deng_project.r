
##############################################################
########-------PRELIMINARY DATA FILTER IN PLINK-------########
##############################################################
# codes for command line
# plink
# --noweb
# --file 3L_thinM
# --pheno phenoCG9186_3L.txt
# --maf 0.05
# --geno 0 
# --mind 0.1
# --prune
# --recodeA
# --out PreDataFilter

##############################################################
#########------------PHENOTYPE CHECK------------#############
##############################################################
pheno <- read.csv('phenoCG9186_3L.txt',sep='')
png("Distribution of CG9186 Expression.png", width=1000, height=500)
hist(as.matrix(pheno[3]), main='', xlab='Distribution of CG9186 Expression Level in all 94 lines',xlim=c(-40,80),nclass=40)
dev.off()

##############################################################
############------------DATA INPUT------------################
##############################################################
predata <- read.csv('PreDataFilter.raw',sep='')
genotypes <- predata[-1:-6]
sample.size <- dim(genotypes)[2]
line.num <- dim(genotypes)[1]
Xa <- as.matrix(genotypes-1)
Xd <- as.matrix(1-abs((genotypes-1)*2))
Y  <- as.matrix(predata[6])
Xmu  <- matrix(1,dim(Y)[1],1)
#Imporve to print out the genotype with image(Xa) here#

##############################################################
############------------PCA ANALYSIS------------##############
##############################################################
X <- t(genotypes)
W<- (X - rowMeans(X)) / sqrt(diag((cov(X))))
#perform a PCA
geno.pc <- princomp(W)
#plot the loadings of the first two PCs
png("Plot of n individuals on the first two PCs.png", width=500, height=500)
plot(geno.pc$loadings[,c(1,2)], main="Plot of n individuals on the loadings of the first two PCs", xlab="PC1", ylab="PC2")
dev.off()

##############################################################
##########---------LINEAR MODEL REGRESSION---------###########
##############################################################
Xz <- rep(1,line.num)
Xz[which(geno.pc$loadings[,1]>0.5)]=-1
pvals.LM <- NULL
# MLE and p value for linear model#
for (i in 1:sample.size){
  LM <- data.frame(Y,Xmu,Xa[,i],Xd[,i],Xz)
  names(LM) <- c('Y','Xmu','Xa','Xd','Xz')
  L0 <- lm(Y ~ Xmu, data=LM )
  L1 <- lm(Y ~ Xmu + Xa + Xd + Xz, data=LM )
  out.aov <- anova(L0, L1)
  pvals.LM <- c(pvals.LM, out.aov$"Pr(>F)"[2])
}
SNP.LM <- which(pvals.LM<1e-4)
SNP.LM.sites <- data.frame(names(predata[SNP.LM+6]),pvals.LM[SNP.LM])
names(SNP.LM.sites) <- c('genetic markers','p-values')
save(SNP.LM.sites,file = 'candidate casual SNP sites with linear regression model.Rdata')

#Mahattan plot and QQ plot#
obs.mlogpval.vec.LM <- sort( -log( pvals.LM ) )
exp.mlogpval.vec.LM <- sort( -log( c(1:length(obs.mlogpval.vec.LM)) * (1/length(obs.mlogpval.vec.LM)) ) )
png("Mahattan plot and QQ plot for GWAS with linear model", width=500, height=800)
par(mfrow=c(2,1))
plot(-log10(pvals.LM),type='l',main='Mahattan Plot under linear model', xlab='genetic markers in analysis', ylab='-log p-values')
abline(4,0)
plot(exp.mlogpval.vec.LM, obs.mlogpval.vec.LM, main="QQ plot under linear model", xlab="expected -log pvalue assuming null hypothesis true for every test", ylab="observed -log p-values",xaxs='i',yaxs='i')
abline(0,1)
dev.off()

##############################################################
##########-------LINEAR MODEL WITH EPISTASIS-------###########
##############################################################
#MLE of epistasis and F-stat#
epistasis <- function(Y,Xa1,Xd1,Xa2,Xd2){
  Xaa <- Xa1*Xa2
  Xa1d2 <- Xa1*Xd2
  Xa2d1 <- Xa2*Xd1
  Xdd <- Xd1*Xd2
  epi.df <- data.frame(Y,Xmu,Xa1,Xd1,Xa2,Xd2,Xaa,Xa1d2,Xa2d1,Xdd,Xz)
  names(epi.df) <- c('Y','Xmu','Xa1','Xd1','Xa2','Xd2','Xaa','Xa1d2','Xa2d1','Xdd','Xz')
  #Consider only epistasis effect here
  L0 <- lm(Y~Xmu+Xa1+Xd1+Xa2+Xd2+Xz,data=epi.df)
  L1 <- lm(Y~Xmu+Xa1+Xd1+Xa2+Xd2+Xaa+Xa1d2+Xa2d1+Xdd+Xz,data=epi.df)
  out.aov <- anova(L0,L1)
  return(out.aov$'Pr(>F)'[2])
}
# main function, run whole matrix#
pvals.epi <- matrix(1:1,sample.size,sample.size)
for (i in 1:dim(Xa)[2]){
  for(j in 1:dim(Xa)[2]){
    if (j<i) {
      pvals.epi[i,j] <- epistasis(Y,Xa[,i],Xd[,i],Xa[,j],Xd[,j])
    }
  }
  #cat('i=',i,'\n')
}


SNP.epi <- which(pvals.epi<1e-8,t)
SNP.epi.1 <- NULL;SNP.epi.2 <- NULL;
for (t in 1:nrow(SNP.epi)){SNP.epi.1 <- c(SNP.epi.1,names(predata[SNP.epi[t,1]+6]));SNP.epi.2 <- c(SNP.epi.2,names(predata[SNP.epi[t,2]+6]))}
SNP.epi.sites <- data.frame(SNP.epi.1,SNP.epi.2,pvals.epi[which(pvals.epi<1e-8)])
names(SNP.epi.sites) <- c('epi SNP site 1','epi SNP site 2','p-values')
save(SNP.epi.sites,file='candidate casual SNP sites with epistasis analysis.Rdata')

png('p-values of pairwise gene epistasis analysis.png',width=1000,height=1000)
image(pvals.epi,xlab='p-values of pairwise gene epistasis analysis',xaxt='n', yaxt='n')
dev.off()

obs.mlogpval.vec.epi <- sort( -log( pvals.epi ) )
exp.mlogpval.vec.epi <- sort( -log( c(1:length(obs.mlogpval.vec.epi)) * (1/length(obs.mlogpval.vec.epi)) ) )
png('QQ plot of epistasis analysis.png',width=1000,height=1000)
plot(exp.mlogpval.vec.epi, obs.mlogpval.vec.epi, main="QQ plot of epistasis analysis", xlab="expected -log pvalue assuming null hypothesis true for every test", ylab="observed -log p-values",xaxs='i',yaxs='i')
abline(0,1)
dev.off()

##############################################################
##########---------MIXED REGRESSION MODEL----------###########
##############################################################
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
  cat('.')
}

SNP.EM <- which(pvals.EM<1e-4)
SNP.EM.sites <- data.frame(names(predata[SNP.EM+6]),pvals.EM[SNP.EM])
names(SNP.EM.sites) <- c('genetic markers','p-values')
save(SNP.EM.sites,file = 'candidate casual SNP sites with mixed regression model.Rdata')

png("Manhattan plot for GWAS with mixed model.png",width=800,height=500)
plot(-log10(pvals.EM), main="Manhattan plot under mixed model", xlab="genetic markers in analysis", ylab="-log10 p-values")
abline(4,0)
dev.off()

obs.mlogpval.vec.EM <- sort( -log( pvals.EM ) )
exp.mlogpval.vec.EM <- sort( -log( c(1:length(obs.mlogpval.vec.EM)) * (1/length(obs.mlogpval.vec.EM)) ) )
png("QQ plot for GWAS with mixed model", width=500, height=500)
plot(exp.mlogpval.vec.EM, obs.mlogpval.vec.EM, main="QQ plot under mixed model", xlab="expected -log pvalue assuming null hypothesis true for every test", ylab="observed -log p-values",xaxs='i',yaxs='i')
abline(0,1)
dev.off()

