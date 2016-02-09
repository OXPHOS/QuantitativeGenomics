MAF <- function(loci){
  maf <-min(table(loci)[1]/sum(table(loci)),table(loci)[2]/sum(table(loci)))
  return(maf)
}

XaXd <- function(loci){
  nt2num <- NULL
  min <- names(which(table(loci)==min(table(loci))))
  nt2num[loci==min] <- 1; nt2num[loci!=min] <- 0
  
  Xa <- NULL
  for (i in 1:(length(loci)/2)){
    Xa <- rbind(Xa,nt2num[2*i-1]+nt2num[2*i]-1)
  }
  return(Xa)
}

estimation_MLE<- function(Xa,Xd,Y){
  X <- cbind(1,Xa,Xd)
  Beta <- solve(t(X)%*%X)%*%(t(X))%*%Y
  Yhat <- X%*%Beta
  return(rbind(Beta,Yhat))
}

calc_stats_pvals<- function(sample.size,Yhat,Y){
  SSM <- sum((Yhat-mean(Y))^2)
  SSE <- sum((Y-Yhat)^2)
  MSM <- SSM/2
  MSE <- SSE/(sample.size-3)
  Fstat <- MSM/MSE
  pvals <- pf(Fstat,2,sample.size-3,lower.tail=FALSE)
  return(pvals)
}

calc_corr <- function(a,b){
  Xa <- NULL; Xb <- NULL
  mina <- names(which(table(a)==min(table(a))))
  minb <- names(which(table(b)==min(table(b))))
  Xa[a==mina] =1;Xa[a!=mina]=0 
  Xb[b==minb] =1;Xb[b!=minb]=0
  r <- cov(Xa,Xb)/(sqrt(var(Xa))*sqrt(var(Xb))) 
  return(r)
}

#QUESTION 1&2
#Data input
phenotypes.read<- read.csv('midterm_phenotypes_fall13.txt',header=F)
phenotypes <- as.matrix(phenotypes.read) 
sample.size <- nrow(phenotypes)  #nrow = dim [1]
genotypes.read <- read.csv('midterm_genotypes_fall13.txt',header=F)
genotypes <- as.matrix(genotypes.read) 
SNP.size <- ncol(genotypes)

#Calculation of MAF
maf<- apply(genotypes,2,MAF)

#Plot output of Q1&2
png('Q1&2.png',height=1000,width=2000,res=300)
par(mfrow=c(1,2))
hist(phenotypes, main='',ylim=c(0,80))
hist(maf,main='',xlab='MAF',xlim=c(0,0.5),ylim=c(0,200))
dev.off()

#QUESTION 3&4
#Transform nucleotides into XaXd
Xa <- apply(genotypes,2,XaXd)
Xd <- 1-abs((Xa*2))

#MLE,pvalue calculation
beta.u <- NULL; beta.a<- NULL; beta.d <- NULL;
pvals.ls <- NULL;beta.ls <- NULL;
for(i in 1:SNP.size){
  beta_Yhat<- estimation_MLE(Xa[,i],Xd[,i],phenotypes)
  beta.u <- c(beta.u,beta_Yhat[1,1])
  beta.a <- c(beta.a,beta_Yhat[2,1])
  beta.d <- c(beta.d,beta_Yhat[3,1])
  pvals <- calc_stats_pvals(sample.size,beta_Yhat[-1:-3,1],phenotypes)
  pvals.ls <- c(pvals.ls,pvals) 
  beta.ls <- rbind(beta.ls,beta_Yhat[1:3,1])
}

#Plot output of beta MLE
png('Q3.png',height=1000,width=2000,res=300)
par(mfrow=c(1,3))
hist(beta.u, main='');hist(beta.a, main='');hist(beta.d, main='')
dev.off()

#Plot output of pvals.
png('Q4.png',height=1000,width=1000,res=300)
hist(pvals.ls, main='',xlab='p-values')
dev.off()

#Mahattan Plot.
png('Q5.png',height=1000,width=3000,res=300)
plot(1:1000,-log10(pvals.ls),'l',ylab='-log10(pvalues)')
dev.off()

#QUESTION 6
alphae2 <- sum(pvals.ls<0.01)
alphae5 <- sum(pvals.ls<0.00001)

#QUESTION 7
peak.position <- which(pvals.ls<0.00001)#which returns column names fit the criteria.
peak<- data.frame(position=peak.position,set=1)
for (i in 2:length(peak.position))
{
  if (peak.position[i] != peak.position[i-1]+1) {
    peak$set[i]=peak$set[i-1]+1
  }else{
    peak$set[i]=peak$set[i-1]
  }
}
num.set <- max(peak$set)

#QUESTION 8
#peak2 <- which(pvals.ls==min(pvals.ls[-which(pvals.ls==min(pvals.ls))]))
pvals.rank = rank(pvals.ls, ties.method='first')
peak2 <- which(pvals.rank==3)
prox.peak = which(peak$position == peak2)
if (prox.peak > 1 && 
      peak$set[prox.peak]==peak$set[prox.peak - 1]){
  corr.left <- calc_corr(genotypes[,peak2-1],genotypes[,peak2])
}
if (prox.peak < nrow(peak) && 
      peak$set[prox.peak]==peak$set[prox.peak + 1]){
  corr.right <- calc_corr(genotypes[,peak2],genotypes[,peak2+1])
}

#QUESTION 9
x <- c(-1:1)
png('Q9.png',height=3000,width=3000,res=300)
par(mfrow=c(1,2))
plot(Xa[,peak2],phenotypes,ylim=range(phenotypes))
par(new=TRUE)
plot(x,beta.ls[peak2,1]+x*beta.ls[peak2,2],'l',ylim=range(phenotypes),axes=FALSE,xlab='',ylab='')
plot(Xd[,peak2],phenotypes,ylim=range(phenotypes))
par(new=TRUE)
plot(x,beta.ls[peak2,1]+x*beta.ls[peak2,3],'l',ylim=range(phenotypes),axes=FALSE,xlab='',ylab='')
dev.off()


  
