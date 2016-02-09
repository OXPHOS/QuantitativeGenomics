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

logistic.IRLS <- function(y, x, d.stop.th = 1e-6, it.max = 100) {
  
  #initialize the beta parameter vector at t=0
  beta.t = rep(0,ncol(x))
  
  for(i in 1:it.max) {
    t<-x%*%beta.t
    #cat (beta.t, "\n")
    gamma_inv0<-exp(t)/(1+exp(t))
    
    #initialize dt
    dt<-2*( sum(y[y==1]*log(y[y==1]/gamma_inv0[y==1])) + sum((1-y[y==0])*log((1-y[y==0])/(1-gamma_inv0[y==0]))) )
    #cat (dt, "\n")
    
    W<-matrix(0,dim(x)[1],dim(x)[1])
    for(k in 1:dim(x)[1]){
      W[k,k] <- gamma_inv0[k]%*%(1-gamma_inv0[k])
    }
    beta.t <- beta.t + solve(t(x)%*%W%*%x)%*%t(x)%*%(y-gamma_inv0)
    #cat (beta.t, "\n")
    #update gamma since it's a function of beta
    t<-x%*%beta.t
    gamma_inv1<-exp(t)/(1+exp(t))
    #cat (gamma_inv, "\n")
    
    #calculate new divergence
    dpt1 <- 2*( sum(y[y==1]*log(y[y==1]/gamma_inv1[y==1])) + sum((1-y[y==0])*log((1-y[y==0])/(1-gamma_inv1[y==0]))) )
    
    #cat(dpt1, "\n")
    absD <- abs(dt - dpt1)
    #cat (absD, "\n")
    
    
    if(absD < d.stop.th) {
      #cat("Convergence at iteration:", i, "at threshold:", d.stop.th, "\n")
      logl<-sum(y*log(gamma_inv1)+(1-y)*log(1-gamma_inv1))
      return(logl)
    }	
  }
  
  return("Threshold convergence not reached")
}

# QUESTION 1 discription of phenotypes
phenotypes<- as.matrix(read.csv('final_phenotypes_fall13.txt',header=F))
sample.size <- nrow(phenotypes) 
png('Phenotype Distribution.png',height=400,width=600)
hist(phenotypes,main='Phenotype Distribution',xlab='')
dev.off()

# QUESTION 2 calculate MAF and filter genotypes
MAF <- function(loci){
  maf <-min(table(loci)[1]/sum(table(loci)),table(loci)[2]/sum(table(loci)))
  return(maf)
}

genotypes <- as.matrix(read.csv('final_genotypes_fall13.txt',header=F))
maf<- apply(genotypes,2,MAF)
genotypes.filtered <- genotypes[,which(maf>0.05)]
SNP.size <- ncol(genotypes.filtered)

#QUESTION 3. parameters estimation and p-values calculation
Xa <- apply(genotypes.filtered,2,XaXd)
Xd <- 1-abs((Xa*2))
y <- phenotypes

pvals.ls<-NULL
logl_HA<-NULL
for(i in 1:SNP.size){
  
  x <- cbind(1, Xa[,i], Xd[,i])
  logl<-logistic.IRLS(y, x)
  logl_HA=c(logl_HA,logl)
}

x <- matrix(1,sample.size,1)
logl_H0<-logistic.IRLS(y, x)

LRT<-2*logl_HA-2*logl_H0 #likelihood ratio test statistic
pval <- pchisq(LRT, 2, lower.tail = F)
pvals.ls<-c(pvals.ls, pval)

#Plot p-values of alternative hypothesis
png('p-value distribution.png',height=300,width=600)
hist(pvals.ls,main='p-value distribution',xlab='p-values')
dev.off()

#QUESTION 4. Mahattan plot and significant SNP sites.
threshold <- 0.05/SNP.size
SNP.sig <- which(pvals.ls<threshold)

png('Mahattan plot.png',height=300,width=600)
plot(-log10(pvals.ls),type = 'l',ylab='-log10 p-values',main='Mahattan plot', xlab='SNP sites')
abline( -log10(threshold),0)
dev.off()

#QUESTION 5 QQ plot
obs.mlogpval.vec <- sort( -log(pvals.ls ) )
exp.mlogpval.vec <- sort( -log( c(1:length(pvals.ls)) * (1/length(pvals.ls)) ) )
png("QQ plot for GWAS with linear model.png", width=500, height=500)
plot(exp.mlogpval.vec, obs.mlogpval.vec, main="QQ plot under linear model", xlab="expected -log pvalue assuming null hypothesis true for every test", ylab="observed -log p-values",xaxs='i',yaxs='i')
abline(0,1)
dev.off()

#QUESTION 6 Calculate correlation between phenotype and the most significant SNP site
pvals.rank = rank(pvals.ls, ties.method='first')
peak <- which(pvals.rank==1)
corXaY <- cor(Xa[,peak],y)
corXdY <- cor(Xd[,peak],y)

#QUESTION 7 & 8 epistasis analysis 
# Set different SNP clusters
peak.df<- data.frame(position=SNP.sig,set=1)
for (i in 2:length(SNP.sig))
{
  if (SNP.sig[i] != SNP.sig[i-1]+1) {
    peak.df$set[i]=peak.df$set[i-1]+1
  }else{
    peak.df$set[i]=peak.df$set[i-1]
  }
}
set.size <- max(peak.df$set)

#Find peaks in each cluster
peaks <- NULL
for (i in 1:set.size){
  peak.cluster <- peak.df$position[peak.df$set==i]
  pvals.rank = rank(pvals.ls[peak.cluster], ties.method='first')
  peaks <- c(peaks,peak.cluster[which(pvals.rank==1)])
}

#Epistasis analysis 1
pvals.epi1 <- matrix(0,set.size,set.size)
x <- matrix(1,sample.size,1)
logl_H0<-logistic.IRLS(y, x)
for (i in 1:(set.size - 1)){
  for (j in (i + 1):set.size){
    xa1 <- Xa[,peaks[i]]; xa2 <- Xa[,peaks[j]]; xd1 <- Xd[,peaks[i]]; xd2 <- Xd[,peaks[j]] 
    x <- cbind(1,xa1,xa2,xd1,xd2)
    logl<-logistic.IRLS(y, x)  
    LRT <- 2*logl-2*logl_H0 #likelihood ratio test statistic
    pvals.epi1[i,j] <- pchisq(LRT, 4, lower.tail = F)
  }
}

#Epistasis analysis 2
pvals.epi2 <- matrix(0,set.size,set.size)
x <- matrix(1,sample.size,1)
logl_H0<-logistic.IRLS(y, x)
for (i in 1:(set.size - 1)){
  for (j in (i + 1):set.size){
    xaa <- Xa[,peaks[i]]*Xa[,peaks[j]]; xdd <- Xd[,peaks[i]]*Xd[,peaks[j]] 
    xa1d2 <- Xa[,peaks[i]]*Xd[,peaks[j]]; xd1a2 <- Xa[,peaks[j]]*Xd[,peaks[i]]
    x <- cbind(1,xaa,xa1d2,xd1a2,xdd)
    logl<-logistic.IRLS(y, x)  
    LRT <- 2*logl-2*logl_H0 #likelihood ratio test statistic
    pvals.epi2[i,j] <- pchisq(LRT, 4, lower.tail = F)   
  }
}

#QUESTION 10
dataset<-as.matrix(read.csv("MetropolisHastings_output_fall13.txt",sep='',header=F))
png('MHalgorithm.png',width=600,height=300)
par(mfrow=c(1,2))
hist(dataset[,1],main='Bayesian analysis of u1',xlab='')
hist(dataset[,2],main='Bayesian analysis of u2',xlab='')
dev.off()
u1 <- mean(dataset[,1])
u2 <- mean(dataset[,2])

