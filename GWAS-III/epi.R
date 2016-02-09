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
