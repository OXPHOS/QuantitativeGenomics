##############################################################
#########------------PHENOTYPE CHECK------------#############
##############################################################
pheno <- read.csv('phenoCG9186_3L.txt',sep='')
png("Distribution of CG9186 Expression.png", width=1000, height=500)
hist(as.matrix(pheno[3]), main='Fig.1 Distribution of CG9186 Expression Level in all 94 lines', xlab='',xlim=c(-40,80),nclass=40)
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
plot(geno.pc$loadings[,c(1,2)], main="Fig.2 Plot of n individuals on the loadings of the first two PCs", xlab="PC1", ylab="PC2")
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
png("Manhattan plot for GWAS with linear model.png",width=800,height=500)
plot(-log10(pvals.LM), main="Fig.3 Manhattan plot under linear model", xlab="genetic markers in analysis", ylab="-log10 p-values",type='p')
abline(4,0)
dev.off()

obs.mlogpval.vec.LM <- sort( -log( pvals.LM ) )
exp.mlogpval.vec.LM <- sort( -log( c(1:length(obs.mlogpval.vec.LM)) * (1/length(obs.mlogpval.vec.LM)) ) )
png("QQ plot for GWAS with linear model.png", width=500, height=500)
plot(exp.mlogpval.vec.LM, obs.mlogpval.vec.LM, main="Fig.4 QQ plot under linear model", xlab="expected -log pvalue assuming null hypothesis true for every test", ylab="observed -log p-values",xaxs='i')
abline(0,1)
dev.off()

##############################################################
##########-------LINEAR MODEL WITH EPISTASIS-------###########
##############################################################
