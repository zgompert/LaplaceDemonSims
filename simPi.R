## tryiing to match params from Chaves et al. 2016, doi:10.1111/mec.13743
## 15 SNPs with pip > 0.1
## ng = 28.3 (16, 56)
## pve = 0.947 (0.899, 0.979)
## pge = 0.775 (0.638, 0.887)

pi<-c(runif(15,0.1,0.6),runif(100,0.05,0.1),runif(300,0.01,0.05),runif(1000,0.001,0.01))
ng<-rep(NA,1000)
for(i in 1:1000){
	ng[i]<-sum(rbinom(n=length(pi),size=1,prob=pi))
}
## this would probably do it
quantile(ng,probs=c(.5,0.025,.975))

p0<-rbeta(length(pi),.6,.6)

pitrue<-rbinom(n=length(pi),size=1,prob=pi)

beta<-rnorm(length(pi),0,1)/sum(pi)

write.table(round(cbind(pi,beta),5),file="genarch_unc.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(round(cbind(pitrue,beta),5),file="genarch_known.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

## turn inital allele freq. into a matrix
p0mat<-matrix(rep(p0,40),nrow=length(p0),ncol=40,byrow=FALSE)

write.table(round(p0mat,5),file="p0.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

## Ne, from https://esajournals.onlinelibrary.wiley.com/doi/abs/10.2307/1940156
## turn inital allele freq. into a matrix
nemat<-matrix(60,nrow=length(p0),ncol=40,byrow=FALSE)
write.table(nemat,file="ne.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)


save(list=ls(),file="gendat.rdat")
