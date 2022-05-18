## create summary figures for finch based simulations

tits<-c("(A) Drift","(B) Unc. selection","(C) Unc. sel. + weather","(D) Drift + unc. sel. + weather",
	"(E) Drift + unc. genetics","(F) Unc. sel. + genetics","(G) Unc. sel. + weather + gen.",
	"(H) Drift + unc. sel. + w. + gen","(I) Predictability")

library(scales)

## change to Mean BV and -1,1.6 for bounds for bv
mkplots<-function(sfl=NA,out="out.pdf",tit=tits,yl=c(0,1),lab="Frequency"){
	cl<-1.6
	ca<-1.1
	cm<-1.5
	pdf(file=out,width=10,height=10)
	par(mfrow=c(3,3))
	par(mar=c(4.5,5.5,2.5,1.5))
	N<-length(sfl)
	vars<-matrix(NA,nrow=N,ncol=35)
	s_dat<-vector("list",N)
	for(j in 1:N){
		s_dat[[j]]<-as.matrix(read.table(sfl[j],header=FALSE))
		plot(1:35,rep(0,35),type='n',xlab="Time",ylab=lab,
		     cex.lab=cl,cex.axis=ca,ylim=yl)
		for(i in 1:99){
			lines(1:35,s_dat[[j]][i,-36],col=alpha("gray",.5),lwd=.6)
		}
		lines(1:35,s_dat[[j]][100,-36],lwd=1)
		title(main=tit[j],cex.main=cm)
		vars[j,]<-apply(s_dat[[j]][,-36],2,var,na.rm=TRUE)
		mtext(round(median(1/vars[j,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
	}
	
	e_vars<-log10(apply(1/vars,1,quantile,probs=c(.5,.25,.75)))
	ub<-max(as.vector(e_vars))*1.05
	x<-barplot(e_vars[1,],names.arg=LETTERS[1:N],xlab="Conditions",ylab="Precision (log10 scale)",cex.lab=cl,cex.axis=ca,
		cex.names=ca,ylim=c(0,ub))
    segments(x,e_vars[2,],x,e_vars[3,])
	title(main=tits[N+1],cex.main=cm)

	dev.off()
	return(e_vars)
}

## bv
sf<-list.files(pattern="soa_")
sf<-sf[c(2,7,9,1,4,5,6,3)]

prec_bv<-mkplots(sfl=sf,out="BvSimsFinchAugmented.pdf",yl=c(-1.0,1.75),lab="Mean BV")
10^prec_bv

#50% 506.7248 188.1906 47.31550 44.49076 13.68209 11.85538 10.879791 12.81528
#25% 441.1097 170.7114 38.42154 39.34778 13.53634 11.52984  9.973983 12.63387
#75% 530.3767 218.6910 69.76206 46.65088 14.77050 12.90895 11.490963 13.05900


