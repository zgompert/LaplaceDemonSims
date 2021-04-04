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
sf<-list.files(pattern="so_")
#sf_t<-sf[8]
sf<-sf[c(2,7,9,1,4,5,6,3)]

prec_bv<-mkplots(sfl=sf,out="BvSimsFinch.pdf",yl=c(-1.0,1.75),lab="Mean BV")
10^prec_bv
#          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
#50%  564.48719 3.568060 4.203582 4.633216 20.90380 5.604278 7.297287 6.343187
#25%   75.49717 3.328521 4.089894 4.477487 15.00014 5.261162 6.954480 6.118540
#75% 2175.73796 4.475068 4.972971 6.682791 87.38361 7.672006 7.615365 7.203490


## 6 snps
kk<-1
sf_all<-vector("list",4)
for(k in c(3,4,39,1193)){
	sf<-list.files(pattern=glob2rx(paste("sop*_",k,".txt",sep="")))
	sf_t<-sf[8]
	sf_all[[kk]]<-sf[c(4,6,3)]
    kk<-kk+1
	#mkplotsSnps(sfl=sf,sft=sf_t,out=paste("Snp",k,"SimsFinch.pdf",sep=""),yl=c(0,1),lab="Frequency")
}


tits<-c("Drift + unc. genetics","Unc. sel. + weather + gen.","Drift + unc. sel. + w. + gen","Predictability")

library(scales)

## change to Mean BV and -1,1.6 for bounds for bv
mkplotsSnps<-function(sfl=NA,out="out.pdf",tit=tits,yl=c(0,1),lab="Frequency"){
	cl<-1.6
	ca<-1.1
	cm<-1.45
	pdf(file=out,width=12,height=12)
	par(mfrow=c(4,4))
	par(mar=c(4.5,5.5,2.5,1.5))
	NN<-length(sfl)
	N<-length(sfl[[1]])
	vars<-vector("list",NN)
	e_vars<-vector("list",NN)
	for(jj in 1:NN){ ## loop over SNPs
	    vars[[jj]]<-matrix(NA,nrow=N,ncol=35)
    	s_dat<-vector("list",N)
    	for(j in 1:N){
    		s_dat[[j]]<-as.matrix(read.table(sfl[[jj]][j],header=FALSE))
    		plot(1:35,rep(0,35),type='n',xlab="Time",ylab=lab,
    		     cex.lab=cl,cex.axis=ca,ylim=yl)
    		for(i in 1:99){
    			lines(1:35,s_dat[[j]][i,-36],col=alpha("gray",.5),lwd=.6)
    		}
    		lines(1:35,s_dat[[j]][100,-36],lwd=1)
    		jx<-j + (jj-1) * (N+1)
    		title(main=paste("(",LETTERS[jx],") ",tit[j],sep=""),cex.main=cm)
    		vars[[jj]][j,]<-apply(s_dat[[j]][,-36],2,var,na.rm=TRUE)
		mtext(round(median(1/vars[[jj]][j,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
    	}
	
	    e_vars[[jj]]<-log10(apply(1/vars[[jj]],1,quantile,probs=c(.5,.25,.75)))
    	ub<-max(as.vector(e_vars[[jj]]))*1.05
    	x<-barplot(e_vars[[jj]][1,],names.arg=LETTERS[((jj-1)*(N +1)+1):((jj-1)*(N+1) +N)],xlab="Conditions",ylab="Precision (log10 scale)",cex.lab=cl,cex.axis=ca,
		cex.names=ca,ylim=c(0,ub))
        segments(x,e_vars[[jj]][2,],x,e_vars[[jj]][3,])
        jx<-(N+1) + (jj-1) * (N+1)
	    title(main=paste("(",LETTERS[jx],") ",tit[j+1],sep=""),cex.main=cm)
    }
	dev.off()
	return(e_vars)
}

prec_snps<-mkplotsSnps(sfl=sf_all,out="SnpSimsFinch.pdf",yl=c(0,1),lab="Frequency")
prec_snps
#[[1]]
#        [,1]      [,2]      [,3]
#50% 1.237148 0.7488244 0.7040656
#25% 0.741145 0.7266163 0.6902307
#75% 1.869948 0.7917364 0.7389513

#[[2]]
#        [,1]      [,2]      [,3]
#50% 1.230229 1.0507972 1.0366754
#25% 1.188322 0.9999504 0.9735922
#75% 1.490610 1.1136724 1.0802171

#[[3]]
#        [,1]     [,2]     [,3]
#50% 1.854983 1.320204 1.545931
#25% 1.153537 1.241413 1.502693
#75% 2.536480 1.352471 1.660692

#[[4]]
#        [,1]     [,2]     [,3]
#50% 1.543812 3.741257 1.651341
#25% 1.337271 3.587720 1.493152
#75% 1.868316 4.072093 1.832537


recip_sqrt_l10<-function(x){
    x<-10^x
    y<-1/x
    y<-sqrt(y)
    return(y)
    }
    
lapply(prec_snps,recip_sqrt_l10)
#[[1]]
#         [,1]      [,2]      [,3]
#50% 0.2406722 0.4222676 0.4445977
#25% 0.4260176 0.4332034 0.4517359
#75% 0.1161518 0.4019128 0.4270949

#[[2]]
#         [,1]      [,2]      [,3]
#50% 0.2425970 0.2982644 0.3031533
#25% 0.2545887 0.3162458 0.3259897
#75% 0.1797608 0.2774366 0.2883311

#[[3]]
#          [,1]      [,2]      [,3]
#50% 0.11817038 0.2187249 0.1686686
#25% 0.26499136 0.2394935 0.1772774
#75% 0.05392125 0.2107485 0.1477930

#[[4]]
#         [,1]        [,2]      [,3]
#50% 0.1690806 0.013470128 0.1493928
#25% 0.2144690 0.016074602 0.1792354
#75% 0.1163702 0.009203505 0.1212639

    
